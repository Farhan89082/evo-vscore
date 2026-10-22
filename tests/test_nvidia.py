import base64
import io
import json
import urllib.error
import zipfile

import numpy as np
import pytest

from evo_vscore.backends import BackendError, NvidiaBackend, get_api_key
from evo_vscore.backends.nvidia import KEY_ENV_VARS


def make_reply(logits, name="unembed.output.npy"):
    """Build a reply in the format the real API returned: JSON -> base64 zip -> npy."""
    buf = io.BytesIO()
    np.save(buf, np.asarray(logits)[None, :, :])  # shape (1, length, vocab)
    zbuf = io.BytesIO()
    with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(name, buf.getvalue())
    return json.dumps({"data": base64.b64encode(zbuf.getvalue()).decode()}).encode()


def http_error(code, headers=None, body=b"nope"):
    import email.message

    msg = email.message.Message()
    for k, v in (headers or {}).items():
        msg[k] = v
    return urllib.error.HTTPError("http://x", code, "err", msg, io.BytesIO(body))


class Fake(NvidiaBackend):
    """Backend whose network layer is replaced by a scripted list of outcomes."""

    def __init__(self, outcomes, **kw):
        self.sleeps = []
        self.now = 0.0

        def fake_sleep(s):
            self.sleeps.append(s)
            self.now += s

        super().__init__("nvapi-SECRET", sleep=fake_sleep, clock=lambda: self.now, **kw)
        self.outcomes = list(outcomes)
        self.bodies = []

    def _post(self, body):
        self.bodies.append(body)
        outcome = self.outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


def test_decode_and_score_roundtrip():
    seq = "ACGT"
    logits = np.full((4, 512), -50.0)
    logits[0, ord("C")] = logits[1, ord("G")] = logits[2, ord("T")] = 50.0
    backend = Fake([make_reply(logits)])
    assert backend.score(seq) == pytest.approx(0.0, abs=1e-9)
    sent = json.loads(backend.bodies[0])
    assert sent == {"sequence": "ACGT", "output_layers": ["unembed"]}


def test_garbage_reply_is_a_clean_error():
    for raw in (b"not json", json.dumps({"x": 1}).encode(), json.dumps({"data": "AAAA"}).encode()):
        with pytest.raises(BackendError):
            NvidiaBackend.decode_response(raw)


def test_wrong_array_shape_rejected():
    buf = io.BytesIO()
    np.save(buf, np.zeros(7))  # a 1-D array is not (length, vocab)
    zbuf = io.BytesIO()
    with zipfile.ZipFile(zbuf, "w") as z:
        z.writestr("x.npy", buf.getvalue())
    raw = json.dumps({"data": base64.b64encode(zbuf.getvalue()).decode()}).encode()
    with pytest.raises(BackendError, match="shape"):
        NvidiaBackend.decode_response(raw)


@pytest.mark.parametrize("code", [401, 403])
def test_bad_key_fails_fast_without_retry(code):
    backend = Fake([http_error(code)])
    with pytest.raises(BackendError, match="rejected the API key"):
        backend.score("ACGT")
    assert len(backend.bodies) == 1 and backend.sleeps == []


def test_credit_and_validation_errors_fail_fast():
    with pytest.raises(BackendError, match="credits"):
        Fake([http_error(402)]).score("ACGT")
    with pytest.raises(BackendError, match="rejected the request"):
        Fake([http_error(422, body=b'{"error":"bad"}')]).score("ACGT")


def test_429_is_retried_with_retry_after_then_succeeds():
    logits = np.zeros((4, 512))
    backend = Fake([http_error(429, {"Retry-After": "7"}), http_error(503), make_reply(logits)])
    assert backend.score("ACGT") == pytest.approx(3 * -np.log(512))
    assert 7.0 in backend.sleeps
    assert len(backend.bodies) == 3


def test_gives_up_after_max_retries():
    backend = Fake([http_error(429)] * 3, max_retries=2)
    with pytest.raises(BackendError, match="Giving up"):
        backend.score("ACGT")
    assert len(backend.bodies) == 3


def test_network_errors_are_retried():
    logits = np.zeros((4, 512))
    backend = Fake([urllib.error.URLError("dns"), TimeoutError("slow"), make_reply(logits)])
    backend.score("ACGT")
    assert len(backend.bodies) == 3


def test_throttle_spaces_requests_by_rpm():
    logits = np.zeros((4, 512))
    backend = Fake([make_reply(logits)] * 3, rpm=30)  # one request per 2 seconds
    for _ in range(3):
        backend.score("ACGT")
    assert sum(backend.sleeps) == pytest.approx(4.0)


def test_key_is_never_in_repr_or_errors():
    backend = Fake([http_error(401)])
    assert "SECRET" not in repr(backend)
    with pytest.raises(BackendError) as info:
        backend.score("ACGT")
    assert "SECRET" not in str(info.value)


def test_get_api_key_lookup_and_missing():
    assert get_api_key({"NVIDIA_API_KEY": " nvapi-a "}) == "nvapi-a"
    assert get_api_key({"NVCF_KEY": "nvapi-b"}) == "nvapi-b"
    assert set(KEY_ENV_VARS) == {"NVIDIA_API_KEY", "NVCF_KEY"}
    with pytest.raises(BackendError, match="No NVIDIA API key"):
        get_api_key({})

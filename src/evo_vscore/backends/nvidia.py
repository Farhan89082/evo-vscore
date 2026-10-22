"""Backend that calls NVIDIA's hosted Evo 2 7B 'forward' endpoint.

Needs a free API key from https://build.nvidia.com (it starts with ``nvapi-``).
The key is read from an environment variable and is never printed or stored.
"""

import base64
import io
import json
import os
import time
import urllib.error
import urllib.request
import zipfile

import numpy as np

from .base import Backend, BackendError, sequence_logprob

DEFAULT_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-7b/forward"
KEY_ENV_VARS = ("NVIDIA_API_KEY", "NVCF_KEY")
_MAX_UNZIPPED_BYTES = 512 * 1024 * 1024


def get_api_key(environ=None):
    environ = os.environ if environ is None else environ
    for var in KEY_ENV_VARS:
        value = environ.get(var, "").strip()
        if value:
            return value
    raise BackendError(
        "No NVIDIA API key found. Get a free key at https://build.nvidia.com and set it "
        "for this terminal session, for example in PowerShell:\n"
        '    $env:NVIDIA_API_KEY = Read-Host "Paste your key"\n'
        "(the variable name can be NVIDIA_API_KEY or NVCF_KEY)."
    )


class NvidiaBackend(Backend):
    name = "nvidia"
    cache_id = "nvidia:evo2-7b:unembed:v1"

    def __init__(self, api_key, url=DEFAULT_URL, rpm=30, max_retries=6, timeout=300,
                 sleep=time.sleep, clock=time.monotonic):
        if rpm <= 0:
            raise ValueError("rpm must be positive")
        self._api_key = api_key
        self.url = url
        self.min_interval = 60.0 / rpm
        self.max_retries = max_retries
        self.timeout = timeout
        self._sleep = sleep
        self._clock = clock
        self._last_request = None
        self.calls = 0

    def __repr__(self):  # never include the key
        return "NvidiaBackend(url=%r, min_interval=%.2fs)" % (self.url, self.min_interval)

    # -- network ------------------------------------------------------------
    def _throttle(self):
        if self._last_request is not None:
            wait = self.min_interval - (self._clock() - self._last_request)
            if wait > 0:
                self._sleep(wait)
        self._last_request = self._clock()

    def _post(self, body):
        request = urllib.request.Request(
            self.url,
            data=body,
            headers={
                "Authorization": "Bearer " + self._api_key,
                "Content-Type": "application/json",
                "User-Agent": "evo-vscore/0.1",
            },
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            return response.read()

    @staticmethod
    def _retry_after(error):
        try:
            value = error.headers.get("Retry-After")
            return float(value) if value else None
        except (AttributeError, TypeError, ValueError):
            return None

    def _request_with_retries(self, body):
        for attempt in range(self.max_retries + 1):
            self._throttle()
            self.calls += 1
            try:
                return self._post(body)
            except urllib.error.HTTPError as error:
                code = error.code
                detail = error.read()[:300].decode("utf-8", "replace")
                if code in (401, 403):
                    raise BackendError(
                        "NVIDIA rejected the API key (HTTP %d). Check that the key is correct "
                        "and still valid on build.nvidia.com." % code
                    )
                if code == 402:
                    raise BackendError(
                        "NVIDIA reports no remaining credits or payment required (HTTP 402). "
                        "Check your account on build.nvidia.com."
                    )
                if code in (400, 413, 422):
                    raise BackendError("NVIDIA rejected the request (HTTP %d): %s" % (code, detail))
                if code != 429 and code < 500:
                    raise BackendError("Unexpected response from NVIDIA (HTTP %d): %s" % (code, detail))
                reason = "rate limited (HTTP 429)" if code == 429 else "server error (HTTP %d)" % code
                wait = self._retry_after(error)
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as error:
                reason = "network problem (%s)" % error
                wait = None
            if attempt == self.max_retries:
                raise BackendError("Giving up after %d attempts: %s" % (attempt + 1, reason))
            if wait is None:
                wait = min(60.0, 2.0 * (2 ** attempt))
            self._sleep(wait)

    # -- response decoding ----------------------------------------------------
    @staticmethod
    def decode_response(raw):
        """Turn the API reply (JSON -> base64 zip -> .npy) into a (length, vocab) array."""
        try:
            payload = json.loads(raw)
            archive = zipfile.ZipFile(io.BytesIO(base64.b64decode(payload["data"])))
            names = [n for n in archive.namelist() if n.endswith(".npy")]
            if not names:
                raise KeyError("no .npy file in the reply")
            info = archive.getinfo(names[0])
            if info.file_size > _MAX_UNZIPPED_BYTES:
                raise ValueError("reply is unexpectedly large")
            array = np.load(io.BytesIO(archive.read(names[0])), allow_pickle=False)
        except BackendError:
            raise
        except Exception as error:
            raise BackendError("Could not read the NVIDIA reply (unexpected format): %s" % error)
        array = np.squeeze(array)
        if array.ndim != 2:
            raise BackendError("Unexpected score array shape %s" % (tuple(array.shape),))
        return array.astype(np.float64)

    def score(self, seq):
        body = json.dumps({"sequence": seq, "output_layers": ["unembed"]}).encode("utf-8")
        logits = self.decode_response(self._request_with_retries(body))
        return sequence_logprob(logits, seq)

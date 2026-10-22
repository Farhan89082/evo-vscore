from .base import Backend, BackendError, sequence_logprob
from .mock import MockBackend
from .nvidia import NvidiaBackend, get_api_key

#: name -> class, used by the command line
BACKENDS = {"mock": MockBackend, "nvidia": NvidiaBackend}

__all__ = [
    "Backend",
    "BackendError",
    "sequence_logprob",
    "MockBackend",
    "NvidiaBackend",
    "get_api_key",
    "BACKENDS",
]

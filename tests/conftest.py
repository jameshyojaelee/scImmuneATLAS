import os
import site
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

try:
    user_site = site.getusersitepackages()
except Exception:
    user_site = []

if isinstance(user_site, str):
    user_sites = [user_site]
else:
    user_sites = list(user_site)

for path in user_sites:
    if path not in sys.path:
        sys.path.append(path)

import pytest

try:
    import ssl  # noqa: F401
except ImportError as exc:  # pragma: no cover - environment dependent
    pytest.skip(
        "OpenSSL support is unavailable (ssl import failed: {}). "
        "Activate the 'immune-atlas' Conda environment or install openssl/cryptography before running tests.".format(exc),
        allow_module_level=True,
    )

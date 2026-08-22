"""Machine-specific paths for the N40 gm/ID pipeline.

Reads key=value pairs from ``machine.env`` in this directory (lines starting
with ``#`` and empty lines ignored). Values may reference ``$HOME`` / other
keys (``$SIMDIR/raw``), resolved against the process environment plus the
config itself. Environment variables take precedence over the file, and the
file takes precedence over the built-in defaults.

All project code must obtain machine paths from here (via ``PATHS`` or
:func:`require`) instead of hardcoding directories or usernames.
"""
import os
import string
from pathlib import Path

_ENV_FILE = Path(__file__).resolve().parent / "machine.env"

_DEFAULTS = {
    "PDK_MODEL_FILE": "",
    "SPECTRE_BIN": "spectre",
    "PYTHON_BIN": "python3",
    "SIMDIR": "$HOME/simulation",
    "RAWDIR": "$HOME/simulation/raw",
    "MATDIR": "$HOME/simulation/out",
    "LOGDIR": "$HOME/simulation/logs",
    "SCRIPTDIR": str(Path(__file__).resolve().parent),
    "TMPDIR": "$HOME/tmp",
}


def _expand(value, namespace):
    """Expand ``$VAR`` / ``${VAR}`` references against ``namespace``."""
    try:
        return string.Template(value).safe_substitute(namespace)
    except Exception:
        return value


def _load():
    ns = dict(os.environ)
    for key, val in _DEFAULTS.items():
        ns[key] = _expand(val, ns)
    if _ENV_FILE.exists():
        for line in _ENV_FILE.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, val = line.split("=", 1)
            ns[key.strip()] = _expand(val.strip(), ns)
    return {key: ns.get(key, _DEFAULTS[key]) for key in _DEFAULTS}


PATHS = _load()


def require(key):
    """Return the value for ``key``, raising if it is empty."""
    val = PATHS[key]
    if not val:
        raise RuntimeError(
            f"{key} is not set in machine.env (copy machine.env.example "
            f"and fill it in)"
        )
    return val

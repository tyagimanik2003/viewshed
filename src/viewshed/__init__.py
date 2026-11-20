"""Correct, reproducible viewsheds from latitude/longitude, anywhere on Earth.

    >>> from viewshed import viewshed
    >>> r = viewshed(77.015, 32.329, radius=5000)        # lon, lat, metres
    >>> print(r.summary())
"""

from .core import INVISIBLE, MAYBE_VISIBLE, NODATA, OUT_OF_RANGE, VISIBLE, ViewshedResult, viewshed
from .ztv import ZTVResult, ztv

__all__ = ["viewshed", "ViewshedResult", "ztv", "ZTVResult",
           "VISIBLE", "INVISIBLE", "MAYBE_VISIBLE", "OUT_OF_RANGE", "NODATA"]
__version__ = "0.1.0.dev0"

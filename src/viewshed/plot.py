"""Quick-look map of a ViewshedResult: hillshade, visible area, range ring, sector."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .core import INVISIBLE, NODATA, OUT_OF_RANGE, VISIBLE, ViewshedResult

VISIBLE_RGB = (0.98, 0.62, 0.11)   # amber: reads as "lit" on grey terrain
SECTOR_RGB = (0.06, 0.35, 0.45)


def hillshade(dem: np.ndarray, cell: float, azimuth: float = 315.0, altitude: float = 45.0) -> np.ndarray:
    z = np.where(dem <= -32767, np.nan, dem).astype(float)
    dy, dx = np.gradient(z, cell)
    slope = np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    az, alt = np.radians(360 - azimuth + 90), np.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.nan_to_num(np.clip(hs, 0, 1), nan=1.0)


def plot(result: ViewshedResult, out: str | Path | None = None, title: str | None = None, ax=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    gt = result.geotransform
    h, w = result.codes.shape
    extent = [gt[0] / 1000, (gt[0] + w * gt[1]) / 1000, (gt[3] + h * gt[5]) / 1000, gt[3] / 1000]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 7.2), dpi=130)
    ax.imshow(hillshade(result.dem, gt[1]), cmap="gray", vmin=0, vmax=1, extent=extent, interpolation="bilinear")

    rgba = np.zeros((h, w, 4))
    rgba[result.codes == VISIBLE] = (*VISIBLE_RGB, 0.62)
    rgba[result.codes == OUT_OF_RANGE] = (1, 1, 1, 0.55)
    rgba[result.codes == NODATA] = (0.8, 0.1, 0.1, 0.5)
    ax.imshow(rgba, extent=extent, interpolation="nearest")

    r = result.radius / 1000
    t = np.linspace(0, 2 * np.pi, 361)
    ax.plot(r * np.sin(t), r * np.cos(t), color="black", lw=0.8, alpha=0.6)
    sa, ea = result.params.get("start_angle"), result.params.get("end_angle")
    if sa is not None:
        for a in (sa, ea):
            ax.plot([0, r * np.sin(np.radians(a))], [0, r * np.cos(np.radians(a))], color=SECTOR_RGB, lw=1.4)
    ax.plot(0, 0, marker="*", ms=14, mec="white", mfc="black", mew=1.2, zorder=5)
    ax.annotate("N", xy=(0.95, 0.93), xycoords="axes fraction", ha="center", fontsize=11, weight="bold")
    ax.annotate("", xy=(0.95, 0.91), xytext=(0.95, 0.83), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="black"))
    ax.set_xlabel("km east of observer")
    ax.set_ylabel("km north of observer")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    lon, lat = result.observer
    ax.set_title(title or f"Viewshed from {lat:.3f}°N {lon:.3f}°E", fontsize=11)
    ax.text(0.01, 0.01, f"{result.visible_km2:.1f} km² visible ({100 * result.visible_fraction:.0f}%)",
            transform=ax.transAxes, fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.legend(handles=[Patch(fc=(*VISIBLE_RGB, 0.62), label="visible"),
                       Patch(fc=(0.55, 0.55, 0.55), label="not visible"),
                       Patch(fc=(0.95, 0.95, 0.95), ec="0.7", label="outside range / sector")],
              loc="upper left", fontsize=8, framealpha=0.9)
    if own and out:
        fig.tight_layout()
        fig.savefig(out)
        plt.close(fig)
    return ax

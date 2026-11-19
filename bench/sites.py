"""ViewshedBench pilot sites, chosen to span relief classes in north India.

Coordinates are approximate town/feature centres (±0.01°). The relief class used in the
analysis is computed from the DEM itself (elevation range within the radius), so the
region labels below are only for reading.
"""
SITES = {
    # high Himalaya
    "Kangra-Lahaul ridge": (77.015, 32.329, "high Himalaya"),
    "Samudra Tapu lake, Lahaul": (77.54, 32.50, "high Himalaya"),
    "Kaza, Spiti": (78.07, 32.23, "high Himalaya"),
    "Padum, Zanskar": (76.88, 33.47, "high Himalaya"),
    "Kalpa, Kinnaur": (78.25, 31.54, "high Himalaya"),
    "Kedarnath, Garhwal": (79.07, 30.73, "high Himalaya"),
    # high plateau / broad valleys
    "Leh, Ladakh": (77.58, 34.16, "Ladakh plateau"),
    "Hanle, Changthang": (78.97, 32.78, "Ladakh plateau"),
    # foothills
    "Dharamshala": (76.32, 32.22, "foothills"),
    "Mandi": (76.93, 31.71, "foothills"),
    "Solan": (77.10, 30.91, "foothills"),
    # plains
    "Kurukshetra": (76.85, 29.97, "plains"),
    "Ludhiana": (75.85, 30.90, "plains"),
    "Bikaner dunes, Thar": (73.31, 28.02, "plains"),
}
NAMED = dict(SITES)

# Stratified random sites (bench/sample_sites.py), appended when present.
import json as _json
from pathlib import Path as _Path

_rand = _Path(__file__).with_name("sites_random.json")
RANDOM = {k: tuple(v) for k, v in _json.loads(_rand.read_text(encoding="utf-8")).items()} if _rand.exists() else {}
SITES.update(RANDOM)

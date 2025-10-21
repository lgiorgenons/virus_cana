"""Tools for downloading and analysing satellite imagery to monitor Sugarcane Wilt Syndrome.

This module provides a small command line utility to:

* query and download Sentinel-2 Level-2A products using the Copernicus Open Access Hub;
* extract the necessary spectral bands; and
* compute vegetation/water stress indicators such as NDVI, NDWI and MSI.

The workflow is inspired by the agronomic discussion available in the
file "Análise Detalhada e Abrangente da Síndrome da Murcha da Cana-de-Açúcar.md".

The implementation favours composable functions so that the script can be
orchestrated from notebooks or larger data pipelines. For the sake of simplicity
we only download the first product that matches the query criteria. The
``sentinelsat`` package is used for the download while ``rasterio`` and ``numpy``
perform raster handling and analytics.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, Optional, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from sentinelsat import SentinelAPI, geojson_to_wkt


_LOGGER = logging.getLogger(__name__)


DEFAULT_API_URL = "https://apihub.copernicus.eu/apihub"


@dataclass
class DownloadConfig:
    """Configuration parameters required to query Sentinel-2 scenes."""

    username: str
    password: str
    api_url: str = DEFAULT_API_URL


def authenticate_from_env(
    prefix: str = "SENTINEL",
    *,
    username: Optional[str] = None,
    password: Optional[str] = None,
    api_url: Optional[str] = None,
) -> DownloadConfig:
    """Build :class:`DownloadConfig` from CLI arguments and environment variables.

    Parameters
    ----------
    prefix:
        Prefix used to read the credentials, e.g. ``SENTINEL_USERNAME``.

    Returns
    -------
    DownloadConfig
        Credential bundle ready for use by :class:`SentinelAPI`.
    """

    env_username = os.environ.get(f"{prefix}_USERNAME")
    env_password = os.environ.get(f"{prefix}_PASSWORD")
    env_api_url = os.environ.get(f"{prefix}_API_URL")

    username = username or env_username
    password = password or env_password
    api_url = api_url or env_api_url or DEFAULT_API_URL

    if not username or not password:
        raise RuntimeError(
            "Missing Sentinel Hub credentials. Set the "
            f"{prefix}_USERNAME and {prefix}_PASSWORD environment variables."
        )

    return DownloadConfig(username=username, password=password, api_url=api_url)


@dataclass
class AreaOfInterest:
    """Simple container for a polygon of interest."""

    geometry: Dict

    @classmethod
    def from_geojson(cls, geojson_path: Path) -> "AreaOfInterest":
        with Path(geojson_path).open("r", encoding="utf-8") as file:
            geometry = json.load(file)
        return cls(geometry=geometry)


def query_latest_product(
    api: SentinelAPI,
    area: AreaOfInterest,
    start_date: str,
    end_date: str,
    cloud_cover: Tuple[int, int] = (0, 20),
) -> Optional[Dict]:
    """Query the latest Sentinel-2 product matching the constraints."""

    footprint = geojson_to_wkt(area.geometry)
    products = api.query(
        footprint,
        date=(start_date, end_date),
        platformname="Sentinel-2",
        producttype="S2MSI2A",
        cloudcoverpercentage=cloud_cover,
    )
    if not products:
        return None

    # Sentinel API returns a dict keyed by UUID. We select the most recent scene.
    return products[next(iter(sorted(products, key=lambda k: products[k]["ingestiondate"], reverse=True)))]


def download_product(api: SentinelAPI, product: Dict, target_dir: Path) -> Path:
    """Download the provided product into *target_dir* and return the path."""

    target_dir.mkdir(parents=True, exist_ok=True)
    product_path = api.download(product["uuid"], directory_path=str(target_dir))
    return Path(product_path["path"])


SENTINEL_BANDS = {
    "B02": "blue",
    "B03": "green",
    "B04": "red",
    "B08": "nir",
    "B11": "swir1",
    "B12": "swir2",
}


def _infer_product_name(product_path: Path, fallback: Optional[str] = None) -> str:
    """Infer a human-friendly product name from a SAFE archive path."""

    stem = product_path.stem
    if stem.endswith(".SAFE"):
        stem = stem[:-5]
    return stem or fallback or product_path.name


def _locate_band(safe_root: Path, band: str) -> Path:
    patterns = [
        f"**/IMG_DATA/*/*_{band}_*.jp2",
        f"**/IMG_DATA/*_{band}_*.jp2",
        f"**/IMG_DATA/**/*_{band}_*.jp2",
    ]
    for pattern in patterns:
        matches = list(safe_root.glob(pattern))
        if matches:
            # Prefer higher resolution files first (10m -> 20m -> 60m)
            matches.sort(key=lambda p: ("10m" not in p.name, "20m" in p.name, p.name))
            return matches[0]
    raise FileNotFoundError(f"Band {band} not found inside {safe_root}")


def extract_bands_from_safe(safe_archive: Path, destination: Path) -> Dict[str, Path]:
    """Extract relevant bands from a SAFE archive or directory.

    The function is tolerant to either a ``.zip`` file or the already extracted
    SAFE directory. Each requested band is converted to GeoTIFF and stored inside
    ``destination`` with simplified naming, e.g. ``nir.tif``. Bands are
    resampled later on demand so that the calling code can control the reference
    resolution.
    """

    destination.mkdir(parents=True, exist_ok=True)

    tmp_dir: Optional[tempfile.TemporaryDirectory[str]] = None
    if safe_archive.suffix == ".zip":
        tmp_dir = tempfile.TemporaryDirectory(prefix="safe_")
        tmp_path = Path(tmp_dir.name)
        _LOGGER.info("Extracting SAFE archive %s", safe_archive)
        with zipfile.ZipFile(safe_archive) as archive:
            archive.extractall(tmp_path)
        safe_root = next(tmp_path.glob("*.SAFE"))
    else:
        safe_root = safe_archive

    extracted: Dict[str, Path] = {}
    for band_id, alias in SENTINEL_BANDS.items():
        try:
            jp2_path = _locate_band(safe_root, band_id)
        except FileNotFoundError:
            _LOGGER.warning("Band %s not found in SAFE structure", band_id)
            continue

        tif_path = destination / f"{alias}.tif"
        with rasterio.open(jp2_path) as src:
            profile = src.profile
            data = src.read(1)

        profile.update(driver="GTiff")
        with rasterio.open(tif_path, "w", **profile) as dst:
            dst.write(data, 1)
        extracted[alias] = tif_path
    if tmp_dir is not None:
        tmp_dir.cleanup()

    return extracted


def _compute_index(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Utility used by spectral indices to avoid division-by-zero."""

    mask = denominator == 0
    denominator = np.where(mask, np.nan, denominator)
    index = numerator / denominator
    return np.where(np.isnan(index), 0, index)


def load_raster(
    path: Path,
    reference_path: Optional[Path] = None,
) -> Tuple[np.ndarray, rasterio.Affine, rasterio.crs.CRS]:
    """Load a single band raster file with optional resampling."""

    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float32)
        transform = src.transform
        crs = src.crs
        height = src.height
        width = src.width

    if reference_path is not None:
        with rasterio.open(reference_path) as ref:
            if (transform != ref.transform) or (height != ref.height) or (width != ref.width):
                destination = np.empty((ref.height, ref.width), dtype=np.float32)
                reproject(
                    source=data,
                    destination=destination,
                    src_transform=transform,
                    src_crs=crs,
                    dst_transform=ref.transform,
                    dst_crs=ref.crs,
                    resampling=Resampling.bilinear,
                )
                data = destination
                transform = ref.transform
                crs = ref.crs

    return data, transform, crs


def compute_ndvi(nir: np.ndarray, red: np.ndarray) -> np.ndarray:
    """Compute the NDVI index."""

    return _compute_index(nir - red, nir + red)


def compute_ndwi(nir: np.ndarray, swir: np.ndarray) -> np.ndarray:
    """Compute the NDWI index."""

    return _compute_index(nir - swir, nir + swir)


def compute_msi(nir: np.ndarray, swir: np.ndarray) -> np.ndarray:
    """Compute the Moisture Stress Index (MSI)."""

    return _compute_index(swir, nir)


def save_raster(array: np.ndarray, template_path: Path, destination: Path) -> Path:
    """Persist an index using the metadata from *template_path*."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(template_path) as src:
        meta = src.meta.copy()
    meta.update(dtype=rasterio.float32, count=1)
    with rasterio.open(destination, "w", **meta) as dst:
        dst.write(array.astype(rasterio.float32), 1)
    return destination


@dataclass(frozen=True)
class IndexSpec:
    """Metadata describing how to compute a spectral index."""

    bands: Tuple[str, ...]
    func: Callable[..., np.ndarray]


INDEX_SPECS: Dict[str, IndexSpec] = {
    "ndvi": IndexSpec(bands=("nir", "red"), func=compute_ndvi),
    "ndwi": IndexSpec(bands=("nir", "swir1"), func=compute_ndwi),
    "msi": IndexSpec(bands=("nir", "swir1"), func=compute_msi),
}


def analyse_scene(
    band_paths: Dict[str, Path],
    output_dir: Path,
    indices: Optional[Iterable[str]] = None,
) -> Dict[str, Path]:
    """Compute spectral indices for the scene and store them in *output_dir*."""

    requested = list(dict.fromkeys(indices)) if indices is not None else list(INDEX_SPECS.keys())
    if not requested:
        raise ValueError("No spectral indices requested.")

    unknown = sorted(set(requested) - INDEX_SPECS.keys())
    if unknown:
        raise ValueError(f"Unsupported indices requested: {', '.join(unknown)}")

    required = set()
    for name in requested:
        required.update(INDEX_SPECS[name].bands)

    missing_bands = required - band_paths.keys()
    if missing_bands:
        raise RuntimeError(f"Missing bands for analysis: {', '.join(sorted(missing_bands))}")

    nir_data, transform, crs = load_raster(band_paths["nir"])
    band_arrays: Dict[str, np.ndarray] = {"nir": nir_data}

    for band in required - {"nir"}:
        data, _, _ = load_raster(band_paths[band], reference_path=band_paths["nir"])
        band_arrays[band] = data

    outputs: Dict[str, Path] = {}
    for name in requested:
        spec = INDEX_SPECS[name]
        arrays = [band_arrays[band] for band in spec.bands]
        result = spec.func(*arrays)
        outputs[name] = save_raster(result, band_paths["nir"], output_dir / f"{name}.tif")
    return outputs


def _parse_arguments(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geojson", type=Path, help="Path to AOI GeoJSON polygon (required when downloading)")
    parser.add_argument("--start-date", help="Start date YYYY-MM-DD (required when downloading)")
    parser.add_argument("--end-date", help="End date YYYY-MM-DD (required when downloading)")
    parser.add_argument("--cloud", type=int, nargs=2, default=(0, 20), help="Acceptable cloud cover percentage range")
    parser.add_argument("--download-dir", type=Path, default=Path("data/raw"), help="Directory to store downloaded products")
    parser.add_argument("--workdir", type=Path, default=Path("data/processed"), help="Directory for processed outputs")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--username", help="Sentinel API username (overrides SENTINEL_USERNAME env var)")
    parser.add_argument("--password", help="Sentinel API password (overrides SENTINEL_PASSWORD env var)")
    parser.add_argument("--api-url", help=f"Sentinel API URL (defaults to {DEFAULT_API_URL})")
    parser.add_argument(
        "--safe-path",
        type=Path,
        help="Existing SAFE archive (.zip) or directory to analyse instead of downloading a new product",
    )
    parser.add_argument(
        "--indices",
        nargs="+",
        choices=sorted(INDEX_SPECS.keys()),
        help="Subset of spectral indices to compute (default: all available)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> None:
    args = _parse_arguments(argv)
    logging.basicConfig(level=args.log_level)

    selected_indices = list(dict.fromkeys(args.indices or list(INDEX_SPECS.keys())))

    if args.safe_path:
        product_path = args.safe_path.expanduser().resolve()
        if not product_path.exists():
            _LOGGER.error("Provided SAFE path does not exist: %s", product_path)
            return
        product_title = _infer_product_name(product_path)
        _LOGGER.info("Using existing SAFE product at %s", product_path)
    else:
        missing_args = [name for name, value in (("geojson", args.geojson), ("start-date", args.start_date), ("end-date", args.end_date)) if value is None]
        if missing_args:
            _LOGGER.error(
                "Missing required arguments for download: %s",
                ", ".join(missing_args),
            )
            return
        config = authenticate_from_env(username=args.username, password=args.password, api_url=args.api_url)
        assert args.geojson is not None
        geojson_path = args.geojson.expanduser().resolve()
        area = AreaOfInterest.from_geojson(geojson_path)

        api = SentinelAPI(config.username, config.password, config.api_url)
        product = query_latest_product(api, area, args.start_date, args.end_date, tuple(args.cloud))  # type: ignore[arg-type]
        if not product:
            _LOGGER.error("No Sentinel-2 product found for the provided parameters.")
            return

        _LOGGER.info("Downloading product %s", product["title"])
        product_path = download_product(api, product, args.download_dir)
        product_title = product.get("title") or _infer_product_name(product_path)

    _LOGGER.info("Extracting spectral bands")
    bands = extract_bands_from_safe(product_path, args.workdir / product_title)  # type: ignore[arg-type]

    _LOGGER.info("Computing spectral indices: %s", ", ".join(selected_indices))
    outputs = analyse_scene(bands, args.workdir / product_title / "indices", indices=selected_indices)

    for name, path in outputs.items():
        _LOGGER.info("Generated %s at %s", name.upper(), path)


if __name__ == "__main__":
    main()

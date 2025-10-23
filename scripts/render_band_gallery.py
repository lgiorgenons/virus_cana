"""Generate an HTML gallery showing all Sentinel-2 bands for a processed product.

The gallery renders each band (coastal, blue, green, red, red-edge, NIR, SWIR) using
a common stretch and optionally restricts the view to the supplied GeoJSON AOI.
This is a quick way to inspect every band without opening a GIS.
"""
from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import base64
import io

import numpy as np
import rasterio
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio.warp import calculate_default_transform, reproject, transform_bounds
from rasterio.windows import from_bounds

TARGET_CRS = "EPSG:4326"

BAND_ORDER = OrderedDict(
    [
        ("coastal", "B01 Coastal/Aerosol"),
        ("blue", "B02 Blue"),
        ("green", "B03 Green"),
        ("red", "B04 Red"),
        ("rededge1", "B05 Red-edge 1"),
        ("rededge2", "B06 Red-edge 2"),
        ("rededge3", "B07 Red-edge 3"),
        ("nir", "B08 NIR"),
        ("rededge4", "B8A Narrow NIR"),
        ("water_vapor", "B09 Water Vapour"),
        ("cirrus", "B10 Cirrus"),
        ("swir1", "B11 SWIR 1"),
        ("swir2", "B12 SWIR 2"),
    ]
)


def _load_geojson(path: Optional[Path]) -> Optional[Dict]:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as file:
        return json.load(file)


def _extract_bbox(geojson: Dict) -> Tuple[float, float, float, float]:
    def extract(geometry: Dict) -> Dict:
        gtype = geometry.get("type")
        if gtype == "FeatureCollection":
            features = geometry.get("features", [])
            if not features:
                raise ValueError("GeoJSON sem features.")
            return extract(features[0])
        if gtype == "Feature":
            return extract(geometry["geometry"])
        return geometry

    coords = []
    geometry = extract(geojson)
    gtype = geometry.get("type")
    if gtype == "Polygon":
        coords = geometry["coordinates"][0]
    elif gtype == "MultiPolygon":
        for poly in geometry["coordinates"]:
            coords.extend(poly[0])
    else:
        raise ValueError(f"Tipo de geometria nao suportado: {gtype}")

    lons = [pt[0] for pt in coords]
    lats = [pt[1] for pt in coords]
    return min(lons), min(lats), max(lons), max(lats)


def _read_band(
    path: Path,
    clip_bounds_wgs84: Optional[Tuple[float, float, float, float]] = None,
    resampling: Resampling = Resampling.bilinear,
) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
    with rasterio.open(path) as src:
        if clip_bounds_wgs84 is not None:
            left, bottom, right, top = transform_bounds(
                TARGET_CRS, src.crs, *clip_bounds_wgs84, densify_pts=21
            )
            left = max(left, src.bounds.left)
            bottom = max(bottom, src.bounds.bottom)
            right = min(right, src.bounds.right)
            top = min(top, src.bounds.top)
            window = from_bounds(left, bottom, right, top, transform=src.transform).round_offsets().round_lengths()
            data = src.read(1, window=window).astype(np.float32)
            src_transform = src.window_transform(window)
            src_bounds = (left, bottom, right, top)
        else:
            data = src.read(1).astype(np.float32)
            src_transform = src.transform
            src_bounds = src.bounds

        dst_transform, width, height = calculate_default_transform(
            src.crs, TARGET_CRS, data.shape[1], data.shape[0], *src_bounds
        )
        destination = np.full((height, width), np.nan, dtype=np.float32)
        reproject(
            source=data,
            destination=destination,
            src_transform=src_transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=TARGET_CRS,
            src_nodata=src.nodata if src.nodata is not None else None,
            dst_nodata=np.nan,
            resampling=resampling,
        )
    bounds = array_bounds(destination.shape[0], destination.shape[1], dst_transform)
    return destination, bounds


def _figure_to_base64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"


def build_gallery(
    product_dir: Path,
    output_html: Path,
    geojson_path: Optional[Path],
    stretch_percentiles: Tuple[float, float] = (2, 98),
) -> Path:
    geojson = _load_geojson(geojson_path)
    clip_bounds = _extract_bbox(geojson) if geojson else None

    band_entries = []
    extent: Optional[Tuple[float, float, float, float]] = None
    for band_key, band_label in BAND_ORDER.items():
        band_path = product_dir / f"{band_key}.tif"
        if not band_path.exists():
            continue

        data, bounds = _read_band(band_path, clip_bounds_wgs84=clip_bounds)
        if extent is None:
            extent = bounds
        valid = np.isfinite(data)
        if not np.any(valid):
            continue

        vmin, vmax = np.nanpercentile(data[valid], stretch_percentiles)
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)

        fig, ax = plt.subplots(figsize=(4, 4))
        if extent is not None:
            min_lon, min_lat, max_lon, max_lat = extent
            ax.imshow(norm(data), cmap="gray", extent=[min_lon, max_lon, min_lat, max_lat], origin="lower")
            ax.set_xlim(min_lon, max_lon)
            ax.set_ylim(min_lat, max_lat)
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.imshow(norm(data), cmap="gray")
        if geojson:
            for geom in geojson.get("features", [geojson]):
                geom_data = geom["geometry"] if geom.get("type") == "Feature" else geom
                xs = [coord[0] for coord in geom_data["coordinates"][0]]
                ys = [coord[1] for coord in geom_data["coordinates"][0]]
                ax.plot(xs, ys, color="cyan", linewidth=1)
        ax.set_title(f"{band_label}\n({band_key}.tif)\n[{vmin:.2f}, {vmax:.2f}]")
        ax.axis("off")

        img_base64 = _figure_to_base64(fig)
        band_entries.append((band_label, band_key, img_base64))

    html_parts = [
        "<html><head><meta charset='utf-8'><title>Sentinel-2 Band Gallery</title>",
        "<style>body { font-family: sans-serif; background:#111; color:#eee;}",
        ".grid { display: flex; flex-wrap: wrap; gap: 12px;}",
        ".band { background:#222; padding: 8px; border-radius: 6px; width: 320px; text-align:center;}",
        ".band img { width:100%; border-radius:4px; }",
        "</style></head><body>",
        f"<h1>Sentinel-2 Band Gallery — {product_dir.name}</h1>",
    ]
    if geojson_path:
        html_parts.append(f"<p>Recorte aplicado com base no GeoJSON: {geojson_path}</p>")
    html_parts.append("<div class='grid'>")
    for label, key, img in band_entries:
        html_parts.append(f"<div class='band'><img src='{img}' alt='{key}'/><div>{label}</div></div>")
    html_parts.append("</div></body></html>")

    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text("\n".join(html_parts), encoding="utf-8")
    return output_html


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--product-dir", type=Path, required=True, help="Diretorio processado com as bandas.")
    parser.add_argument("--output", type=Path, default=Path("mapas/band_gallery.html"), help="Arquivo HTML de saída.")
    parser.add_argument("--geojson", type=Path, help="GeoJSON para recorte (opcional).")
    parser.add_argument("--stretch", type=float, nargs=2, default=(2, 98), help="Percentis para stretch (ex.: 2 98).")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    output = build_gallery(
        product_dir=args.product_dir.expanduser().resolve(),
        output_html=args.output.expanduser().resolve(),
        geojson_path=args.geojson.expanduser().resolve() if args.geojson else None,
        stretch_percentiles=tuple(args.stretch),
    )
    print(f"Galeria salva em {output}")


if __name__ == "__main__":
    main()

"""Generate a side-by-side map comparing true-color imagery with a spectral index.

Left pane shows an online basemap (Esri World Imagery by default).
Right pane overlays the spectral index rendered from a GeoTIFF.
Both panes are synchronized for zoom/pan to aid visual inspection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import folium
import numpy as np
import rasterio
from branca.colormap import LinearColormap
from folium.plugins import DualMap
from matplotlib import colormaps, colors
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio.warp import calculate_default_transform, reproject
from scipy.ndimage import gaussian_filter

TARGET_CRS = "EPSG:4326"


def _load_raster(path: Path) -> Tuple[np.ndarray, rasterio.Affine, Tuple[float, float, float, float]]:
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            data = np.where(data == nodata, np.nan, data)
        elif np.ma.isMaskedArray(data):
            data = data.filled(np.nan)

        dst_transform, width, height = calculate_default_transform(
            src.crs, TARGET_CRS, src.width, src.height, *src.bounds
        )
        destination = np.full((height, width), np.nan, dtype=np.float32)

        reproject(
            source=data,
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=TARGET_CRS,
            src_nodata=nodata if nodata is not None else None,
            dst_nodata=np.nan,
            resampling=Resampling.nearest,
        )

    bounds = array_bounds(destination.shape[0], destination.shape[1], dst_transform)
    return destination, dst_transform, bounds


def _load_geojson(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as file:
        return json.load(file)


def _extract_bounds(geojson_data: Dict) -> Optional[Tuple[float, float, float, float]]:
    def extract(geometry: Dict) -> Dict:
        if geometry.get("type") == "FeatureCollection":
            features = geometry.get("features", [])
            if not features:
                raise ValueError("GeoJSON vazio.")
            return extract(features[0])
        if geometry.get("type") == "Feature":
            return extract(geometry["geometry"])
        return geometry

    geometry = extract(geojson_data)
    coords = geometry.get("coordinates")
    gtype = geometry.get("type")
    if gtype == "Polygon":
        points = coords[0]
    elif gtype == "MultiPolygon":
        points = [pt for polygon in coords for pt in polygon[0]]
    else:
        return None

    lons = [pt[0] for pt in points]
    lats = [pt[1] for pt in points]
    return min(lons), min(lats), max(lons), max(lats)


def _apply_unsharp_mask(array: np.ndarray, radius: float, amount: float) -> np.ndarray:
    mask = np.isfinite(array)
    if not np.any(mask):
        return array

    filled = np.where(mask, array, float(np.nanmean(array[mask])))
    blurred = gaussian_filter(filled, sigma=radius, mode="nearest")
    sharpened = filled + amount * (filled - blurred)
    return np.where(mask, sharpened, np.nan)


def _prepare_rgba(
    data: np.ndarray,
    cmap_name: str,
    opacity: float,
    vmin: Optional[float],
    vmax: Optional[float],
) -> Tuple[np.ndarray, float, float]:
    valid = np.isfinite(data)
    if not np.any(valid):
        raise RuntimeError("Raster nao possui pixels validos.")

    finite = data[valid]
    min_value = vmin if vmin is not None else float(np.nanpercentile(finite, 2))
    max_value = vmax if vmax is not None else float(np.nanpercentile(finite, 98))
    if np.isclose(min_value, max_value):
        max_value = min_value + 1e-3

    cmap = colormaps[cmap_name]
    normaliser = colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
    rgba = cmap(normaliser(data))
    rgba[..., 3] = np.where(valid, opacity, 0.0)
    image = (rgba * 255).astype(np.uint8)
    return image, min_value, max_value


def build_comparison_map(
    index_path: Path,
    output_path: Path,
    geojson_paths: Iterable[Path],
    cmap_name: str = "RdYlGn",
    opacity: float = 0.75,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    sharpen: bool = False,
    sharpen_radius: float = 1.2,
    sharpen_amount: float = 1.5,
    basemap: str = "OpenStreetMap",
    attr: Optional[str] = None,
    max_zoom: int = 19,
) -> Path:
    data, transform, bounds = _load_raster(index_path)
    if sharpen:
        data = _apply_unsharp_mask(data, radius=sharpen_radius, amount=sharpen_amount)
    image, min_value, max_value = _prepare_rgba(data, cmap_name, opacity, vmin, vmax)

    overlays = [ _load_geojson(path) for path in geojson_paths ]
    geom_bounds = [_extract_bounds(o) for o in overlays]
    geom_bounds = [b for b in geom_bounds if b is not None]

    min_lon, min_lat, max_lon, max_lat = bounds
    if geom_bounds:
        min_lon = min(min_lon, min(b[0] for b in geom_bounds))
        min_lat = min(min_lat, min(b[1] for b in geom_bounds))
        max_lon = max(max_lon, max(b[2] for b in geom_bounds))
        max_lat = max(max_lat, max(b[3] for b in geom_bounds))

    centre_lat = (min_lat + max_lat) / 2
    centre_lon = (min_lon + max_lon) / 2

    dual_map = DualMap(
        location=[centre_lat, centre_lon],
        zoom_start=14,
        tiles=None,
        max_zoom=max_zoom,
    )

    # Left map: basemap only (true color)
    folium.TileLayer(tiles=basemap, attr=attr, name="Basemap (left)", control=False).add_to(dual_map.m1)
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri World Imagery",
        name="Esri World Imagery",
        overlay=False,
        control=True,
    ).add_to(dual_map.m1)

    # Right map: index overlay + basemap
    folium.TileLayer(tiles=basemap, attr=attr, name="Basemap (right)", control=False).add_to(dual_map.m2)
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri World Imagery",
        name="Esri World Imagery",
        overlay=False,
        control=True,
    ).add_to(dual_map.m2)

    folium.raster_layers.ImageOverlay(
        image=image,
        bounds=[[min_lat, min_lon], [max_lat, max_lon]],
        opacity=1.0,
        name=index_path.stem.upper(),
    ).add_to(dual_map.m2)

    if overlays:
        for data_geo in overlays:
            folium.GeoJson(data=data_geo, name="Area de interesse", style_function=lambda _: {"fillOpacity": 0}).add_to(dual_map.m1)
            folium.GeoJson(data=data_geo, name="Area de interesse", style_function=lambda _: {"fillOpacity": 0}).add_to(dual_map.m2)

    colorbar = LinearColormap(
        [colors.rgb2hex(colormaps[cmap_name](x)) for x in np.linspace(0, 1, 10)],
        vmin=min_value,
        vmax=max_value,
    )
    colorbar.caption = f"{index_path.stem.upper()} (min={min_value:.3f}, max={max_value:.3f})"
    colorbar.add_to(dual_map.m2)

    dual_map.fit_bounds([[min_lat, min_lon], [max_lat, max_lon]])
    folium.LayerControl(position="topright").add_to(dual_map.m1)
    folium.LayerControl(position="topright").add_to(dual_map.m2)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    dual_map.save(str(output_path))
    return output_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True, help="GeoTIFF com o indice (NDVI, NDRE, etc.).")
    parser.add_argument("--output", type=Path, default=Path("mapas/comparativo.html"), help="Arquivo HTML de saída.")
    parser.add_argument("--geojson", type=Path, nargs="*", default=[], help="Arquivos GeoJSON para destacar no mapa.")
    parser.add_argument("--colormap", default="RdYlGn", help="Colormap Matplotlib para renderizar o indice.")
    parser.add_argument("--opacity", type=float, default=0.8, help="Opacidade aplicada ao indice (0-1).")
    parser.add_argument("--vmin", type=float, help="Valor mínimo para normalização.")
    parser.add_argument("--vmax", type=float, help="Valor máximo para normalização.")
    parser.add_argument("--sharpen", action="store_true", help="Aplica unsharp mask ao raster do indice.")
    parser.add_argument("--sharpen-radius", type=float, default=1.2, help="Raio (sigma) da unsharp mask.")
    parser.add_argument("--sharpen-amount", type=float, default=1.5, help="Intensidade do realce de nitidez.")
    parser.add_argument("--tiles", default="OpenStreetMap", help="Tileset base (OpenStreetMap, Stamen Terrain, etc.).")
    parser.add_argument("--tile-attr", default=None, help="Atribuição do tileset base.")
    parser.add_argument("--max-zoom", type=int, default=19, help="Zoom máximo permitido no mapa.")
    return parser.parse_args(argv)

def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    output = build_comparison_map(
        index_path=args.index.expanduser().resolve(),
        output_path=args.output.expanduser().resolve(),
        geojson_paths=[path.expanduser().resolve() for path in args.geojson],
        cmap_name=args.colormap,
        opacity=args.opacity,
        vmin=args.vmin,
        vmax=args.vmax,
        sharpen=args.sharpen,
        sharpen_radius=args.sharpen_radius,
        sharpen_amount=args.sharpen_amount,
        basemap=args.tiles,
        attr=args.tile_attr,
        max_zoom=args.max_zoom,
    )
    print(f"Mapa comparativo salvo em {output}")

if __name__ == "__main__":
    main()

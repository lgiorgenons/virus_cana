"""Generate an interactive HTML map overlaying a spectral índice GeoTIFF.

The script converts the raster to EPSG:4326, applies a colour map and serves the
result as a folium overlay so that it can be inspected in any browser.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Dict, Any

import folium
import numpy as np
import rasterio
from affine import Affine
from branca.colormap import LinearColormap
from matplotlib import colormaps, colors
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import array_bounds
from rasterio.warp import calculate_default_transform, reproject, transform_bounds
from rasterio.windows import from_bounds
from scipy.ndimage import gaussian_filter


DEFAULT_CMAP = "RdYlGn"
TARGET_CRS = "EPSG:4326"


def _load_raster(path: Path, clip_bounds_wgs84: Optional[Tuple[float, float, float, float]] = None) -> Tuple[np.ndarray, rasterio.Affine, Sequence[float]]:
    """Load a raster and reproject it to WGS84 (EPSG:4326)."""

    with rasterio.open(path) as src:
        if clip_bounds_wgs84 is not None:
            left, bottom, right, top = transform_bounds(TARGET_CRS, src.crs, *clip_bounds_wgs84, densify_pts=21)
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

        nodata = src.nodata
        if nodata is not None:
            data = np.where(data == nodata, np.nan, data)
        elif np.ma.isMaskedArray(data):
            data = data.filled(np.nan)

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
            src_nodata=nodata if nodata is not None else None,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )

    bounds = array_bounds(destination.shape[0], destination.shape[1], dst_transform)
    return destination, dst_transform, bounds


def _generate_rgba(
    data: np.ndarray,
    cmap_name: str,
    vmin: Optional[float],
    vmax: Optional[float],
    opacity: float,
) -> Tuple[np.ndarray, float, float]:
    """Convert the raster into an RGBA image and return color limits."""

    valid = np.isfinite(data)
    if not np.any(valid):
        raise RuntimeError("Raster nao possui pixels validos para renderizacao.")

    finite_values = data[valid]
    min_value = vmin if vmin is not None else float(np.nanpercentile(finite_values, 2))
    max_value = vmax if vmax is not None else float(np.nanpercentile(finite_values, 98))
    if np.isclose(min_value, max_value):
        max_value = min_value + 1e-3

    cmap = colormaps[cmap_name]
    normaliser = colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
    rgba = cmap(normaliser(data))
    rgba[..., 3] = np.where(valid, opacity, 0.0)
    image = (rgba * 255).astype(np.uint8)
    return image, min_value, max_value


def _load_geojson(geojson_path: Path) -> Dict[str, Any]:
    with geojson_path.open("r", encoding="utf-8") as file:
        return json.load(file)


def _iterate_geometries(geometry: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    gtype = geometry.get("type")
    if gtype == "FeatureCollection":
        for feature in geometry.get("features", []):
            yield from _iterate_geometries(feature)
    elif gtype == "Feature":
        yield from _iterate_geometries(geometry["geometry"])
    elif gtype == "Polygon":
        yield geometry
    elif gtype == "MultiPolygon":
        for polygon in geometry.get("coordinates", []):
            yield {"type": "Polygon", "coordinates": polygon}


def _add_geojson_layer(map_obj: folium.Map, geojson_data: Dict[str, Any]) -> None:
    folium.GeoJson(data=geojson_data, name="Area de interesse", style_function=lambda _: {"fillOpacity": 0}).add_to(map_obj)


def _extract_geometry_bounds(geojson_data: Dict[str, Any]) -> Optional[Tuple[float, float, float, float]]:
    def extract_geometry(geometry: Dict[str, Any]) -> Dict[str, Any]:
        if geometry.get("type") == "FeatureCollection":
            features = geometry.get("features", [])
            if not features:
                raise ValueError("GeoJSON sem features.")
            return extract_geometry(features[0])
        if geometry.get("type") == "Feature":
            return extract_geometry(geometry["geometry"])
        return geometry

    geometry = extract_geometry(geojson_data)
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
    min_lon, max_lon = min(lons), max(lons)
    min_lat, max_lat = min(lats), max(lats)
    return min_lon, min_lat, max_lon, max_lat


def _apply_unsharp_mask(data: np.ndarray, radius: float, amount: float) -> np.ndarray:
    mask = np.isfinite(data)
    if not np.any(mask):
        return data

    valid_values = data[mask]
    fill_value = float(np.nanmean(valid_values)) if valid_values.size else 0.0
    filled = np.where(mask, data, fill_value)
    blurred = gaussian_filter(filled, sigma=radius, mode="nearest")
    sharpened = filled + amount * (filled - blurred)
    return np.where(mask, sharpened, np.nan)


def _apply_smoothing(data: np.ndarray, radius: float) -> np.ndarray:
    if radius <= 0:
        return data
    mask = np.isfinite(data)
    if not np.any(mask):
        return data
    fill_value = float(np.nanmean(data[mask])) if mask.any() else 0.0
    filled = np.where(mask, data, fill_value)
    smoothed = gaussian_filter(filled, sigma=radius, mode="nearest")
    return np.where(mask, smoothed, np.nan)


def _upsample_raster(
    data: np.ndarray,
    transform: rasterio.Affine,
    scale_factor: float,
    nodata: float = np.nan,
) -> Tuple[np.ndarray, rasterio.Affine]:
    if scale_factor <= 1.0:
        return data, transform

    new_height = int(data.shape[0] * scale_factor)
    new_width = int(data.shape[1] * scale_factor)
    new_transform = transform * Affine.scale(1 / scale_factor, 1 / scale_factor)

    destination = np.full((new_height, new_width), np.nan, dtype=np.float32)
    reproject(
        source=data,
        destination=destination,
        src_transform=transform,
        src_crs=TARGET_CRS,
        dst_transform=new_transform,
        dst_crs=TARGET_CRS,
        src_nodata=nodata,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )
    return destination, new_transform


def build_map(
    index_path: Path,
    output_path: Path,
    cmap_name: str = DEFAULT_CMAP,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    opacity: float = 0.75,
    overlays: Optional[Iterable[Path]] = None,
    tiles: str = "CartoDB positron",
    tile_attr: Optional[str] = None,
    padding_factor: float = 0.3,
    upsample: float = 1.0,
    sharpen: bool = False,
    sharpen_radius: float = 1.0,
    sharpen_amount: float = 1.2,
    smooth_radius: float = 0.0,
) -> Path:

    overlay_geojsons = []
    clip_bounds = None
    if overlays:
        for geojson_path in overlays:
            overlay_geojsons.append(_load_geojson(geojson_path))

        geom_bounds = [_extract_geometry_bounds(geojson_data) for geojson_data in overlay_geojsons]
        geom_bounds = [b for b in geom_bounds if b is not None]
        if geom_bounds:
            min_lon_geo = min(b[0] for b in geom_bounds)
            min_lat_geo = min(b[1] for b in geom_bounds)
            max_lon_geo = max(b[2] for b in geom_bounds)
            max_lat_geo = max(b[3] for b in geom_bounds)
            width = max_lon_geo - min_lon_geo
            height = max_lat_geo - min_lat_geo
            pad_lon = width * padding_factor / 2
            pad_lat = height * padding_factor / 2
            clip_bounds = (
                min_lon_geo - pad_lon,
                min_lat_geo - pad_lat,
                max_lon_geo + pad_lon,
                max_lat_geo + pad_lat,
            )

    data, transform, bounds = _load_raster(index_path, clip_bounds_wgs84=clip_bounds)
    if sharpen:
        data = _apply_unsharp_mask(data, radius=sharpen_radius, amount=sharpen_amount)

    if overlay_geojsons:
        shapes = []
        for geojson_data in overlay_geojsons:
            for geom in _iterate_geometries(geojson_data):
                shapes.append((geom, 1))
        if shapes:
            mask = rasterize(
                shapes=shapes,
                out_shape=data.shape,
                transform=transform,
                fill=0,
                all_touched=False,
                dtype=np.uint8,
            )
            data = np.where(mask == 1, data, np.nan)

    data, transform = _upsample_raster(data, transform, upsample)
    data = _apply_smoothing(data, smooth_radius)

    if overlay_geojsons and upsample > 1.0:
        shapes = []
        for geojson_data in overlay_geojsons:
            for geom in _iterate_geometries(geojson_data):
                shapes.append((geom, 1))
        if shapes:
            mask = rasterize(
                shapes=shapes,
                out_shape=data.shape,
                transform=transform,
                fill=0,
                all_touched=False,
                dtype=np.uint8,
            )
            data = np.where(mask == 1, data, np.nan)

    image, min_value, max_value = _generate_rgba(data, cmap_name, vmin, vmax, opacity)

    min_lon, min_lat, max_lon, max_lat = bounds

    if clip_bounds is not None:
        min_lon, min_lat, max_lon, max_lat = clip_bounds

    centre_lat = (min_lat + max_lat) / 2
    centre_lon = (min_lon + max_lon) / 2

    if tiles.lower() == "none":
        base_map = folium.Map(location=[centre_lat, centre_lon], zoom_start=11, tiles=None)
    else:
        base_map = folium.Map(location=[centre_lat, centre_lon], zoom_start=11, tiles=tiles, attr=tile_attr)
        folium.TileLayer(
            tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
            attr="Esri World Imagery",
            name="Esri World Imagery",
            overlay=False,
            control=True,
        ).add_to(base_map)
    folium.raster_layers.ImageOverlay(
        image=image,
        bounds=[[min_lat, min_lon], [max_lat, max_lon]],
        opacity=1.0,
        name=index_path.stem,
    ).add_to(base_map)

    linear = LinearColormap(
        [colors.rgb2hex(colormaps[cmap_name](x)) for x in np.linspace(0, 1, 10)],
        vmin=min_value,
        vmax=max_value,
    )
    linear.caption = f"{index_path.stem} (min={min_value:.3f}, max={max_value:.3f})"
    linear.add_to(base_map)

    for geojson_data in overlay_geojsons:
        _add_geojson_layer(base_map, geojson_data)

    folium.LayerControl().add_to(base_map)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_map.save(str(output_path))
    return output_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True, help="Arquivo GeoTIFF do indice (NDVI/NDWI/MSI).")
    parser.add_argument("--output", type=Path, default=Path("mapa_indice.html"), help="Arquivo HTML de saida.")
    parser.add_argument("--colormap", default=DEFAULT_CMAP, help="Nome do colormap Matplotlib (default: RdYlGn).")
    parser.add_argument("--vmin", type=float, help="Valor minimo para o Stretch.")
    parser.add_argument("--vmax", type=float, help="Valor maximo para o Stretch.")
    parser.add_argument("--opacity", type=float, default=0.75, help="Opacidade da camada raster (0-1).")
    parser.add_argument("--geojson", type=Path, nargs="*", help="Arquivos GeoJSON adicionais para sobrepor no mapa.")
    parser.add_argument(
        "--tiles",
        default="CartoDB positron",
        help="Camada base do folium (ex.: 'OpenStreetMap'). Use 'none' para modo offline.",
    )
    parser.add_argument(
        "--tile-attr",
        default=None,
        help="Atribuicao personalizada para a camada base (quando aplicavel).",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=0.3,
        help="Fator de expansao do envelope do GeoJSON (default: 0.3).",
    )
    parser.add_argument("--sharpen", action="store_true", help="Aplica filtro de nitidez (unsharp mask).")
    parser.add_argument("--sharpen-radius", type=float, default=1.0, help="Raio (sigma) usado na gaussiana do filtro.")
    parser.add_argument("--sharpen-amount", type=float, default=1.3, help="Intensidade do realce de nitidez.")
    parser.add_argument("--upsample", type=float, default=1.0, help="Fator de upsample para suavizar pixels (ex.: 4).")
    parser.add_argument("--smooth-radius", type=float, default=0.0, help="Aplica suavizacao gaussiana apos o upsample.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    output = build_map(
        index_path=args.index.expanduser().resolve(),
        output_path=args.output.expanduser().resolve(),
        cmap_name=args.colormap,
        vmin=args.vmin,
        vmax=args.vmax,
        opacity=args.opacity,
        overlays=[path.expanduser().resolve() for path in (args.geojson or [])],
        tiles=args.tiles,
        tile_attr=args.tile_attr,
        padding_factor=args.padding,
        upsample=args.upsample,
        sharpen=args.sharpen,
        sharpen_radius=args.sharpen_radius,
        sharpen_amount=args.sharpen_amount,
        smooth_radius=args.smooth_radius,
    )
    print(f"Mapa salvo em {output}")


if __name__ == "__main__":
    main()

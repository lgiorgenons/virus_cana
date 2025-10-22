"""Renderiza um mapa interativo com composicao RGB (true color) a partir de bandas Sentinel-2.

As bandas B04 (vermelho), B03 (verde) e B02 (azul) devem estar disponiveis como GeoTIFFs,
como gerado pelo pipeline em `data/processed/<produto>/`. O resultado e um HTML com visualizacao
estatica (tiles opcionais) adequada para navegacao offline ou online.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Dict, Any

import folium
import numpy as np
import rasterio
from branca.colormap import LinearColormap
from matplotlib import colormaps
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio.warp import calculate_default_transform, reproject
from scipy.ndimage import gaussian_filter

TARGET_CRS = "EPSG:4326"


def _reproject_to_wgs84(
    band_path: Path,
    dst_transform: Optional[rasterio.Affine] = None,
    dst_width: Optional[int] = None,
    dst_height: Optional[int] = None,
) -> Tuple[np.ndarray, rasterio.Affine]:
    """Reprojeta uma banda para WGS84 utilizando os parametros de destino, se fornecidos."""

    with rasterio.open(band_path) as src:
        if dst_transform is None or dst_width is None or dst_height is None:
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src.crs, TARGET_CRS, src.width, src.height, *src.bounds
            )

        destination = np.empty((dst_height, dst_width), dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=TARGET_CRS,
            resampling=Resampling.bilinear,
        )

    return destination, dst_transform


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
    return min(lons), min(lats), max(lons), max(lats)


def _load_geojson(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as file:
        return json.load(file)


def _apply_unsharp_mask(array: np.ndarray, radius: float, amount: float) -> np.ndarray:
    mask = np.isfinite(array)
    if not np.any(mask):
        return array

    valid = array[mask]
    fill_value = float(np.nanmean(valid)) if valid.size else 0.0
    filled = np.where(mask, array, fill_value)
    blurred = gaussian_filter(filled, sigma=radius, mode="nearest")
    sharpened = filled + amount * (filled - blurred)
    return np.where(mask, sharpened, np.nan)


def _stretch(array: np.ndarray, lower: float = 2, upper: float = 98) -> np.ndarray:
    finite = np.isfinite(array)
    if not np.any(finite):
        raise RuntimeError("Banda sem valores validos.")

    vmin = np.percentile(array[finite], lower)
    vmax = np.percentile(array[finite], upper)
    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-3
    stretched = np.clip((array - vmin) / (vmax - vmin), 0, 1)
    stretched[~finite] = 0
    return stretched.astype(np.float32)


def _create_rgb_image(red: np.ndarray, green: np.ndarray, blue: np.ndarray) -> np.ndarray:
    r = _stretch(red)
    g = _stretch(green)
    b = _stretch(blue)
    rgb = np.stack([r, g, b], axis=-1)
    return (rgb * 255).astype(np.uint8)


def build_truecolor_map(
    red_path: Path,
    green_path: Path,
    blue_path: Path,
    output_path: Path,
    overlays: Optional[Iterable[Path]] = None,
    tiles: str = "CartoDB positron",
    tile_attr: Optional[str] = None,
    padding_factor: float = 0.3,
    sharpen: bool = False,
    sharpen_radius: float = 1.0,
    sharpen_amount: float = 1.2,
) -> Path:
    red_array, transform = _reproject_to_wgs84(red_path)
    height, width = red_array.shape
    green_array, _ = _reproject_to_wgs84(green_path, transform, width, height)
    blue_array, _ = _reproject_to_wgs84(blue_path, transform, width, height)

    if sharpen:
        red_array = _apply_unsharp_mask(red_array, sharpen_radius, sharpen_amount)
        green_array = _apply_unsharp_mask(green_array, sharpen_radius, sharpen_amount)
        blue_array = _apply_unsharp_mask(blue_array, sharpen_radius, sharpen_amount)

    bounds = array_bounds(height, width, transform)
    min_lon, min_lat, max_lon, max_lat = bounds

    geojson_data = [_load_geojson(path) for path in overlays] if overlays else []
    if geojson_data:
        geom_bounds = [_extract_geometry_bounds(data) for data in geojson_data]
        geom_bounds = [b for b in geom_bounds if b is not None]
        if geom_bounds:
            min_lon_geo = min(b[0] for b in geom_bounds)
            min_lat_geo = min(b[1] for b in geom_bounds)
            max_lon_geo = max(b[2] for b in geom_bounds)
            max_lat_geo = max(b[3] for b in geom_bounds)
            width_geo = max_lon_geo - min_lon_geo
            height_geo = max_lat_geo - min_lat_geo
            pad_lon = width_geo * padding_factor / 2
            pad_lat = height_geo * padding_factor / 2
            min_lon = min_lon_geo - pad_lon
            max_lon = max_lon_geo + pad_lon
            min_lat = min_lat_geo - pad_lat
            max_lat = max_lat_geo + pad_lat

    centre_lat = (min_lat + max_lat) / 2
    centre_lon = (min_lon + max_lon) / 2

    if tiles.lower() == "none":
        base_map = folium.Map(location=[centre_lat, centre_lon], zoom_start=12, tiles=None)
    else:
        base_map = folium.Map(location=[centre_lat, centre_lon], zoom_start=12, tiles=tiles, attr=tile_attr)
        folium.TileLayer(
            tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
            attr="Esri World Imagery",
            name="Esri World Imagery",
            overlay=False,
            control=True,
        ).add_to(base_map)

    image = _create_rgb_image(red_array, green_array, blue_array)
    folium.raster_layers.ImageOverlay(
        image=image,
        bounds=[[min_lat, min_lon], [max_lat, max_lon]],
        opacity=1.0,
        name="True color",
    ).add_to(base_map)

    if geojson_data:
        for data in geojson_data:
            folium.GeoJson(data=data, name="Area de interesse", style_function=lambda _: {"fillOpacity": 0}).add_to(base_map)

    # Opcional: legenda dummy para reforcar que se trata de composicao RGB
    colormap = LinearColormap(["#000000", "#FFFFFF"], vmin=0, vmax=255)
    colormap.caption = "Composicao RGB (stretch 2-98%)"
    colormap.add_to(base_map)

    folium.LayerControl().add_to(base_map)
    base_map.fit_bounds([[min_lat, min_lon], [max_lat, max_lon]])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    base_map.save(str(output_path))
    return output_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--red", type=Path, required=True, help="Caminho para banda vermelha (B04).")
    parser.add_argument("--green", type=Path, required=True, help="Caminho para banda verde (B03).")
    parser.add_argument("--blue", type=Path, required=True, help="Caminho para banda azul (B02).")
    parser.add_argument("--output", type=Path, default=Path("mapas/truecolor.html"), help="Arquivo HTML de saida.")
    parser.add_argument("--geojson", type=Path, nargs="*", help="Arquivos GeoJSON para sobrepor no mapa.")
    parser.add_argument("--tiles", default="CartoDB positron", help="Camada base do folium (ou 'none').")
    parser.add_argument("--tile-attr", default=None, help="Atribuicao para a camada base.")
    parser.add_argument("--padding", type=float, default=0.3, help="Fator de expansao do envelope do GeoJSON.")
    parser.add_argument("--sharpen", action="store_true", help="Aplica unsharp mask para realcar detalhes.")
    parser.add_argument("--sharpen-radius", type=float, default=1.0, help="Raio (sigma) da gaussiana na nitidez.")
    parser.add_argument("--sharpen-amount", type=float, default=1.2, help="Intensidade do realce de nitidez.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    overlays = [path.expanduser().resolve() for path in (args.geojson or [])]
    output = build_truecolor_map(
        red_path=args.red.expanduser().resolve(),
        green_path=args.green.expanduser().resolve(),
        blue_path=args.blue.expanduser().resolve(),
        output_path=args.output.expanduser().resolve(),
        overlays=overlays,
        tiles=args.tiles,
        tile_attr=args.tile_attr,
        padding_factor=args.padding,
        sharpen=args.sharpen,
        sharpen_radius=args.sharpen_radius,
        sharpen_amount=args.sharpen_amount,
    )
    print(f"Mapa true color salvo em {output}")


if __name__ == "__main__":
    main()

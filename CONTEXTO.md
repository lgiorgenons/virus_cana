# Contexto do Projeto

## O que foi entregue
- Script `scripts/satellite_pipeline.py` atualizado para consumir o catalogo OData do Copernicus Data Space, autenticar via OAuth2 (senha + client id), reutilizar produtos SAFE locais e selecionar dinamicamente os indices espectrais disponiveis.
- Extração padronizada das bandas visíveis, red-edge, NIR, SWIR e bandas atmosféricas (B01, B09, B10) e cálculo automático de índices (NDVI, NDWI, MSI, EVI, NDRE, NDMI, NDRE1-4, CI_REDEDGE, SIPI) reaproveitando reprojeções em disco.
- Scripts `scripts/render_index_map.py` e `scripts/render_truecolor_map.py` para gerar mapas HTML (offline/online) com nitidez ajustável e camada de fundo Esri World Imagery opcional.
- README detalhando setup, autenticação, execução do pipeline, visualização (NDVI, EVI, NDRE, NDMI, true color) e solução de problemas.

## Fluxo atual
1. Opcionalmente reutiliza um SAFE existente (`--safe-path`) ou baixa o produto mais recente com base em poligono, intervalo de datas e cobertura de nuvens.
2. Extrai bandas relevantes (visível, red-edge, NIR, SWIR) e as salva como GeoTIFF reutilizáveis.
3. Calcula os índices solicitados (`--indices`) e grava os resultados em `<workdir>/<produto>/indices/<indice>.tif` (NDVI, NDWI, MSI, EVI, NDRE, NDMI, NDRE1-4, CI_REDEDGE, SIPI).
4. Gera mapas HTML interativos (NDVI, EVI, NDRE, NDMI, true color) com nitidez ajustável; quando online, oferece a camada Esri World Imagery como fundo.

## Pendencias e proximos passos
- Estender o conjunto de indicadores (ex.: SIPI, NDVIre, MCARI2) e estatisticas temporais para diagnósticos precoces.
- Criar um fluxo automatizado (cron/notebook) para baixar novos produtos periodicamente e comparar a evolucao temporal dos indices.
- Incorporar analises agronomicas adicionais (relatorios graficos ou alertas) usando os rasters gerados como entrada.
- Definir estrategia de armazenamento para os SAFE baixados (limpeza periodica ou upload para bucket dedicado).
- Integrar sensores complementares (Landsat/ECOSTRESS, SMAP/CHIRPS, Sentinel-1) e basemaps de maior resolucao quando houver conectividade, mantendo claro o limite de 10 m do Sentinel-2.

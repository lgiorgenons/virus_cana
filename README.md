# Sentinel-2 Sugarcane Health Toolkit

Python utilities to download Sentinel-2 Level-2A imagery and generate vegetation stress indicators that support sugarcane wilt monitoring. The workflow follows the recommendations documented in `Analise Detalhada e Abrangente da Sindrome da Murcha da Cana-de-Acucar.md`.

## Requisitos
- Python 3.10 ou superior
- Conta ativa no [Copernicus Open Access Hub](https://scihub.copernicus.eu/)
- Utilitarios do sistema: `unzip` **nao** e necessario (a extracao acontece via `zipfile`)

Instale as dependencias em um ambiente virtual:

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Autenticacao
Por padrao o script le as credenciais das variaveis de ambiente:

- `SENTINEL_USERNAME`
- `SENTINEL_PASSWORD`
- `SENTINEL_API_URL` (opcional, padrao `https://catalogue.dataspace.copernicus.eu/odata/v1/`)
- `SENTINEL_TOKEN_URL` (opcional, padrao `https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token`)
- `SENTINEL_CLIENT_ID` (opcional, padrao `cdse-public`)

Tambem e possivel fornece-las pela linha de comando via `--username`, `--password`, `--api-url`, `--token-url` e `--client-id`, o que substitui o valor das variaveis. O script solicita um token OAuth2 automaticamente antes de consultar o catalogo OData e reutiliza a mesma sessao para o download.

### Habilitando o acesso OData
Antes da primeira execucao, ative a API OData para o seu usuario no portal do Copernicus Data Space:

1. Acesse `https://dataspace.copernicus.eu/`, clique em **Login → Account Console** e entre com suas credenciais.
2. Em **Personal info**, marque os termos de uso/privacidade e salve.
3. Abra **Applications** e habilite o aplicativo **OData API** (status *In use*).
4. Opcional: gere um client secret apenas se for usar um cliente proprio; para este script o `client_id` padrao `cdse-public` ja funciona.
5. Aguarde alguns minutos, faca logoff/logon e teste:
   ```bash
   curl -I -u "$SENTINEL_USERNAME:$SENTINEL_PASSWORD" \
     'https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$top=1'
   ```
   A chamada deve retornar `HTTP/1.1 200 OK`. Se receber `401/403`, repita o passo 3 ou acione o suporte.

## Uso rapido
Baixe o produto mais recente para um poligono GeoJSON e gere todos os indices (NDVI, NDWI, MSI):

```bash
python scripts/satellite_pipeline.py \
  --geojson dados/aoi.geojson \
  --start-date 2024-04-01 \
  --end-date 2024-04-30 \
  --cloud 0 20 \
  --download-dir data/raw \
  --workdir data/processed \
  --log-level INFO
```

Parametros importantes:

- `--cloud` controla a faixa aceitavel de cobertura de nuvens (percentual minimo e maximo).
- `--indices` aceita uma lista (`ndvi`, `ndwi`, `msi`) caso deseje gerar apenas alguns produtos.
- `--safe-path` permite pular o download e reutilizar um arquivo `.SAFE` (zip) ou diretorio ja existente.

### Exemplo (reaproveitando download)
```bash
python scripts/satellite_pipeline.py \
  --safe-path data/raw/S2A_MSIL2A_*.SAFE.zip \
  --workdir data/processed \
  --indices ndvi ndwi msi evi ndre ndmi ndre1 ndre2 ndre3 ndre4 ci_rededge sipi
```

## Saidas
Os arquivos GeoTIFF dos indices sao salvos em:

```
<workdir>/<produto>/indices/<indice>.tif
```

Os nomes dos indices:

- `ndvi`: vigor vegetativo.
- `ndwi`: status hidrico.
- `msi`: indicador de estresse por falta de umidade.
- `evi`: vigor em dosses densos com correção atmosférica.
- `ndre`: sensível à clorofila (banda red-edge).
- `ndmi`: umidade da vegetação usando NIR e SWIR.
- `ndre1`/`ndre2`/`ndre3`/`ndre4`: variantes NDRE usando cada banda red-edge (B05-B08A) para detectar clorofila em diferentes profundidades do dossel.
- `ci_rededge`: índice de clorofila baseado na razão `nir/rededge4 - 1`.
- `sipi`: Structure Insensitive Pigment Index, bom para avaliar relações clorofila/carotenoides.

### Visualizacao rapida dos indices
Gere um mapa interativo (HTML) sobrepondo um indice ao mapa de fundo:

```bash
python scripts/render_index_map.py \
  --index data/processed/<produto>/indices/ndvi.tif \
  --geojson dados/map.geojson \
  --upsample 12 --smooth-radius 1 \
  --sharpen --sharpen-radius 1.2 --sharpen-amount 1.5 \
  --output mapas/ndvi.html
```

Substitua `ndvi.tif` por `evi.tif`, `ndre.tif` ou `ndmi.tif` para visualizar os demais indices derivados do Sentinel-2.

A visualizacao usa, por padrao, o mapa `CartoDB positron` (necessita internet). Em ambiente offline execute com `--tiles none` para carregar apenas o raster.

A pagina resultante pode ser aberta diretamente no navegador. Ajuste `--colormap`, `--vmin`, `--vmax` e os parametros de nitidez (`--sharpen-radius`, `--sharpen-amount`) conforme necessario.
Use `--upsample` (ex.: 10-12) e `--smooth-radius` para suavizar pixels e aproximar a escala visual do basemap.

Para uma composicao RGB (true color) com nitidez maxima:

```bash
python scripts/render_truecolor_map.py \
  --red data/processed/<produto>/red.tif \
  --green data/processed/<produto>/green.tif \
  --blue data/processed/<produto>/blue.tif \
  --geojson dados/map.geojson \
  --sharpen \
  --tiles OpenStreetMap \
  --output mapas/truecolor.html
```

> 💡 **Limite de resolucao:** os produtos Sentinel-2 oferecem pixels de 10 m; os filtros de nitidez ajudam a realcar contrastes, mas nao criam detalhes alem do que o sensor realmente capturou. Para imagens sub-metricas, utilize basemaps externos (Esri, Mapbox, etc.) ou dados comerciais.

Para alternar entre uso offline e um basemap de alta resolucao (Esri World Imagery), utilize o seletor do mapa gerado. Se quiser forcar modo offline, use `--tiles none`; para uma visualizacao mais nitida online, use qualquer camada suportada (`--tiles OpenStreetMap`) e ative a camada “Esri World Imagery” diretamente no mapa. (O script adiciona automaticamente a camada Esri quando um tileset online e usado.)

## Estrutura sugerida
```
data/
- raw/           # Produtos SAFE baixados
- processed/     # Bandas extraidas e indices calculados

As bandas salvas incluem `coastal`, `blue`, `green`, `red`, `rededge1-4`, `nir`, `water_vapor`, `cirrus`, `swir1` e `swir2`, permitindo calcular indices baseados nas bandas red-edge e SWIR do Sentinel-2.
```

## Solucao de problemas
- **Credenciais ausentes**: confirme variaveis de ambiente ou passe `--username`/`--password`.
- **Erro 401/403 ao consultar OData**: confira se o aplicativo *OData API* foi habilitado no portal do Copernicus Data Space e se o usuario aceitou os termos de uso; apos a ativacao pode levar alguns minutos para propagar.
- **Erro DAT-ZIP-111 / 422 ao baixar**: normalmente indica que o endpoint de download recebeu um UUID em formato incorreto; confirme que esta usando `Products(<uuid>)/$value` (sem aspas) ou simplesmente deixe o script cuidar do download.
- **Nenhum produto encontrado**: ajuste intervalo de datas, cobertura de nuvens ou geometria do poligono.
- **Bandas ausentes**: verifique se o SAFE contem as bandas B04, B08 e B11; algumas cenas incompletas nao incluem todas.
- **Erros de CRS/resolucao**: o script reprojeta bandas automaticamente para a resolucao da banda NIR (10 m). Use um SAFE integro para evitar incompatibilidades.

## Proximos passos
- Automatizar o pipeline (cron/notebook) para baixar cenas periodicamente e gerar series temporais (NDVI, EVI, NDRE, NDMI).
- Enriquecer os indicadores com novos indices (SIPI, NDVIre, MCARI2) e estatisticas de anomalia por talhao.
- Integrar sensores adicionais (Landsat/ECOSTRESS para temperatura, SMAP/CHIRPS para umidade/chuva, Sentinel-1 para estrutural) e combinar com o Sentinel-2.
- Definir governanca de armazenamento (rotina de limpeza dos SAFEs ou upload para bucket dedicado) e entrega de relatorios/alertas agronomicos automaticamente.

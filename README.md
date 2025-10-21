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
- `SENTINEL_API_URL` (opcional, padrao `https://apihub.copernicus.eu/apihub`)

Tambem e possivel fornece-las pela linha de comando via `--username`, `--password` e `--api-url`, o que substitui o valor das variaveis.

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
  --indices ndvi ndwi
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

## Estrutura sugerida
```
data/
- raw/           # Produtos SAFE baixados
- processed/     # Bandas extraidas e indices calculados
```

## Solucao de problemas
- **Credenciais ausentes**: confirme variaveis de ambiente ou passe `--username`/`--password`.
- **Nenhum produto encontrado**: ajuste intervalo de datas, cobertura de nuvens ou geometria do poligono.
- **Bandas ausentes**: verifique se o SAFE contem as bandas B04, B08 e B11; algumas cenas incompletas nao incluem todas.
- **Erros de CRS/resolucao**: o script reprojeta bandas automaticamente para a resolucao da banda NIR (10 m). Use um SAFE integro para evitar incompatibilidades.

## Proximos passos
- Automatizar a execucao via cron/notebook integrando as saidas com analises agronomicas.
- Adicionar metricas complementares (ex.: NDRE) conforme novas bandas sejam necessarias.

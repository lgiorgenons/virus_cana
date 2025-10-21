# Contexto do Projeto

## O que foi entregue
- Script `scripts/satellite_pipeline.py` atualizado para aceitar credenciais via CLI, reutilizar produtos SAFE locais e selecionar dinamicamente os indices espectrais disponiveis.
- Registro padronizado de indices (NDVI, NDWI, MSI) com requisitos de bandas e reaproveitamento de dados reprojetados.
- README com instrucoes de setup, autenticacao, exemplos de uso on-line/off-line e solucao de problemas.

## Fluxo atual
1. Opcionalmente reutiliza um SAFE existente (`--safe-path`) ou baixa o produto mais recente com base em poligono, intervalo de datas e cobertura de nuvens.
2. Extrai bandas relevantes para o monitoramento da murcha (B04, B08, B11, B12) e as salva como GeoTIFF.
3. Calcula os indices solicitados (`--indices`) e grava os resultados em `<workdir>/<produto>/indices/<indice>.tif`.

## Pendencias e proximos passos
- Validar o pipeline end-to-end com credenciais reais e um poligono de teste, confirmando a geracao correta dos arquivos.
- Revisar se ha indices adicionais de interesse (ex.: NDRE, SIPI) e adiciona-los ao dicionario `INDEX_SPECS`.
- Criar um fluxo automatizado (cron/notebook) para baixar novos produtos periodicamente e comparar a evolucao temporal dos indices.
- Incorporar analises agronomicas adicionais (relatorios graficos ou alertas) usando os rasters gerados como entrada.

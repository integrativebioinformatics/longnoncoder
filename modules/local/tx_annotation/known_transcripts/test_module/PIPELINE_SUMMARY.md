# 🧬 Pipeline de Teste - Módulo de Anotação de Transcritos

## ✅ Pipeline Completo Criado!

Criei um pipeline completo e funcional para testar seu módulo de anotação de transcritos.

## 📦 O que foi incluído:

### 🔧 **Arquivos Principais**
- **`main.nf`** - Processo Nextflow com container e conda
- **`known_transcripts.R`** - Script R com interface CLI
- **`test_workflow.nf`** - Workflow de teste
- **`nextflow.config`** - Configuração completa

### 📊 **Dados de Teste Realistas**
- **10 transcritos Ensembl** com diferentes biotipos
- **Contagens simuladas** para 3 amostras
- **Arquivo GTF** com anotações correspondentes

### 🚀 **Scripts de Execução**
- **`run_test.sh`** - Teste completo com validação
- **`quick_test.sh`** - Teste rápido em modo stub

### 📚 **Documentação**
- **`README.md`** - Guia completo de uso
- **`PIPELINE_SUMMARY.md`** - Este resumo

## 🎯 **Como Usar**

### Teste Rápido (Recomendado primeiro)
```bash
cd test_pipeline
./quick_test.sh
```

### Teste Completo
```bash
cd test_pipeline
./run_test.sh
```

### Teste Manual
```bash
cd test_pipeline
nextflow run test_workflow.nf -profile conda
```

## 📋 **Outputs Esperados**

O pipeline gerará **11 arquivos**:
- 5 arquivos CSV de metadados
- 3 arquivos GTF filtrados  
- 2 arquivos de contagens filtradas
- 1 arquivo de versões

## 🔍 **Validação**

O script de teste automaticamente:
- ✅ Executa o pipeline
- ✅ Verifica se todos os arquivos foram gerados
- ✅ Conta linhas nos arquivos principais
- ✅ Reporta o status de cada output

## 🌟 **Características do Pipeline**

- **Flexível**: Suporta conda, docker e singularity
- **Robusto**: Tratamento de erros e casos extremos
- **Documentado**: README completo e comentários
- **Testável**: Dados de exemplo e scripts de teste
- **Configurável**: Parâmetros ajustáveis via config

## 🚀 **Pronto para Usar!**

Seu módulo está agora:
- ✅ Completamente funcional
- ✅ Bem documentado
- ✅ Pronto para testes
- ✅ Integração em pipelines maiores

Execute `./quick_test.sh` para começar!
# Pipeline de Teste - Módulo de Anotação de Transcritos

Este é um pipeline completo para testar o módulo de anotação de transcritos que você desenvolveu.

## 📁 Estrutura do Pipeline

```
test_pipeline/
├── main.nf                           # Processo principal do módulo
├── test_workflow.nf                  # Workflow de teste
├── known_transcripts.R          # Script R com interface CLI
├── nextflow.config                  # Configuração do pipeline
├── run_test.sh                      # Script para executar o teste
├── README.md                        # Esta documentação
└── test_data/                       # Dados de teste
    ├── BambuOutput_counts_transcript.txt
    ├── BambuOutput_counts_gene.txt
    └── BambuOutput_extended_annotations.gtf
```

## 🚀 Como Executar o Teste

### Opção 1: Script Automático (Recomendado)
```bash
cd test_pipeline
./run_test.sh
```

### Opção 2: Comando Manual
```bash
cd test_pipeline
nextflow run test_workflow.nf -profile conda
```

### Opção 3: Com Docker/Singularity
```bash
cd test_pipeline
nextflow run test_workflow.nf -profile docker
# ou
nextflow run test_workflow.nf -profile singularity
```

## 📊 Dados de Teste

Os dados de teste incluem:

- **10 transcritos Ensembl** (ENST IDs) com diferentes biotipos:
  - lncRNAs (ENST00000473358.1, ENST00000469289.1)
  - Protein-coding (ENST00000488147.1, ENST00000518655.2)
  - Outros tipos (miRNA, pseudogenes)

- **Contagens simuladas** para 3 amostras (sample1, sample2, sample3)

- **Arquivo GTF** com anotações correspondentes

## 📋 Outputs Esperados

O módulo deve gerar os seguintes arquivos:

### Metadados (CSV)
- `annotated_transcriptome_metadata.csv` - Metadados completos do transcriptoma
- `annotated_lncRNAs_metadata.csv` - Metadados específicos de lncRNAs
- `annotated_lncRNAs_exonlength.csv` - Comprimentos de exons de lncRNAs
- `annotated_protein-coding_metadata.csv` - Metadados de transcritos codificantes
- `annotated_protein-coding_exonlength.csv` - Comprimentos de exons codificantes

### Contagens Filtradas (CSV)
- `bambu_annotated_transcriptome_tx_counts.csv` - Contagens de transcritos anotados
- `bambu_annotated_transcriptome_gene_counts.csv` - Contagens de genes anotados

### Arquivos GTF
- `bambu_annotated_transcriptome.gtf` - Transcriptoma anotado completo
- `bambu_annotated_lncRNAs.gtf` - Apenas lncRNAs
- `bambu_annotated_mRNAs.gtf` - Apenas transcritos codificantes

### Versões
- `versions.yml` - Versões dos softwares utilizados

## ⚙️ Configurações

### Perfis Disponíveis
- `conda` - Usa ambientes conda (padrão)
- `docker` - Usa containers Docker
- `singularity` - Usa containers Singularity
- `test` - Recursos reduzidos para teste rápido

### Parâmetros Configuráveis
```bash
# Versão do Ensembl (padrão: 113)
--ensembl_version 112

# Diretório de saída (padrão: ./results)
--outdir /path/to/output
```

## 🔧 Troubleshooting

### Problemas Comuns

1. **Erro de conexão biomaRt**
   - Verifique sua conexão com a internet
   - O biomaRt pode estar temporariamente indisponível

2. **Falta de pacotes R**
   - Certifique-se de que o conda está instalado
   - Use o perfil correto (`-profile conda`)

3. **Memória insuficiente**
   - Use o perfil test: `-profile test`
   - Ou ajuste a memória no nextflow.config

### Logs e Debugging
```bash
# Ver logs detalhados
nextflow run test_workflow.nf -profile conda -with-trace -with-report

# Modo debug
nextflow run test_workflow.nf -profile conda --debug
```

## 📈 Validação dos Resultados

Após a execução, verifique:

1. **Todos os arquivos foram gerados** (10 arquivos + versions.yml)
2. **Arquivos não estão vazios** (exceto se não houver dados do tipo específico)
3. **Metadados contêm informações do biomaRt** (chromosome_name, gene_biotype, etc.)
4. **GTFs contêm apenas transcritos anotados** (com IDs Ensembl)

## 🎯 Próximos Passos

Após validar que o módulo funciona:

1. **Integre em seu pipeline principal**
2. **Ajuste parâmetros** conforme necessário
3. **Teste com dados reais** do seu projeto
4. **Configure recursos** adequados para seus dados

## 📞 Suporte

Se encontrar problemas:
1. Verifique os logs do Nextflow (`.nextflow.log`)
2. Confirme que todos os arquivos de entrada existem
3. Teste com o perfil `test` para recursos reduzidos
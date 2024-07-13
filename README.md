# The transcriptional landscape of PPRE responsive genes in right ventricular heart nuclei between control and PKP2-cardiomyopathy
## PPRE responsive genes
### Datasets
These datasets are used to establish a PPRE responsive geneset. Specific datasets on PPARs and RXRs, specific to heart and liver were found. 
| Data type  | Target          | Database                | Source ID   | Publication                                                       |
|------------|-----------------|-------------------------|-------------|-------------------------------------------------------------------|
| RNA-seq    | RNA             | GEO Expression Omnibus  | GSE160987   | Wickramasinghe NM et al., 2022 PMID: 35325615                     |
| ATAC-seq   | Open chromatin  |                         | GSE178984   | Wickramasinghe NM et al., 2022 PMID: 35325615                     |
| ChIP-seq   | PPARα           |                         | GSE244905   | Xie et al., 2023 PMID: 37994595                                   |
|            | PPARα           |                         | GSM4748812  | Mishra S et al., 2021 PMID: 34618683                              |
|            | PPARα           |                         | GSM4748813  | Mishra S et al., 2021 PMID: 34618683                              |
|            | PPARγ           |                         | GSM624141   | Pott S et al., 2012 PMID: 23118933                                |
|            | RXR             |                         | GSM624142   | Pott S et al., 2012 PMID: 23118933                                |
|            | RXR             |                         | GSE188994   | Paredes A et al., 2023 PMID: 37225978                             |
|            | H3K27ac         | In house (Cardiomnics Lab) | -           | https://www.biorxiv.org/content/10.1101/2020.11.30.402792v1.full  |
| Motif      | PPARα-RXRα      | JASPAR                  | MA11481.1   | Rauluseviciute I et al., 2023 doi:10.1093/nar/gkad1059            |
|            |                 |                         | MA00651.1   | Rauluseviciute I et al., 2023 doi:10.1093/nar/gkad1059            |
| GO         | Proteins        | AmiGO                   | GO:0006635  |                                                                   |

### Analysis
Collaboration with Tim van der Wiel

## snRNAseq analysis
### Retrieved data
snRNAseq: [DCM heart atlas](https://github.com/heiniglab/DCM_heart_cell_atlas) (Reichert et al., 2022)

### Analysis
Reanalysis is performed using Seurat v5



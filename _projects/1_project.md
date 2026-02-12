---
layout: page
title: 10X multiome analysis in R
description: with background image
img: assets/img/Common_marmoset_(Callithrix_jacchus).jpeg
importance: 1
category: workbench
related_publications: true
---

The following workflow is an analysis of publicly available raw data from the Allen Institute's [Brain Knowledge Platform](https://brain-map.org/bkp). These analyses were carried out on a laptop with 8 cores and 64GB of memory, so to address computational resource and storage limitations, I analyzed 10X Multiomic data generated from brainstem tissue in the common marmoset *Callithrix jacchus*.

In addition to the [Seurat](https://satijalab.org/seurat/articles/get_started_v5_new) and [Signac](https://stuartlab.org/signac/articles/overview) documentation, I also found the following resources incredibly helpful: \
<https://ngs101.com/tutorials/> \
<https://github.com/mousepixels/sanbomics_scripts>

## Contents

*   [System setup](#system-setup)
*   [Preprocessing](#preprocessing)
*   [Normalization](#normalization)
*   [Celltype annotation with MapMyCells](#cell-annotation)
*   [Differential expression analysis](#differential-expression-analysis)
*   [Differential accessibility analysis](#differential-accessibility-analysis)
*   [Addendum: Celltype prediction with SingleR](#addendum-celltype-prediction-with-primate-data-and-singler)

## System setup

**Note**: At the time of this exercise the developer version of Seurat was required to fix a [bug](https://github.com/satijalab/seurat/issues/10180) encountered in downstream plotting.

~~~
# analysis
library(Seurat) # v5.3.1.9001
library(Signac)
library(presto) # much faster Wilcoxon rank sum test
library(SingleR)

# object manipulation
library(anndataR)
library(scCustomize)
library(reticulate)

# annotation
library(rtracklayer)
library(GenomeInfoDb)
library(biomaRt)

# parallelization
library(future)
plan(multisession, workers = 2) # parallelization
options(future.globals.maxSize= 11000*1024^2) # e.g. 6000gb*1024^2

# visualization
library(ggplot2)
library(pheatmap)
library(VennDiagram)
~~~

## Preprocessing

Loading metadata for the raw counts:
```
meta <- read.csv(
  file = '2021-08-25_Marm038_brainstem_sl9_rxn1.atac.per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
```

For the custom annotation, I had to use the same _C. jacchus_ genome build as what the reads were aligned to (NCBI RefSeq Assembly ID: [GCF_009663435.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009663435.1/) ('caljac4' in the UCSC browser).
```
annotations <- rtracklayer::import('Callithrix_jacchus_cj1700_1/ncbi_dataset/data/GCF_009663435.1/genomic.gtf')
```

Some functions to be used downstream expect certain column names in the Assay objects and will call errors without them, so they're added here.
```
colnames(mcols(annotations))[colnames(mcols(annotations)) == "gene"] <- "gene_name" # gene_name required for CreateChromatinAssay() annotation parameter
annotations$tx_id <- annotations$transcript_id # tx_id column potentially required to prevent CoveragePlot() error
annotations$gene_name <- annotations$gene # gene_name col required for CreateChromatinAssay() annotation argument
colnames(mcols(annotations))
```

Seurat's `CreateChromatinAssay()` also requires UCSC format for the chromosome labels of the reference (currently RefSeq), so here they're changed by vector mapping with the Alias.txt downloaded from the UCSC Genome Browser.
```
alias_map_df <- read.csv("Callithrix_jacchus_cj1700_1/GCF_009663435.1.chromAlias.txt", 
                         header = TRUE, sep = "\t"
)
head(alias_map_df)
colnames(alias_map_df)[1] <- "refseq"
head(alias_map_df)
#       refseq assembly    genbank ncbi ucsc
#1 NC_025586.1       MT              MT chrM
#2 NC_048383.1     chr1 CM018917.1    1 chr1
#3 NC_048384.1     chr2 CM018918.1    2 chr2
#4 NC_048385.1     chr3 CM018919.1    3 chr3
#5 NC_048386.1     chr4 CM018920.1    4 chr4
#6 NC_048387.1     chr5 CM018921.1    5 chr5

alias_map_vector <- setNames(alias_map_df$ucsc, alias_map_df$refseq)
alias_map_vector <- alias_map_vector[!is.na(alias_map_vector) & alias_map_vector != ""] # remove NAs and blanks in case present
head(alias_map_vector)
#NC_025586.1 NC_048383.1 NC_048384.1 NC_048385.1 NC_048386.1 NC_048387.1 
#     "chrM"      "chr1"      "chr2"      "chr3"      "chr4"      "chr5" 
```

Next, a vector of the current chromosomes/scaffolds in the genome that need to be changed is constructed, mapped to the names to names to UCSC names, and assigned to the genome/annotation object.
```
current_seqlevels <- seqlevels(annotations)
new_seqlevels <- alias_map_vector[as.character(current_seqlevels)]
seqlevels(annotations) <- new_seqlevels
```

The presence of NA's in the 'gene_biotype' column of the annotation object will call an error, so rows with them are removed.
```
annotations <- annotations[!is.na(annotations$gene_biotype)] # NAs will introduce error during TSSEnrichment
seqlevelsStyle(annotations) # confirmed that new style is UCSC
```

## Create multiome Seurat object

The raw RNA and ATAC feature matrices were individually loaded and filtered to remove all droplets deemed cell non-calls by the CellRanger pipeline.
```
counts <- Read10X_h5(filename = "2021-08-25_Marm038_brainstem_sl9_rxn1.atac.raw_feature_bc_matrix.h5")

multiome_data <- CreateSeuratObject(counts = counts$`Gene Expression`, meta.data = meta, assay = "RNA")
multiome_data <- subset(multiome_data, subset = is_cell == 1)

atac <- CreateSeuratObject(counts = counts$Peaks, meta.data = meta, assay = "ATAC")
atac <- subset(atac, subset = is_cell == 1)
```

Chromatin assay object creation with the ATAC feature matrix requires both the fragments and the annotation object. Raw objects are no longer needed and can be removed to free up memory.
```
multiome_data[["ATAC"]] <- CreateChromatinAssay(
  counts = LayerData(atac), 
  sep = c(":", "-"),
  genome = "calJac4",
  annotation = annotations,
  fragments = './2021-08-25_Marm038_brainstem_sl9_rxn1.atac_fragments.tsv.gz',
)
rm(atac, counts, current_seqlevels, new_seqlevels)
gc()
```

Per-cell quality control metrics were calculated and visualized to inform additional filtering.
```
DefaultAssay(multiome_data) <- "ATAC"
multiome_data <- NucleosomeSignal(object = multiome_data)
multiome_data <- TSSEnrichment(object = multiome_data, fast = FALSE)

DensityScatter(multiome_data, x="nCount_ATAC", y="TSS.enrichment", 
               log_x=TRUE, quantiles=TRUE)
VlnPlot(
  object = multiome_data,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = FALSE,
  ncol = 3
)
```

Here I filter cells with cuttoffs for counts, TSS enrichment, and nucleosome signal based on the bottom and top quantiles shown in the density scatterplot.

```
multiome_data_filtered <- subset(
  x = multiome_data,
  subset = nCount_ATAC < 74000 &
    nCount_ATAC > 2000 &
    TSS.enrichment > 2.7 &
    nucleosome_signal < 1.75
)
rm(multiome_data)
gc()
```

## Normalization

**Note**: MapMyCells does not require normalization for cell annotation, and developers urge caution in interpreting results from log-normalized data, so if performing this step first, it's important to reference the correct layer within the Seurat object when formatting its input.

```
DefaultAssay(multiome_data_filtered) <- "RNA"
multiome_data_filtered <- NormalizeData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- FindVariableFeatures(multiome_data_filtered, assay="RNA", nfeatures = 5000)
multiome_data_filtered <- ScaleData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- RunPCA(multiome_data_filtered)

DefaultAssay(multiome_data_filtered) <- "ATAC"
multiome_data_filtered <- FindTopFeatures(multiome_data_filtered, min.cutoff = 5)
multiome_data_filtered <- RunTFIDF(multiome_data_filtered)
multiome_data_filtered <- RunSVD(multiome_data_filtered)
```

## Cell annotation

Fomatting an anndata (.h5ad) or .csv file for MapMyCells input:
```
as.anndata(
  multiome_data_filtered,
  file_path = "",
  file_name= "multiome_data_filtered.h5ad",
  assay = "RNA",
  main_layer = "counts",
  other_layers = NULL,
  transer_dimreduc = FALSE,
  verbose = TRUE,
)

SaveSeuratRds(object = multiome_data_filtered, file = "multiome_data_filtered.Rds")
```

Loading results for cell annotation using the [MapMyCells]("https://brain-map.org/bkp/analyze/mapmycells") browser platform with the 10x Whole Mouse Brain Taxonomy (CCN20230722) as the reference and adding to the multiome object. Seurat rownames were renamed as the cell_ID for integration with the Seurat object.
```
mapmycells_WMB <- read.csv(file = "multiome_data_filtered_10xWholeMouseBrain(CCN20230722)/multiome_data_filtered_10xWholeMouseBrain(CCN20230722).csv",
                             skip = 4, header = TRUE)https://brain-map.org/bkp/analyze/mapmycells
colnames(mapmycells_WMB)
# [1] "cell_id"                             "class_label"                        
# [3] "class_name"                          "class_bootstrapping_probability"    
# [5] "subclass_label"                      "subclass_name"                      
# [7] "subclass_bootstrapping_probability"  "supertype_label"                    
# [9] "supertype_name"                      "supertype_bootstrapping_probability"
#[11] "cluster_label"                       "cluster_name"                       
#[13] "cluster_alias"                       "cluster_bootstrapping_probability" 

rownames(mapmycells_WMB) <- mapmycells_WMB$cell_id
```

From the distribution of celltype counts (N=16227), over 60% of the cells are identified as oligodendrocytes or oligodendrocyte precursor cells ("31 OPC-Oligo"), with the 2nd and 3rd largest clusters consisting of glutamatergic neurons ("23 P Glut") and astrocytes or ependymal cells ("30 Astro-Epen"), respectively.
```
data.frame(Cell_counts=head(sort(table(mapmycells_WMB$class_name),decreasing=T),34))
# most well represented cells
# Cell_counts.Var1 Cell_counts.Freq
# 1       31 OPC-Oligo            10258
# 2          23 P Glut             2550
# 3      30 Astro-Epen             1407
# 4          34 Immune              861
# 5         20 MB GABA              269
# 6         24 MY Glut              222
# 7         19 MB Glut              155
# 8          26 P GABA              134
# 9        33 Vascular               67
# 10     01 IT-ET Glut               62
# 11        21 MB Dopa               60
# 12    05 OB-IMN GABA               54
# 13   13 CNU-HYa Glut               22
# 14        12 HY GABA               20
# 15        14 HY Glut               14
# 16        29 CB Glut               13
# 17    04 DG-IMN Glut                9
# 18   11 CNU-HYa GABA                9
# 19 02 NP-CT-L6b Glut                8
# 20   06 CTX-CGE GABA                8
# 21     22 MB-HB Sero                7
# 22        27 MY GABA                5
# 23   07 CTX-MGE GABA                4
# 24     17 MH-LH Glut                3
# 25        28 CB GABA                2
# 26            32 OEC                2
# 27     03 OB-CR Glut                1
# 28   08 CNU-MGE GABA                1
```

## Differential expression analysis

I examined differential expression between pairs of these celltype clusters using Wilcoxon Rank Sum test with `FindMarkers()`. 
```
DefaultAssay(data) <- "RNA"
Idents(data) <- "class_name"

# oligodendrocytes vs glutamatergic neurons
oligo_pglut.de.markers <- FindMarkers(
  data,
  ident.1 = "31 OPC-Oligo",
  ident.2 = "23 P Glut",
  test.use = 'wilcox',
  min.pct = 0.1
)
head(oligo_pglut.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# EPHA7            0  -4.990820 0.043 0.962         0
# ADCY8            0  -5.121001 0.052 0.971         0
# SRRM3            0  -4.881668 0.032 0.945         0
# LOC103787440     0  -5.257025 0.028 0.939         0
# TENM4            0  -4.886707 0.085 0.988         0
# LOC103796371     0  -4.610896 0.056 0.959         0

range(oligo_pglut.de.markers$avg_log2FC)
# [1] -8.254263  5.958942
write.csv(oligo_pglut.de.markers, file = "oligo_pglut_de_markers.csv")
```
```
# glutamatergic neurons vs astrocytes
pglut_astro.de.markers <- FindMarkers(
  data,
  ident.1 = "23 P Glut",
  ident.2 = "30 Astro-Epen",
  test.use = 'wilcox',
  min.pct = 0.1
)
head(pglut_astro.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# LOC103787440     0   5.352007 0.939 0.045         0
# LOC100413948     0   5.000467 0.993 0.100         0
# KSR2             0   4.767807 0.971 0.095         0
# SCN3A            0   4.328221 0.957 0.088         0
# ANO4             0   4.711078 0.939 0.075         0
# PLPPR1           0   4.557623 0.921 0.059         0

range(pglut_astro.de.markers$avg_log2FC)
# [1] -8.625408  6.539726
write.csv(pglut_astro.de.markers, file = "pglut_astro_de_markers.csv")
```
```
# astrocytes vs oligodendrocytes
astro_oligo.de.markers <- FindMarkers(
  data,
  ident.1 = "30 Astro-Epen",
  ident.2 = "31 OPC-Oligo",
  test.use = 'wilcox',
  min.pct = 0.1
)
head(astro_oligo.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# BMPR1B       0   6.470252 0.924 0.036         0
# RFX4         0   6.113532 0.896 0.033         0
# GLI3         0   7.028026 0.866 0.013         0
# SLC4A4       0   6.065063 0.962 0.113         0
# ARHGAP26     0   4.884257 0.935 0.096         0
# GASK1A       0   6.239898 0.865 0.030         0
range(astro_oligo.de.markers$avg_log2FC)
# [1] -5.352012  9.268554
write.csv(astro_oligo.de.markers, file = "astro_oligo_de_markers.csv")
```

## Differential accessibility analysis

Differential chromatin accessibility analysis was conducted using the same approach on the ATAC assay: 

```
DefaultAssay(data) <- "ATAC"
data <- SortIdents(data)

# find peaks associated with gene activities
library(presto)

oligo_pglut.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "31 OPC-Oligo",
  ident.2 = "23 P Glut",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(oligo_pglut.da.peaks, file = "oligo_pglut_da_peaks.csv")

astro_oligo.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "30 Astro-Epen",
  ident.2 = "31 OPC-Oligo",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(astro_oligo.da.peaks, file = "astro_oligo_da_peaks.csv")

pglut_astro.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "23 P Glut",
  ident.2 = "30 Astro-Epen",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(pglut_astro.da.peaks, file = "pglut_astro_da_peaks.csv")
```

## Linking peaks to genes

Correlations between gene expression and accessible chromatin regions can be characterized using `LinkPeaks()` within a specified distance from the transciption start sites identified in pre-processing. Here, peaks associated with a set of marker genes associated with the major brainstem celltypes (oligodendrocytes, OPCs, glutamatergic neurons, astrocytes) found with queries from the [CellxGene]("https://cellxgene.cziscience.com/") Gene Expression tool. This can also be run on the entire gene expression set if time permits (runtime on complete dataset with current resources: ~9 hrs).

```
all_markers <- c("ST18", "CTNNA3", "RNF220", "PIP4K2A", "MBP", "TMEM144", "PDE4B", 
             "PHLPP1", "EDIL3", "PRUNE2", "ELMO1", "TTLL7", "C10orf90", "ENPP2", 
             "DOCK5", "MAP7", "MAN2A1", "QKI", "TMTC2", "SHTN1", "PLCL1", "PLEKHH1", 
             "TF", "PXK", "SLC24A2", "PTPRZ1", "TNR", "LHFPL3", "PCDH15", "EPN2", 
             "SOX6", "SCD5", "CA10", "LUZP2", "VCAN", "LRRC4C", "NXPH1", "NOVA1", "XYLT1", 
             "MMP16", "KAT2B", "MEGF11", "SEMA5A", "SMOC1", "SLC35F1", "GALNT13", "AGAP1", 
             "DGKG", "CDH20", "ARPP21", "DLGAP2", "RALYL", "R3HDM1", "KHDRBS2", "NELL2", "SH3GL2", 
             "SATB2", "KHDRBS3", "SV2B", "HECW1", "GRIN2A", "KCNQ5", "RYR2", "CHRM3", "KCNIP4", "LDB2", 
             "PTPRD", "PRKCB", "IQCJ-SCHIP1", "CHN1", "PHACTR1", "CNKSR2", "TAFA1", "FSTL4",
             "SLC1A2", "NPAS3", "DTNA", "GPC5", "PITPNC1", "GPM6A", "BMPR1B", "RFX4", "MSI2",
             "ADGRV1", "GLIS3", "RYR3", "SLC4A4", "TRPS1", "ATP13A4", "ATP1A2", "ADCY2", "NTM", "ZNRF3", 
             "AQP4", "OBI1-AS1", "LRIG1", "PREX2")

data <- LinkPeaks(
  object = data,
  peak.assay = "ATAC",
  expression.assay = "GENE_ACTIVITIES", # or "RNA"
  peak.slot = "counts",
  method = "pearson",
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = all_markers,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
)
```
Of the 97 queried marker genes, 67 were found to be proximally associated with peaks in the ATAC assay: 
```
unique(Links(data[['ATAC']])$gene)
# [1] "PRUNE2"  "GLIS3"   "SH3GL2"  "SLC24A2" "ZNRF3"   "FSTL4"   "MAN2A1"  "EDIL3"  
# [9] "VCAN"    "ADCY2"   "GPM6A"   "TMEM144" "BMPR1B"  "SCD5"    "SLC4A4"  "LDB2"   
# [17] "SLC35F1" "MAP7"    "QKI"     "EPN2"    "PITPNC1" "MSI2"    "SV2B"    "CHN1"   
# [25] "R3HDM1"  "PLCL1"   "SATB2"   "AGAP1"   "PIP4K2A" "RNF220"  "TTLL7"   "HECW1"  
# [33] "NXPH1"   "PTPRZ1"  "NELL2"   "TMTC2"   "RFX4"    "MEGF11"  "NOVA1"   "PLEKHH1"
# [41] "SMOC1"   "SLC1A2"  "GRIN2A"  "XYLT1"   "PRKCB"   "SHTN1"   "DLGAP2"  "DOCK5"  
# [49] "AQP4"    "DTNA"    "CDH20"   "PHLPP1"  "MBP"     "DGKG"    "ATP13A4" "PXK"    
# [57] "LRIG1"   "ST18"    "PREX2"   "MMP16"   "KHDRBS3" "ENPP2"   "TRPS1"   "TF"     
# [65] "KAT2B"   "ARPP21"  "ATP1A2"  "TNR"     "CNKSR2" 
```

Only two of these genes (SLC4A4, RFX4) were also identified as differentially expressed between oligodendrocyte+OPCs and astrocytes+ependymal cell clusters. The incongruence between highly expressed genes revealed in the RNA assay and those found to be associated with peaks in the ATAC assay is to an extent expected since cell-specific gene expression might be mediated by cellular context and enhancers well outside the specified search proximity.

`CoveragePlot()` visualizes differential accessibility between celltypes and linked peaks associated with a gene or region of interest.

```
idents.plot <- c("31 OPC-Oligo", "23 P Glut", "30 Astro-Epen")

CoveragePlot(
  object = data,
  assay = "ATAC",
  expression.assay = "SCT",
  region = "GRIN2A", # schizophrenia associated gene
  extend.upstream = 1000,
  extend.downstream = 1000,
  idents = idents.plot,
  annotation = TRUE
)
```

# Addendum: Celltype prediction with primate data and SingleR

Cell annotation can be drastically impacted by number of factors, including but not limited to the quality of the reference itself. Using MapMyCells, I predicted celltypes using a bootstrap-like algorithm to compare expression profiles in primate tissue with a well annotated wild type/healthy mouse atlas. However, evolutionary distance between reference and query should also be considered. The emergence of the mouse lemur _Microcebus murinus_ as a model primate system enables a lightweight test of the effects of using a reference dataset of a more closely related system on cell annotation and downstream DE and DA analyses. 

I downloaded RNA-seq data derived from brainstem and cortex tissue from the [Tabula murnius]("https://tabula-microcebus.sf.czbiohub.org/whereisthedata") and converted the .h5ad to a Seurat object.
```
library(anndataR)
library(Seurat)
library(SingleR)

# reference data preparation
reference <- read_h5ad("microcebus_murinus_brainstem_and_cortex.h5ad")
reference <- reference$as_Seurat()
```

Similar to the annotation construction and chromosome renaming prior to generating the multiome Seurat object, here the reference data's rownames needed to be converted from ENSEMBL gene identifiers to gene names to match the query data. This was accomplished by querying the ENSEMBL database with ```biomaRt```.
```
listMarts() 
#               biomart                version
#1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 115
#2   ENSEMBL_MART_MOUSE      Mouse strains 115
#3     ENSEMBL_MART_SNP  Ensembl Variation 115
#4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 115

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
mouselemur_mart <- useDataset(dataset = 'mmurinus_gene_ensembl', mart = ensembl, verbose = TRUE)
attributes <- c("ensembl_gene_id", "external_gene_name")
gene_id_data <- getBM(attributes = attributes, mart = mouselemur_mart)
head(gene_id_data)

id_to_symbol <- setNames(gene_id_data$external_gene_name, gene_id_data$ensembl_gene_id)

# get ref's counts matrix and assign new rownames
rownames(reference) <- id_to_symbol[rownames(reference)]
reference <- subset(reference, features = rownames(reference)[rownames(reference) != ""])
```
A Wilcoxon Rank Sum Test was conducted with `SingleR()`, first using log-normalized data. The result is a large dataframe, which was saved as an R object.
```
data <- LoadSeuratRds("multiome_data_final.Rds")
norm_counts <- LayerData(data, assay = "RNA", layer = 'data') # LogNormalized
raw_counts <- LayerData(data, assay = "RNA", layer = 'counts')

# run SingleR
cellannotation_microcebus <- SingleR(test = norm_counts,
                  ref = reference[["RNA"]]$X, 
                  labels = reference$cell_type,
                  de.method = 'wilcox')
saveRDS(cellannotation_microcebus,file = "singleR_microcebus_cellannotation.Rds")
```
Cell x gene expression patterns were visualized with functions from the `pheatmap` library.
```
plotScoreHeatmap(cellannotation_microcebus)
plotDeltaDistribution(cellannotation_microcebus, ncol = 4, dots.on.top = FALSE)
```
![Venn diagram showing cell annotations and quality scores from SingleR run](figs/04_SingleR_heatmap_microebus_annotations.png)

These annotations were added to the original (and now very large) Seurat object's metadata to faciliate comparison of celltype annotations.
```
rownames(cellannotation_microcebus)[1:5] # make sure you have cell IDs
data <- Seurat::AddMetaData(data, cellannotation_microcebus$pruned.labels, col.name = 'SingleR_annt_microcebus_ref')
saveRDS(data,file = "multiome_data_final.Rds")
```
Expectedly, the majority of the cells are identified as oligodendroctyes, but 38% of the cells are assigned to either neurons (not subdivided as in MapMyCells) or astrocytes. The MapMyCells WMB annotations and the SingleR _Microcebus_ annotations do not have perfectly congruent celltypes, but between both strategies there is similarity in number of assignments to either oligodendrocytes or OPCs (MMC-WMB = 9098, SingleR-Mm = 10258), as well as the labeling of either astrocytes or ependymal cells (MMC-WMB = 1407, SingleR-Mm = 1588). Furthermore, 303 cells were also pruned during the SingleR prediction due to low quality annotation.

```
microcebus_annotations <- readRDS("singleR_microcebus_cellannotation.Rds")
data <- readRDS("multiome_data_final.Rds")
data <- SetIdent(data, value = 'SingleR_annt_microcebus_ref')
DimPlot(data, reduction = "umap", repel = T, label = T) #+ NoLegend()
# Error in data.frame(..., check.names = FALSE) : 
#   arguments imply differing number of rows: 16227, 15924
# In addition: Warning message:
#   Removing 303 cells missing data for vars requested 
# referenced here: https://github.com/satijalab/seurat/issues/10180

data.frame(Cell_counts=head(sort(table(microcebus_annotatoins$labels),decreasing=T),34))
# Cell_counts.Var1 Cell_counts.Freq
# 1                 oligodendrocyte             8202
# 2                          neuron             4533
# 3                       astrocyte             1578
# 4  oligodendrocyte precursor cell              896
# 5                GABAergic neuron              351
# 6                      macrophage              249
# 7            glutamatergic neuron              183
# 8                    stromal cell               91
# 9      capillary endothelial cell               30
# 10       myelinating Schwann cell               21
# 11   non-myelinating Schwann cell               18
# 12                         B cell               14
# 13                     lymphocyte               13
# 14                       pericyte               13
# 15                 ependymal cell               10
# 16 choroid plexus epithelial cell                9
# 17            leptomeningeal cell                8
# 18                        unknown                8
```
Visualization of shared features between the _Microcebus murinus_ data and _Callithrix jacchus_ genome reveals an 840 gene discrepancy between the references:
```
library(VennDiagram)
venn.diagram(
  disable.logging = TRUE,
  margin =.1,
  x = list(microcebus_genes, callithtrix_genes, rna_genes, gene_activities),
  category.names = c("Microcebus" , "Callithrix", "RNA Assay", "Gene Activities Assay"),
  filename = 'microcebus_vs_callithrix_venn_diagramm_multiassay.png',
  output=TRUE
)
```
While irrelevant for this experiment due to their absence in the query data, this result highlights the nontriviality of reference choice for both cluster identification and downstream analyses. Addiionally, the near 9000 genes shared by the RNA expression assay and the _Callithrix_ reference genome but not inferred from `GeneActivity()` identificaiton further highlights the signal discordance between direct and indirect gene expression quantification.

## References

### Data

Dataset: Fenna M. Krienen, Kirsten M. Levandowski, Heather Zaniewski, Ricardo C.H. del Rosario, Margaret E. Schroeder, Melissa Goldman, Alyssa Lutservitz, Qiangge Zhang, Katelyn X. Li, Victoria F. Beja-Glasser, Jitendra Sharma, Tay Won Shin, Abigail Mauermann, Alec Wysoker, James Nemesh, Seva Kashin, Josselyn Vergara, Gabriele Chelini, Jordane Dimidschstein, Sabina Berretta, Ed Boyden, Steven A. McCarroll, Guoping Feng (2023) A Marmoset Brain Cell Census Reveals Persistent Influence of Developmental Origin on Neurons. Available from https://assets.nemoarchive.org/dat-1je0mn3

UCSC Genome Browser assembly ID: "calJac4". Perez et al. The UCSC Genome Browser database: 2025 update. Nucleic Acids Research 2025 PMID: 39460617, DOI: 10.1093/nar/gkae974 

### Relevant Literature

Ezran, Camille, et al. "Mouse lemur cell atlas informs primate genes, physiology and disease." Nature 644.8075 (2025): 185-196. 

Gontarz, Paul, et al. "Comparison of differential accessibility analysis strategies for ATAC-seq data." Scientific reports 10.1 (2020): 10150.

Krienen, Fenna M., et al. "A marmoset brain cell census reveals regional specialization of cellular identities." Science advances 9.41 (2023): eadk3986. 

Teo, Alan Yue Yang, et al. "Best practices for differential accessibility analysis in single-cell epigenomics." Nature communications 15.1 (2024): 8805. 

Yan, Feng, et al. "From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis." Genome biology 21.1 (2020): 22. 

Yao, Zizhen, et al. "A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain." Nature 624.7991 (2023): 317-332. 




To give your project a background in the portfolio page, just add the img tag to the front matter like so:

    ---
    layout: page
    title: project
    description: a project with a background image
    img: /assets/img/12.jpg
    ---

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/3.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Caption photos easily. On the left, a road goes through a tunnel. Middle, leaves artistically fall in a hipster photoshoot. Right, in another hipster photoshoot, a lumberjack grasps a handful of pine needles.
</div>
<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    This image can also have a caption. It's like magic.
</div>

You can also put regular text between your rows of images, even citations {% cite einstein1950meaning %}.
Say you wanted to write a bit about your project before you posted the rest of the images.
You describe how you toiled, sweated, _bled_ for your project, and then... you reveal its glory in the next row of images.

<div class="row justify-content-sm-center">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm-4 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/11.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    You can also have artistically styled 2/3 + 1/3 images, like these.
</div>

The code is simple.
Just wrap your images with `<div class="col-sm">` and place them inside `<div class="row">` (read more about the <a href="https://getbootstrap.com/docs/4.4/layout/grid/">Bootstrap Grid</a> system).
To make images responsive, add `img-fluid` class to each; for rounded corners and shadows use `rounded` and `z-depth-1` classes.
Here's the code for the last row of images above:

{% raw %}

```html
<div class="row justify-content-sm-center">
  <div class="col-sm-8 mt-3 mt-md-0">
    {% include figure.liquid path="assets/img/6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
  </div>
  <div class="col-sm-4 mt-3 mt-md-0">
    {% include figure.liquid path="assets/img/11.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
  </div>
</div>
```

{% endraw %}

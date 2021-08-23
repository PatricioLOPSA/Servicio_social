---
title: "Analizando distintos rangos taxonómicos"
author: "Patricio López Sánchez"
date: "8/6/2021"
output:
  html_document:
    keep_md: yes
---




```r
library(HMP2Data)
library(phyloseq)
library(dplyr)
library(stringr)
library(DESeq2)
library(edgeR)
library(eulerr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(PCAtools)
library(infotheo)
```


Repetimos los mismos pasos de _wrangling_.


```r
T2D = T2D16S()
T2D_otu = otu_table(T2D) %>% as.data.frame()
T2D_samples = sample_data(T2D) %>% as.data.frame()
T2D_tax = tax_table(T2D)
T2D_gut_IDs = subset(rownames(T2D_samples), T2D_samples$sample_body_site == "feces")
T2D_otu_gut = select(T2D_otu, T2D_gut_IDs)
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(T2D_gut_IDs)` instead of `T2D_gut_IDs` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

```r
IDs_no_tail = str_sub(T2D_gut_IDs, 1, -7)
IDs_no_tailhead = str_sub(IDs_no_tail,29,-1)

head(as.data.frame(T2D_gut_IDs))
```

```
##                                      T2D_gut_IDs
## 1   HMP2_J00825_1_ST_T0_B0_0120_ZN9YTFN-01_AA31J
## 2 HMP2_J00826_1_ST_T0_B0_0120_ZN9YTFN-1011_AA31J
## 3 HMP2_J00827_1_ST_T0_B0_0120_ZN9YTFN-1012_AA31J
## 4 HMP2_J00828_1_ST_T0_B0_0120_ZN9YTFN-1013_AA31J
## 5 HMP2_J00829_1_ST_T0_B0_0120_ZN9YTFN-1014_AA31J
## 6 HMP2_J00830_1_ST_T0_B0_0120_ZN9YTFN-1015_AA31J
```

```r
head(as.data.frame(IDs_no_tail))
```

```
##                                IDs_no_tail
## 1   HMP2_J00825_1_ST_T0_B0_0120_ZN9YTFN-01
## 2 HMP2_J00826_1_ST_T0_B0_0120_ZN9YTFN-1011
## 3 HMP2_J00827_1_ST_T0_B0_0120_ZN9YTFN-1012
## 4 HMP2_J00828_1_ST_T0_B0_0120_ZN9YTFN-1013
## 5 HMP2_J00829_1_ST_T0_B0_0120_ZN9YTFN-1014
## 6 HMP2_J00830_1_ST_T0_B0_0120_ZN9YTFN-1015
```

```r
head(as.data.frame(IDs_no_tailhead))
```

```
##   IDs_no_tailhead
## 1      ZN9YTFN-01
## 2    ZN9YTFN-1011
## 3    ZN9YTFN-1012
## 4    ZN9YTFN-1013
## 5    ZN9YTFN-1014
## 6    ZN9YTFN-1015
```

```r
OTUs_clinical = T2D_otu_gut
colnames(OTUs_clinical) <- NULL
colnames(OTUs_clinical) <- IDs_no_tailhead

df_clinico = read.csv("clinical_data", sep = "\t")

salud = filter(df_clinico, df_clinico$CL4 == "Healthy")
infeccion = filter(df_clinico, 
                   df_clinico$CL4 == "Infection" | df_clinico$CL4 == "Infection_L")


OTUs_clinical_unique = OTUs_clinical[!duplicated(names(OTUs_clinical))]
OTUs_clinical_IDs = colnames(OTUs_clinical_unique)


salud_IDs = salud$VisitID
salud_IDs_flt = subset(salud_IDs, salud_IDs %in% OTUs_clinical_IDs)


infeccion_IDs = infeccion$VisitID
infeccion_IDs_flt = subset(infeccion_IDs, infeccion_IDs %in% OTUs_clinical_IDs)

#Obtenemos la matriz de OTUs que contiene ambas condiciones anotadas

OTUs_clinical_salud_flt = select(OTUs_clinical_unique, salud_IDs_flt)
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(salud_IDs_flt)` instead of `salud_IDs_flt` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

```r
OTUs_clinical_infeccion_flt = select(OTUs_clinical_unique, infeccion_IDs_flt)
```

```
## Note: Using an external vector in selections is ambiguous.
## ℹ Use `all_of(infeccion_IDs_flt)` instead of `infeccion_IDs_flt` to silence this message.
## ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This message is displayed once per session.
```

```r
OTUs_clinical_complete = cbind(OTUs_clinical_salud_flt, OTUs_clinical_infeccion_flt)
OTUs_clinical_complete_nonzero = OTUs_clinical_complete + 1


coldata_salud = cbind(salud_IDs_flt, rep("Saludable", length(salud_IDs_flt))) %>% as.data.frame()
coldata_infeccion = cbind(infeccion_IDs_flt, rep("Infectados", length(infeccion_IDs_flt))) %>% as.data.frame()

colnames(coldata_salud) = c("Visit_ID","condition")
colnames(coldata_infeccion) = c("Visit_ID","condition")

#unimos los dos dfs.

coldata = rbind(coldata_salud, coldata_infeccion)
coldata$condition = factor(coldata$condition)
```


## Creando el objeto de Phyloseq

El siguiente objetivo es crear un objeto de phyloseq que contenga la información de la tabla de taxonomías que le corresponden a cada OTU, al igual que la nueva tabla de metadatos que le corresponen a cada muestra. Para crear el objeto de phyloseq, necesitamos lo siguiente.

- Tabla de OTUs como matriz.
- Tabla de metadatos de las muestras como data frame
- Tabla de taxonomia de los OTUs como matriz



```r
#OTUs como matriz
OTUs_clinical_complete = as.matrix(OTUs_clinical_complete)

#sample data como DF
coldata_mod = coldata[,-1] %>% as.data.frame()
rownames(coldata_mod) = coldata[,1]
colnames(coldata_mod) = "condition"

#Tabla de taxonomia como matriz ya la tenemos en el objeto original T2D

#Creamos objeto phyloseq desginando las matrices como componenten del objeto phyloseq
OTU = otu_table(OTUs_clinical_complete, taxa_are_rows = TRUE)
samp = sample_data(coldata_mod)



T2D_phseq = phyloseq(OTU, T2D_tax, samp)
T2D_phseq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 12062 taxa and 445 samples ]
## sample_data() Sample Data:       [ 445 samples by 1 sample variables ]
## tax_table()   Taxonomy Table:    [ 12062 taxa by 7 taxonomic ranks ]
```

Filtramos OTUs con una gran prevalencia de 0 a lo largo de todas las muestras.



```r
# Filtramos OTUs sin conteos
T2D_phseq = filter_taxa(T2D_phseq, function(x) sum(x) > (0), TRUE)

T2D_phseq = filter_taxa(T2D_phseq, function(x) sum( x != 0 ) > (0.1*length(x)), TRUE)

T2D_phseq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1477 taxa and 445 samples ]
## sample_data() Sample Data:       [ 445 samples by 1 sample variables ]
## tax_table()   Taxonomy Table:    [ 1477 taxa by 7 taxonomic ranks ]
```








```r
#intentamos combinar OTUs de un mismo nivel taxonómico

T2D_genus = tax_glom(T2D_phseq, taxrank = "Genus")
T2D_family = tax_glom(T2D_phseq, taxrank = "Family")
```


```r
otus_genus <- otu_table(T2D_phseq) %>% as.data.frame()
otus_genus <- otus_genus +1
genus_dds <- DESeqDataSetFromMatrix(otus_genus, colData = coldata_mod, design = ~ condition)
```

```
## converting counts to integer mode
```

```r
genus_vst <- varianceStabilizingTransformation(genus_dds, blind = TRUE)
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```

```r
genus_vst_matriz <- assay(genus_vst)

pheatmap(genus_vst_matriz, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = F)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-6-1.png)<!-- -->



## Expresión diferencial con edgeR

Definimos la funcion que convierte los objetos de Phyloseq a DGE de edgeR


```r
phyloseq_to_edgeR = function(physeq, group, method="TMM", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}
```



```r
dge = phyloseq_to_edgeR(T2D_phseq, group="condition")
# Perform binary test
et = exactTest(dge)
summary(decideTestsDGE(et, p.value = 0.01))
```

```
##        Saludable-Infectados
## Down                    108
## NotSig                 1178
## Up                      191
```

```r
# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.001
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(T2D_phseq)[rownames(sigtab), ], "matrix"))
dim(sigtab)
```

```
## [1] 205  18
```

```r
head(sigtab)
```

```
##         Kingdom         Phylum          Class            Order
## 179585 Bacteria     Firmicutes     Clostridia    Clostridiales
## 188262 Bacteria     Firmicutes     Clostridia    Clostridiales
## 293883 Bacteria     Firmicutes     Clostridia    Clostridiales
## 287691 Bacteria     Firmicutes     Clostridia    Clostridiales
## 186090 Bacteria     Firmicutes     Clostridia    Clostridiales
## 329688 Bacteria Actinobacteria Coriobacteriia Coriobacteriales
##                   Family                 Genus   Species     logFC   logCPM
## 179585   Ruminococcaceae                  <NA>      <NA> -2.749790 9.924167
## 188262              <NA>                  <NA>      <NA> -2.135351 8.342829
## 293883   Veillonellaceae Phascolarctobacterium      <NA>  3.348296 9.086197
## 287691    Clostridiaceae                  <NA>      <NA> -2.099671 9.565036
## 186090   Ruminococcaceae          Ruminococcus      <NA>  3.358083 9.625383
## 329688 Coriobacteriaceae           Collinsella stercoris  3.657399 9.964333
##              PValue          FDR  Kingdom         Phylum          Class
## 179585 9.750090e-29 1.440088e-25 Bacteria     Firmicutes     Clostridia
## 188262 2.042298e-26 1.508237e-23 Bacteria     Firmicutes     Clostridia
## 293883 3.331432e-24 1.640175e-21 Bacteria     Firmicutes     Clostridia
## 287691 2.373716e-21 8.751586e-19 Bacteria     Firmicutes     Clostridia
## 186090 2.962622e-21 8.751586e-19 Bacteria     Firmicutes     Clostridia
## 329688 4.830961e-21 9.518507e-19 Bacteria Actinobacteria Coriobacteriia
##                   Order            Family                 Genus   Species
## 179585    Clostridiales   Ruminococcaceae                  <NA>      <NA>
## 188262    Clostridiales              <NA>                  <NA>      <NA>
## 293883    Clostridiales   Veillonellaceae Phascolarctobacterium      <NA>
## 287691    Clostridiales    Clostridiaceae                  <NA>      <NA>
## 186090    Clostridiales   Ruminococcaceae          Ruminococcus      <NA>
## 329688 Coriobacteriales Coriobacteriaceae           Collinsella stercoris
```


```r
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


Volcano plot



```r
#Primero categorizamos si los OTUS se encuentran UP o DOWN
volcano_plot = function(edgeR_DEgenes, logfc, pvalue){

  edgeR_DEgenestax = edgeR_DEgenes$genes
  edgeR_DEgenesTable=edgeR_DEgenes$table
edgeR_DEgenesTable$p.adj = p.adjust(edgeR_DEgenesTable$PValue, "fdr")
edgeR_DEgenesTable$DA = "NS"
edgeR_DEgenesTable$DA[edgeR_DEgenesTable$logFC >(logfc) & edgeR_DEgenesTable$p.adj <(pvalue)] = "UP"
edgeR_DEgenesTable$DA[edgeR_DEgenesTable$logFC < -(logfc) & edgeR_DEgenesTable$p.adj <(pvalue)] = "DOWN"
edgeR_DEgenesTable$DA = as.factor(edgeR_DEgenesTable$DA)
edgeR_DEgenesTable$delabel <- NA
edgeR_DEgenesTable$delabel[edgeR_DEgenesTable$DA != "NS"] <- edgeR_DEgenestax$Phylum[edgeR_DEgenesTable$DA != "NS"]


p <- ggplot(data=edgeR_DEgenesTable, aes(x=logFC, y=-log10(p.adj), col=DA, label=delabel)) + geom_point(size= 1) + theme_minimal() + geom_text_repel()

p2 = p + geom_vline(xintercept=c(-logfc, logfc), col="red") +
        geom_hline(yintercept=-log10(pvalue), col="red")

mycolors <- c("lightgreen", "coral", "darkgray")
names(mycolors) <- c("DOWN", "UP", "NS")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3 = p2 + scale_color_manual(values=mycolors) + ylim(0,40)
return(p3)
}
```



```r
volcano_plot(et, 1, 0.01)
```

```
## Warning: Removed 1329 rows containing missing values (geom_text_repel).
```

```
## Warning: ggrepel: 130 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-11-1.png)<!-- -->






```r
DAotus <- et$table
DAotus$p.adj <- p.adjust(DAotus$PValue,method = "fdr") 

DAotus_flt <- filter(DAotus, DAotus$p.adj < 0.01)

DAotus_table <- filter(OTUs_clinical_complete_nonzero,
                       rownames(OTUs_clinical_complete_nonzero) %in% rownames(DAotus_flt))

DAtax <- filter(as.data.frame(T2D_tax), rownames(T2D_tax) %in% rownames(DAotus_table))
```


Ahora convertimos a un objeto de phyloseq que solo contenga los OTUs y la tabla de taxonomías de los OTUs DA. Haremos el merge a nivel de genus y familia y orden.



```r
OTUS <- otu_table(as.matrix(DAotus_table), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(DAtax))
SAMP <- sample_data(coldata_mod)

DA<- phyloseq(OTUS, TAX, SAMP )

DA_genus <- tax_glom(DA, taxrank = "Genus")
DA_family <- tax_glom(DA, taxrank = "Family")
DA_order <- tax_glom(DA, taxrank = "Order")

DA_genus
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 29 taxa and 445 samples ]
## sample_data() Sample Data:       [ 445 samples by 1 sample variables ]
## tax_table()   Taxonomy Table:    [ 29 taxa by 7 taxonomic ranks ]
```

```r
DA_family
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 22 taxa and 445 samples ]
## sample_data() Sample Data:       [ 445 samples by 1 sample variables ]
## tax_table()   Taxonomy Table:    [ 22 taxa by 7 taxonomic ranks ]
```

```r
DA_order
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 15 taxa and 445 samples ]
## sample_data() Sample Data:       [ 445 samples by 1 sample variables ]
## tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]
```



```r
species_dds <- DESeqDataSetFromMatrix(as.data.frame(otu_table(DA)), colData = coldata_mod, design = ~ condition)
```

```
## converting counts to integer mode
```

```r
genus_dds <- DESeqDataSetFromMatrix(as.data.frame(otu_table(DA_genus)), colData = coldata_mod, design = ~ condition)
```

```
## converting counts to integer mode
```

```r
family_dds <- DESeqDataSetFromMatrix(as.data.frame(otu_table(DA_family)), colData = coldata_mod, design = ~ condition)
```

```
## converting counts to integer mode
```

```r
order_dds <-  DESeqDataSetFromMatrix(as.data.frame(otu_table(DA_order)), colData = coldata_mod, design = ~ condition)
```

```
## converting counts to integer mode
```

```r
species_vst <- varianceStabilizingTransformation(species_dds, blind = TRUE) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```

```r
genus_vst <- varianceStabilizingTransformation(genus_dds, blind = TRUE) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```

```r
family_vst <- varianceStabilizingTransformation(family_dds, blind = TRUE) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```

```r
order_vst <- varianceStabilizingTransformation(order_dds, blind = TRUE) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```


## Heatmaps de OTUs DA a diferentes niveles taxonómicos


```r
pheatmap(species_vst, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = T)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```r
pheatmap(genus_vst, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = T)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

```r
pheatmap(family_vst, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = T)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-15-3.png)<!-- -->

```r
pheatmap(order_vst, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = T)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-15-4.png)<!-- -->

```r
cor_genus = cor(genus_vst, method = "pearson")

pheatmap(cor_genus, fontsize = 7, border_color = NA, annotation_col = coldata_mod, show_rownames = F, show_colnames = F, labels_row = NULL, annotation_legend = T, cluster_cols = T)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-15-5.png)<!-- -->


Probamos con PCA


```r
pca_vst = pca(family_vst, metadata = as.matrix(coldata_mod))
screeplot(pca_vst, getComponents(pca_vst, 1:5))
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
pairsplot(pca_vst, axisLabSize = 3 , colby = 'condition', components = getComponents(pca_vst, 1:5))
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

```r
biplot(pca_vst, x="PC3", y = "PC1" ,colby = 'condition', showLoadings = F, lab = NULL, pointSize = 1.5)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-16-3.png)<!-- -->

```r
biplot(pca_vst, x="PC2", y = "PC1" ,colby = 'condition', showLoadings = F, lab = NULL, pointSize = 1.5)
```

![](Abundancia_diferencial_3_files/figure-html/unnamed-chunk-16-4.png)<!-- -->







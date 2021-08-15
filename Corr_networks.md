---
title: "Correlación entre OTUs"
author: "Patricio López Sánchez"
date: "8/12/2021"
output:
  html_document:
    keep_md: yes
---




# Correlación entre OTUs


```r
library(HMP2Data)
library(phyloseq)
library(dplyr)
library(stringr)
library(DESeq2)
library(edgeR)
library(eulerr)
library(ggplot2)
library(pheatmap)
library(igraph)
library(Hmisc)
```


Primero necesitamos obtener las matricesde abundandia que corresponen a infectados y a saludables por separado.


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


Realizamos un paso de filtrado para remover OTUs sin conteos o que fueron muy raros a lo largo de todas las muestras de interés. Esto nos ayuda a remover un poco del exceso de ceros en la matriz, y evitar que al calcular la correlación entre OTUs obtengamos una $S = 0$. Para usaremos las herramientas de phyloseq, por lo que necesitamos  obetener al objeto de phyloseq para cada condición.


```r
#OTUs como matriz
OTUs_sld = as.matrix(OTUs_clinical_salud_flt)
OTUs_inf = as.matrix(OTUs_clinical_infeccion_flt)


#Tabla de taxonomia como matriz ya la tenemos en el objeto original T2D

#Creamos objeto phyloseq desginando las matrices como componenten del objeto phyloseq
OTUsalud = otu_table(OTUs_sld, taxa_are_rows = TRUE)
OTUinf = otu_table(OTUs_inf, taxa_are_rows = T)




sld_phseq = phyloseq(OTUsalud, T2D_tax)
inf_phseq = phyloseq(OTUinf, T2D_tax)
```



```r
sld_phseq_flt = filter_taxa(sld_phseq, function(x) sum(x) > (0), TRUE)
sld_phseq_flt
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 5691 taxa and 343 samples ]
## tax_table()   Taxonomy Table:    [ 5691 taxa by 7 taxonomic ranks ]
```

```r
inf_phseq_flt = filter_taxa(inf_phseq, function(x) sum(x) > (0), TRUE )
inf_phseq_flt
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4296 taxa and 102 samples ]
## tax_table()   Taxonomy Table:    [ 4296 taxa by 7 taxonomic ranks ]
```





Normalizamos de manera adecuada para calcular correlación


```r
inf_nonzero = (otu_table(inf_phseq_flt) +1) %>% as.data.frame()
sld_nonzero = (otu_table(sld_phseq_flt) +1) %>%  as.data.frame()

inf_dds = DESeqDataSetFromMatrix(inf_nonzero, colData = coldata_infeccion, design = ~ 1)
```

```
## converting counts to integer mode
```

```r
sld_dds = DESeqDataSetFromMatrix(sld_nonzero, colData = coldata_salud, design = ~1)
```

```
## converting counts to integer mode
```

```r
inf_vst = varianceStabilizingTransformation(inf_dds, blind = T) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```

```r
sld_vst = varianceStabilizingTransformation(sld_dds, blind = T) %>% assay()
```

```
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
```


Se necesitan transponer las matrices para calcular la correlación entre OTUs.


```r
t_infectados = t(inf_vst)
t_saludable = t(sld_vst)


p_infectados <- rcorr(t_infectados, type = "pearson")
p_salud <- rcorr(t_infectados, type = "pearson")

s_infectados <- rcorr(t_infectados, type = "spearman")
s_salud <- rcorr(t_infectados, type = "spearman")
```

Aplanamos la matriz de correlación con sus valores P.



```r
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}
```



```r
prs_infect = flattenCorrMatrix(p_infectados$r, p_infectados$P)
prs_salud = flattenCorrMatrix(p_salud$r, p_salud$P)

sp_infect = flattenCorrMatrix(s_infectados$r, s_infectados$P)
sp_salud = flattenCorrMatrix(s_salud$r, s_salud$P)


head(prs_infect)
```


Corregimos los valores de P para pruebas múltiples



```r
prs_infect$p.adj = p.adjust(prs_infect$p, method = "BH")
prs_salud$p.adj = p.adjust(prs_salud$p, method = "BH")

sp_infect$p.adj = p.adjust(sp_infect$p, method = "BH")
sp_salud$p.adj = p.adjust(sp_salud$p, method = "BH")

#write.csv(prs_infect, "prs_infect.csv")
#write.csv(prs_salud, "prs_salud.csv")
#write.csv(sp_infect, "sp_infect.csv")
#write.csv(sp_salud, "sp_salud.csv")
```





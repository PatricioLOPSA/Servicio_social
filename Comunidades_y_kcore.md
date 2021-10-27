---
title: "Algoritmos de comunidades"
author: "Patricio López Sánchez"
date: "10/2/2021"
output: 
  html_document:
    keep_md: TRUE
---





```r
library(vroom)
library(igraph)
library(ggraph)
library(ggplot2)
library(graphlayouts)
library(HMP2Data)
library(phyloseq)
library(dplyr)
library(RColorBrewer)
library(ggforce)
library(pheatmap)
```




```r
inf = vroom("adj_mat_inf.csv") %>% as.matrix()
inf = inf[,-1]
infn = colnames(inf)

rownames(inf) = infn
diag(inf) = 0


sld = vroom("adj_mat_sld.csv") %>% as.matrix()
sld = sld[,-1]
sldn = colnames(sld)
rownames(sld) = sldn
diag(sld) = 0
```

Cargamos matrices de adyacencia como objetos de igraph.


```r
g_inf = graph_from_adjacency_matrix(inf, mode = "undirected", weighted = T)
g_sld = graph_from_adjacency_matrix(sld, mode = "undirected", weighted = T)
```


Cargamos información taxonómica de los OTUs a los mismos objetos de igraph.


```r
T2D = T2D16S()
tax = tax_table(T2D) %>% as.data.frame()

tax_sld = filter(tax, rownames(tax) %in% V(g_sld)$name)
tax_inf = filter(tax, rownames(tax) %in% V(g_inf)$name)

V(g_inf)$Phylum = tax_inf$Phylum
V(g_sld)$Phylum = tax_sld$Phylum

V(g_inf)$Class = tax_inf$Class
V(g_sld)$Class = tax_sld$Class


V(g_inf)$Order= tax_inf$Order
V(g_sld)$Order = tax_sld$Order

V(g_inf)$Family = tax_inf$Family
V(g_sld)$Family = tax_sld$Family


V(g_inf)$Genus = tax_inf$Genus
V(g_sld)$Genus = tax_sld$Genus
```


# Descomposición en k-nucleos

Queremos explorar si existe algún orden o patrón preferencial de los OTUs dentro de los nucleos (cores) de la red.


```r
g_sld_core = g_sld
isolated = which(degree(g_sld_core)==0)
g_sld_core = delete_vertices(g_sld_core, isolated)
E(g_sld_core)$corr = E(g_sld_core)$weight
E(g_sld_core)$peso = abs(E(g_sld_core)$weight)
E(g_sld_core)$weight = abs(E(g_sld_core)$weight)
V(g_sld_core)$size = 1.5
V(g_sld_core)$label = NA
E(g_sld_core)$color = ifelse(E(g_sld_core)$corr >0, "darkseagreen3","red")
E(g_sld_core)$signo = ifelse(E(g_sld_core)$corr >0, "positivo","negativo")
E(g_sld_core)$width = abs(E(g_sld_core)$corr) / 2



coreness_sld = coreness(g_sld_core)
V(g_sld_core)$core = coreness_sld
V(g_sld_core)$color =  V(g_sld_core)$core
V(g_sld_core)$frame.color =  "darkgrey"
V(g_sld_core)$grado = degree(g_sld_core)

g_sld_core$palette = categorical_pal(8)



l = layout_with_centrality(g_sld_core, V(g_sld_core)$core)

# Debido a que igraph solo acepta 8 colores categóricos, el color de los nodos
# indica diferentes cores, pero se repiten (hay 41 cores). Por lo tanto en esta
# red se resalta más el signo de las interacciones
plot(g_sld_core, layout = l)
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
# En esta viz. observamos solo a las interacciones negativas
g_sld_neg = g_sld_core
E(g_sld_neg)$color = ifelse(E(g_sld_core)$corr >0, NA,"red")
l_neg = layout_with_centrality(g_sld_neg, V(g_sld_neg)$core)

plot(g_sld_neg, layout=l_neg)
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
# En esta viz. el color de los nodos va acorde con la familia a la que pertenece
# cada OTU.
ggraph(g_sld_core,layout = "centrality", centrality = coreness(g_sld_core))+
  draw_circle(use = "cent")+
   annotate_circle(coreness(g_sld_core),format="",pos="bottom", col = "red") +
  geom_edge_link0(width=0.02,colour="lightgrey")+
  geom_node_point(aes(colour=Family),size=0.5)+ 
  theme_graph()
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

La descomposición de los k-cores de la red de saludables indica que no existe realmente una arreglo preferencial de los OTUs en alguno de los cores a nivel de familia, sin embargo pordríamos analizar esto más a detalle con el análisis de enriquecimiento.


Realizamos el mismo procedimiento para la red de infectados.

```r
g_inf_core = g_inf
V(g_inf_core)$label = NA

isolated = which(degree(g_inf_core)==0)
g_inf_core= delete_vertices(g_inf_core, isolated)
E(g_inf_core)$corr = E(g_inf_core)$weight
E(g_inf_core)$peso = abs(E(g_inf_core)$weight)
E(g_inf_core)$weight = abs(E(g_inf_core)$weight)
V(g_inf_core)$size = 1.5
E(g_inf_core)$color = ifelse(E(g_inf_core)$corr >0, "darkseagreen3","red")
E(g_inf_core)$signo = ifelse(E(g_inf_core)$corr >0, "positivo","negativo")
E(g_inf_core)$width = abs(E(g_inf_core)$corr) / 2

coreness_inf = coreness(g_inf_core)
V(g_inf_core)$core = coreness_inf
V(g_inf_core)$color =  V(g_inf_core)$core
V(g_inf_core)$frame.color =  "darkgrey"
V(g_inf_core)$grado = degree(g_inf_core) 
V(g_inf_core)$grado_normalizado = degree(g_inf_core) %>% CoDiNA::normalize()
g_inf_core$palette = categorical_pal(8)



li = layout_with_centrality(g_inf_core, V(g_inf_core)$core)


plot(g_inf_core, layout = li )
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
g_inf_neg = g_inf_core
E(g_inf_neg)$color = ifelse(E(g_inf_neg)$corr >0, NA,"red")
li_neg = layout_with_centrality(g_inf_neg, V(g_inf_neg)$core)

plot(g_inf_neg, layout = li_neg )
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
ggraph(g_inf_core,layout = "centrality", centrality = coreness_inf)+
  draw_circle(use = "cent")+
   annotate_circle(coreness_inf,format="",pos="top", col = "red") +
  geom_edge_link0(width=0.02, colour= "lightgrey")+
  geom_node_point(aes(color= Family),size=0.3)+
  theme_graph(background = "white")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

Algo interesante de ambas descomposiciones, es que las interacciones con signo negativo se centran en los cores más profundos de la red: en la red de saludables se concentran en la subred 38-core, mientras que en la de infectados se concentran en la subred 35-core. Sin embargo, parece que las interacciones negativas son mucho más centrales en la red de infectados que en la de saludables. Además, se puede observar que hay una mayor población de OTUs en los cores más centrales (a partir del 33-core aproximadamente) en la red de saludables contra la de infectados.


# Detección de comunidades y enriquecimiento


```r
back_sld = layout_as_backbone(g_sld_core)
back_inf = layout_as_backbone(g_inf_core)


lou_sld = cluster_louvain(g_sld_core, weights = NULL)
lou_inf = cluster_louvain(g_inf_core, weights = NULL)


infomap_sld = cluster_infomap(g_sld_core)
infomap_inf = cluster_infomap(g_inf_core)

wlkt_sld = walktrap.community(g_sld_core)

wlkt_inf = walktrap.community(g_inf_core)


V(g_sld_core)$louvain = as.factor(lou_sld$membership) 
V(g_inf_core)$louvain = as.factor(lou_inf$membership)
V(g_sld_core)$infomap = as.factor(infomap_sld$membership) 
V(g_inf_core)$infomap = as.factor(infomap_inf$membership)
V(g_sld_core)$wlktrp = as.factor(wlkt_sld$membership) 
V(g_inf_core)$wlktrp = as.factor(wlkt_inf$membership)

V(g_inf_core)$color =  "lightgrey"
V(g_sld_core)$color =  "lightgrey"


max(infomap_inf$modularity)
```

```
## [1] 0.5758855
```

```r
max(infomap_sld$modularity)
```

```
## [1] 0.6117096
```

```r
max(lou_inf$modularity)
```

```
## [1] 0.6162394
```

```r
max(lou_sld$modularity)
```

```
## [1] 0.64039
```

```r
max(wlkt_inf$modularity)
```

```
## [1] 0.5770054
```

```r
max(wlkt_sld$modularity)
```

```
## [1] 0.6168495
```

Se utilizaron tres diferentes algoritmos para la identificación de comunidades en ambas redes: Infomap, Louvain y Walktrap. Los tres obtuvieron valores similares de modularidad, sin embargo, el valor más alto para ambas redes se obtuvo por Louvain. Además, Infomap y Wlktrap identificaban un número mucho más grade de clusters. Debido a lo anterior, se utilizó Louvain para el resto de los análisis.


## Red de saludables


```r
#De igual manera en este plot solo se resaltan los signos de interacción
plot(g_sld_core, layout=back_sld$xy)
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
# Se resalta a qué familia pertenece cada OTUs

ggraph(g_sld_core, layout =  back_sld$xy)+ 
    geom_edge_link0(aes(width=peso,color=signo), show.legend = T, alpha=1)+
  scale_edge_width(range = c(0.01, 0.2))+
  geom_node_point(aes(colour=as.factor(Family), alpha=grado, size=grado))+ 
  theme_graph(background = "white")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
# Se resalta a qué módulo pertenece cada OTUs

ggraph(g_sld_core, layout =  back_sld$xy)+ 
  geom_edge_link0(aes(width=peso,color=signo), show.legend = T, alpha=1)+
  scale_edge_width(range = c(0.01, 0.2))+
  geom_node_point(aes(colour=as.factor(louvain), alpha=grado, size=grado))+ 
theme_graph(background = "white")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```r
# Observamos cuáles son los módulos más prominentes

DF_lou_sld = cbind.data.frame(V(g_sld_core)$name, as.factor(V(g_sld_core)$louvain))

colnames(DF_lou_sld) = c("OTU", "Community")

DF_lou_inf = cbind.data.frame(V(g_inf_core)$name, as.factor(V(g_inf_core)$louvain))

colnames(DF_lou_inf) = c("OTU", "Community")


library(ggplot2)


ggplot(DF_lou_sld) +
 aes(x = Community) +
 geom_bar(fill = "#112446") +
 labs(x = "Módulos", y = "Número de nodos", 
 title = "Tamaño de módulos en Red Saludable") +
 theme_minimal() +
 theme(plot.title = element_text(size = 19L), 
 axis.title.y = element_text(size = 15L), axis.title.x = element_text(size = 15L))
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

```r
ggplot(DF_lou_inf) +
 aes(x = Community) +
 geom_bar(fill = "violet") +
 labs(x = "Módulos", y = "Número de nodos", 
 title = "Tamaño de módulos en Red Infectados") +
 theme_minimal() +
 theme(plot.title = element_text(size = 19L), 
 axis.title.y = element_text(size = 15L), axis.title.x = element_text(size = 15L))
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

```r
ggplot(DF_lou_sld) +
 aes(x = Community) +
 geom_bar(fill = "lightblue") +
  labs(x="Módulo")+
 theme_bw()
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-8-6.png)<!-- -->


## Red de infectados


```r
#interacciones con signo
plot(g_inf_core, layout=back_inf$xy)
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
# Se resalta a qué familia pertenece cada OTUs

ggraph(g_inf_core, layout =  back_inf$xy)+ 
  geom_edge_link0(aes(width=peso,colour=signo), alpha=1)+
  scale_edge_width(range = c(0.01, 0.2))+
 geom_node_point(aes(colour=as.factor(Family), alpha=grado, size=grado))+ 
  theme_graph(background = "white")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

```r
# Se resalta la comunidad de cada OTU

ggraph(g_inf_core, layout =  back_inf$xy)+ 
  geom_edge_link0(aes(width=peso,colour=signo), alpha=1)+
  scale_edge_width(range = c(0.01, 0.2))+
  geom_node_point(aes(colour=as.factor(louvain), alpha=grado, size=grado))+ 
  theme_graph(background = "white")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-9-3.png)<!-- -->



# Densidad de módulos



```r
mi1 <- filter(DF_lou_inf, DF_lou_inf$Community==1)
mi2 <- filter(DF_lou_inf, DF_lou_inf$Community==2)
mi3 <- filter(DF_lou_inf, DF_lou_inf$Community==3)
mi4 <- filter(DF_lou_inf, DF_lou_inf$Community==4)
mi5 <- filter(DF_lou_inf, DF_lou_inf$Community==5)
mi6 <- filter(DF_lou_inf, DF_lou_inf$Community==6)
mi7 <- filter(DF_lou_inf, DF_lou_inf$Community==7)
mi8 <- filter(DF_lou_inf, DF_lou_inf$Community==8)
mi9 <- filter(DF_lou_inf, DF_lou_inf$Community==9)
mi10 <- filter(DF_lou_inf, DF_lou_inf$Community==10)
mi11 <- filter(DF_lou_inf, DF_lou_inf$Community==11)
mi12 <- filter(DF_lou_inf, DF_lou_inf$Community==12)
mi13 <- filter(DF_lou_inf, DF_lou_inf$Community==13)

di1 <- induced_subgraph(g_inf_core, mi1$OTU) %>% gsize() 
di2 <- induced_subgraph(g_inf_core, mi2$OTU) %>% gsize() 
di3 <- induced_subgraph(g_inf_core, mi3$OTU) %>% gsize() 
di4 <- induced_subgraph(g_inf_core, mi4$OTU) %>% gsize() 
di5 <- induced_subgraph(g_inf_core, mi5$OTU) %>% gsize() 
di6 <- induced_subgraph(g_inf_core, mi6$OTU) %>% gsize() 
di7 <- induced_subgraph(g_inf_core, mi7$OTU) %>% gsize() 
di8 <- induced_subgraph(g_inf_core, mi8$OTU) %>% gsize() 
di9 <- induced_subgraph(g_inf_core, mi9$OTU) %>% gsize() 
di10 <- induced_subgraph(g_inf_core, mi10$OTU) %>% gsize() 
di11 <- induced_subgraph(g_inf_core, mi11$OTU) %>% gsize() 
di12 <- induced_subgraph(g_inf_core, mi12$OTU) %>% gsize() 
di13 <- induced_subgraph(g_inf_core, mi13$OTU) %>% gsize() 

Densidades_inf <- c(di1,di2,di3,di4,di5,di6,di7,di8,di9,di10,di11,di12,di13)

size_inf <- table(DF_lou_inf$Community) %>% as.data.frame()
Dens_vs_size_inf <- cbind.data.frame(size_inf,Densidades_inf)
Dens_vs_size_inf
```

```
##    Var1 Freq Densidades_inf
## 1     1  138           1945
## 2     2  113           1484
## 3     3  111            788
## 4     4  244           2879
## 5     5  144           2119
## 6     6    9             22
## 7     7    5              7
## 8     8   41            185
## 9     9  126           2869
## 10   10  189           1865
## 11   11  144           1810
## 12   12    2              1
## 13   13  118            940
```



```r
mi1 <- filter(DF_lou_sld, DF_lou_sld$Community==1)
mi2 <- filter(DF_lou_sld, DF_lou_sld$Community==2)
mi3 <- filter(DF_lou_sld, DF_lou_sld$Community==3)
mi4 <- filter(DF_lou_sld, DF_lou_sld$Community==4)
mi5 <- filter(DF_lou_sld, DF_lou_sld$Community==5)
mi6 <- filter(DF_lou_sld, DF_lou_sld$Community==6)
mi7 <- filter(DF_lou_sld, DF_lou_sld$Community==7)
mi8 <- filter(DF_lou_sld, DF_lou_sld$Community==8)
mi9 <- filter(DF_lou_sld, DF_lou_sld$Community==9)
mi10 <- filter(DF_lou_sld, DF_lou_sld$Community==10)
mi11 <- filter(DF_lou_sld, DF_lou_sld$Community==11)
mi12 <- filter(DF_lou_sld, DF_lou_sld$Community==12)
mi13 <- filter(DF_lou_sld, DF_lou_sld$Community==13)
mi14 <- filter(DF_lou_sld, DF_lou_sld$Community==14)
mi15 <- filter(DF_lou_sld, DF_lou_sld$Community==15)
mi16 <- filter(DF_lou_sld, DF_lou_sld$Community==16)
mi17 <- filter(DF_lou_sld, DF_lou_sld$Community==17)
mi18 <- filter(DF_lou_sld, DF_lou_sld$Community==18)
mi19 <- filter(DF_lou_sld, DF_lou_sld$Community==19)



ds1 <- induced_subgraph(g_sld_core, mi1$OTU) %>% gsize() 
ds2 <- induced_subgraph(g_sld_core, mi2$OTU) %>% gsize() 
ds3 <- induced_subgraph(g_sld_core, mi3$OTU) %>% gsize() 
ds4 <- induced_subgraph(g_sld_core, mi4$OTU) %>% gsize() 
ds5 <- induced_subgraph(g_sld_core, mi5$OTU) %>% gsize() 
ds6 <- induced_subgraph(g_sld_core, mi6$OTU) %>% gsize() 
ds7 <- induced_subgraph(g_sld_core, mi7$OTU) %>% gsize() 
ds8 <- induced_subgraph(g_sld_core, mi8$OTU) %>% gsize() 
ds9 <- induced_subgraph(g_sld_core, mi9$OTU) %>% gsize() 
ds10 <- induced_subgraph(g_sld_core, mi10$OTU) %>% gsize() 
ds11 <- induced_subgraph(g_sld_core, mi11$OTU) %>% gsize() 
ds12 <- induced_subgraph(g_sld_core, mi12$OTU) %>% gsize() 
ds13 <- induced_subgraph(g_sld_core, mi13$OTU) %>% gsize() 
ds14 <- induced_subgraph(g_sld_core, mi14$OTU) %>% gsize() 
ds15 <- induced_subgraph(g_sld_core, mi15$OTU) %>% gsize() 
ds16 <- induced_subgraph(g_sld_core, mi16$OTU) %>% gsize() 
ds17 <- induced_subgraph(g_sld_core, mi17$OTU) %>% gsize() 
ds18 <- induced_subgraph(g_sld_core, mi18$OTU) %>% gsize() 
ds19 <- induced_subgraph(g_sld_core, mi19$OTU) %>% gsize() 


Densidades_sld <- c(ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,ds9,ds10,ds11,ds12,ds13,ds14,ds15,ds16,ds17,ds18,ds19)

size_sld <- table(DF_lou_sld$Community) %>% as.data.frame()
Dens_vs_size_sld <- cbind.data.frame(size_sld,Densidades_sld)
Dens_vs_size_sld
```

```
##    Var1 Freq Densidades_sld
## 1     1   81           1754
## 2     2  258           4686
## 3     3  233           3210
## 4     4  123           1004
## 5     5  304           6214
## 6     6   37            221
## 7     7   23            160
## 8     8   35            169
## 9     9   30             67
## 10   10   35            120
## 11   11   32             89
## 12   12   30            105
## 13   13   42            454
## 14   14   24             77
## 15   15   98            376
## 16   16   36            108
## 17   17    5              8
## 18   18   15             29
## 19   19    2              1
```


```r
ggplot(Dens_vs_size_inf) +
 aes(x =log10(Freq), y = log10(Densidades_inf)) +
 geom_point(shape = "circle", size = 3, 
 colour = "#112446") +
 theme_minimal()
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
ggplot(Dens_vs_size_sld) +
 aes(x = log10(Freq), y = log10(Densidades_sld)) +
 geom_point(shape = "circle", size = 3, 
 colour = "green") +
 theme_minimal()
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-12-2.png)<!-- -->










# Enriquecimiento



```r
# data wrangling. Queremos dos DFs con la informacion de las comunidades, OTUs y Taxonomias

tax_core_inf =  filter(tax_inf, rownames(tax_inf) %in% DF_lou_inf$OTU)
tax_core_sld =  filter(tax_sld, rownames(tax_sld) %in% DF_lou_sld$OTU)

rownames(tax_core_inf) = NULL
rownames(tax_core_sld) = NULL


DF_inf_taxmod = cbind.data.frame(DF_lou_inf, tax_core_inf)
DF_sld_taxmod = cbind.data.frame(DF_lou_sld, tax_core_sld)
```





```r
# Creamos funcion que haga la distribucion hipergeometrica y calcule valor de p
# de una prueba exacta de fisher de una sola cola

# Le damos los datos necesarios para crear una tabla de contingencia de 2x2

hypergeometric = function(inmodulo,outmodulo,ntax,taxinmodulo){

m = inmodulo # Numero de OTUs en modulo definido
n = outmodulo# OTUS fuera de módulo definido
k = ntax# Numero de taxonomia definida
x = taxinmodulo:m # OTUs de taxonomia definida dentro de módulo

probabilities <- dhyper(x, m, n, k, log = FALSE)
p.val=sum(probabilities)

return(p.val)
}
```




```r
Families = levels(as.factor(tax_inf$Family))


enrichment_OTU = function(taxtype, DFred, modulo){
  
  preinmodulo = filter(DFred, DFred$Community==modulo)$OTU
  inmodulo = length(preinmodulo)
  
  outmodulo = length(DFred$OTU) - inmodulo
  
  prentax = filter(DFred, DFred$Family==taxtype)$OTU
  ntax = length(prentax)
  
  taxinmodulo = length(intersect(preinmodulo,prentax))
  
 p.val = hypergeometric(inmodulo, outmodulo, ntax, taxinmodulo)
 
 return(p.val)
  
}
```



```r
modulos_sld =  length(levels(as.factor(DF_sld_taxmod$Community)))
modulos_inf = length(levels(as.factor(DF_inf_taxmod$Community)))



sld_enrichment_mod_loop = matrix(NA, ncol = modulos_sld, nrow = length(Families))

inf_enrichment_mod_loop = matrix(NA, ncol = modulos_inf, nrow = length(Families))

for (i in seq_along(Families)) {
  sld_enrichment_mod_loop[i] = Families[i]
  
  for (j in seq_along(1:modulos_sld)) {
    
    sld_enrichment_mod_loop[i,j] = enrichment_OTU(Families[i], DF_sld_taxmod, j)
    
    
    
  }

  
  
}

sld_enrichment_mod = cbind.data.frame(Families, sld_enrichment_mod_loop)


for (i in seq_along(Families)) {
  inf_enrichment_mod_loop[i] = Families[i]
  
  for (j in seq_along(1:modulos_inf)) {
    
    inf_enrichment_mod_loop[i,j] = enrichment_OTU(Families[i], DF_inf_taxmod, j)
  
  }
  
  
}

inf_enrichment_mod = cbind.data.frame(Families, inf_enrichment_mod_loop)
```


```r
matrix_inf_num = inf_enrichment_mod[,-1] %>% as.matrix() %>% as.numeric() 
matrix_inf = matrix(matrix_inf_num, ncol = 13, nrow = 33) 


df_p_inf = as.data.frame(matrix_inf[,2])




pheatmap((matrix_inf), cluster_cols = F, cluster_rows = F, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="OrRd")))(5),legend=T,breaks = c(0.01, 0.02,0.03,0.04,0.05,0.06), 
         cellwidth = 15, cellheight = 15, border_color = "slategrey", labels_row =inf_enrichment_mod$Families, labels_col = colnames(inf_enrichment_mod[-1]),
         angle_col = 0, , legend_breaks = c(0.06,0.05,0.04,0.03,0.02,0.01),legend_labels = c("P > 0.05","0.05","0.04","0.03","0.02","0.01"))
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-17-1.png)<!-- -->




```r
matrix_sld = sld_enrichment_mod[,-1] %>% as.matrix() %>% as.numeric() 
matrix_sld = matrix(matrix_sld, ncol = 19, nrow = 33) 


pprueba = p.adjust(matrix_sld[3,], method = "fdr")

pheatmap((matrix_sld), cluster_cols = F, cluster_rows = F, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="GnBu")))(5),legend=T,breaks = c(0.01, 0.02,0.03,0.04,0.05,0.06), 
         cellwidth = 15, cellheight = 15, border_color = "slategrey", labels_row =sld_enrichment_mod$Families, labels_col = colnames(sld_enrichment_mod[-1]),
         angle_col = 0, , legend_breaks = c(0.06,0.05,0.04,0.03,0.02,0.01),legend_labels = c("P > 0.05","0.05","0.04","0.03","0.02","0.01"), main = "Familias sobrerepresentadas en módulos de red ")
```

![](Comunidades_y_kcore_files/figure-html/unnamed-chunk-18-1.png)<!-- -->



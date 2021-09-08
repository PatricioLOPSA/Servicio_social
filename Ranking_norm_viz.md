---
title: "Redes de relevancia"
author: "Patricio López Sánchez"
date: "8/12/2021"
output:
  html_document:
    keep_md: yes
---




# Rankeando aristas en redes


```r
library(dplyr)
library(ggplot2)
library(pheatmap)
library(igraph)
library(vroom)
```



```r
net_inf = vroom("prs_infect.csv", col_names = TRUE)
```

```
## New names:
## * `` -> ...1
```

```
## Rows: 990528 Columns: 6
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## dbl (6): ...1, row, column, cor, p, p.adj
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
net_inf = net_inf[,-1]
head(net_inf)
```

```
## # A tibble: 6 × 5
##      row column      cor           p     p.adj
##    <dbl>  <dbl>    <dbl>       <dbl>     <dbl>
## 1 198059 985239 0.480    0.000000332 0.0000274
## 2 198059 840914 0.212    0.0321      0.160    
## 3 985239 840914 0.0657   0.512       0.729    
## 4 198059 174924 0.000205 0.998       0.999    
## 5 985239 174924 0.0878   0.380       0.625    
## 6 840914 174924 0.0970   0.332       0.583
```

```r
net_sld = vroom("prs_salud.csv", col_names = T)
```

```
## New names:
## * `` -> ...1
```

```
## Rows: 1076778 Columns: 6
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## dbl (6): ...1, row, column, cor, p, p.adj
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
net_sld = net_sld[-1]
head(net_sld)
```

```
## # A tibble: 6 × 5
##      row column      cor        p         p.adj
##    <dbl>  <dbl>    <dbl>    <dbl>         <dbl>
## 1 198059 985239  0.351   2.24e-11 0.00000000211
## 2 198059 840914  0.0547  3.13e- 1 0.532        
## 3 985239 840914  0.00175 9.74e- 1 0.988        
## 4 198059 174924  0.0443  4.13e- 1 0.626        
## 5 985239 174924  0.120   2.60e- 2 0.107        
## 6 840914 174924 -0.0650  2.30e- 1 0.442
```



```r
net_inf_fdr = filter(net_inf, net_inf$p.adj < 0.01)
net_inf_fdr2 = filter(net_inf, net_inf$p.adj < 0.001)
```


Podemos usar el número de aristas de la red de infectados como punto de corte para la red de saludables, que manera que tengamos el mismo número de aristas en ambas redes



```r
net_sld_ord = net_sld[order(-abs(net_sld$cor)),]

net_sld_rank = net_sld_ord[1:55857, 1:5]
net_sld_rank2 = net_sld_ord[1:27747, 1:5]

# Checamos los valores de p que tienen cada interacción rankeada

max(net_sld_rank$p.adj)
```

```
## [1] 0.0003949434
```

```r
max(net_sld_rank2$p.adj)
```

```
## [1] 5.327337e-06
```

```r
# En ambos casos, todas las interacciones tienen fdr < 0.001.
```



Exploramos las distribuciones de ambas redes.


```r
plot(density(net_inf_fdr$cor))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
plot(density(net_inf_fdr2$cor))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
plot(density(net_sld_rank$cor))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
plot(density(net_sld_rank2$cor))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-5-4.png)<!-- -->


## Normalización de correlación


Se normalizarán los valores de correlación de manera que el valor más chico sea 0, y el más grande 1. Para lograr esto se utilizará una normalización min-max, y se hará de manera independiente las correlaciones negativas de las positivas.


```r
minmax <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}
```



Separamos correlaciones positivas y negativas


```r
inf_pos = filter(net_inf_fdr, net_inf_fdr$cor >= 0)
inf_neg = filter(net_inf_fdr, net_inf_fdr$cor < 0)

sld_pos = filter(net_sld_rank, net_sld_rank$cor >= 0)
sld_neg = filter(net_sld_rank, net_sld_rank$cor < 0)

## Para el caso p < 0.001

inf_pos2 = filter(net_inf_fdr2, net_inf_fdr2$cor >= 0)
inf_neg2 = filter(net_inf_fdr2, net_inf_fdr2$cor < 0)

sld_pos2 = filter(net_sld_rank2, net_sld_rank2$cor >= 0)
sld_neg2 = filter(net_sld_rank2, net_sld_rank2$cor < 0)
```


Normalizamos las correlaciones y unimos los data frames



```r
inf_pos$norm = minmax(inf_pos$cor)
inf_neg$norm = minmax(abs(inf_neg$cor)) * -1

sld_pos$norm = minmax(sld_pos$cor)
sld_neg$norm = minmax(abs(sld_neg$cor)) * -1


## Para p < 0.001

inf_pos2$norm = minmax(inf_pos2$cor)
inf_neg2$norm = minmax(abs(inf_neg2$cor)) * -1

sld_pos2$norm = minmax(sld_pos2$cor)
sld_neg2$norm = minmax(abs(sld_neg2$cor)) * -1

# Unimos los Data frames

net_inf_norm = rbind.data.frame(inf_pos, inf_neg)
net_sld_norm = rbind.data.frame(sld_pos, sld_neg)

net_inf_norm2 = rbind.data.frame(inf_pos2, inf_neg2)
net_sld_norm2 = rbind.data.frame(sld_pos2, sld_neg2)
```


Checamos density plots de las nuevas redes.


```r
plot(density(net_inf_norm$norm))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
plot(density(net_sld_norm$norm))
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-9-2.png)<!-- -->


El siguiente paso es crear el objeto de *igraph* y corregir para nodos borrados durante los pasos de filtrado de correlación y rankeo.



```r
g_inf = graph_from_data_frame(net_inf_norm, directed = F)
g_sld = graph_from_data_frame(net_sld_norm, directed = F)

v_inf = V(g_inf)$name
v_sld = V(g_sld)$name

g_cominf = graph_from_data_frame(net_inf, directed = F)
g_comsld = graph_from_data_frame(net_sld, directed = F)

v_cominf = V(g_cominf)$name
v_comsld = V(g_comsld)$name

"%ni%" <- Negate("%in%")

v_missinf = subset(v_cominf, v_cominf %ni% v_inf) 
v_misssld = subset(v_comsld, v_comsld %ni% v_sld) 


g_inf_com = g_inf + add_vertices(g_inf, length(v_missinf), attr=list(name=v_missinf))
g_sld_com = g_sld + add_vertices(g_sld, length(v_misssld), attr = list(name=v_misssld))
```




# Distribución del grado y comparando el rankeo del nodos según el grado


Checamos cómo se distribuye el grado del nodo en ambas condiciones.


```r
deg_inf = degree(g_inf_com) %>% as.data.frame()
deg_sld = degree(g_sld_com) %>% as.data.frame()



ggplot(deg_inf) +
 aes(x = .) +
 geom_histogram(bins = 100L, fill = "#6EC3D2") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
ggplot(deg_sld) +
 aes(x = .) +
 geom_histogram(bins = 100L, fill = "pink") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-11-2.png)<!-- -->


Repetimos los mismos pasos para la red p < 0.001


```r
g_inf2 = graph_from_data_frame(net_inf_norm2, directed = F)
g_sld2 = graph_from_data_frame(net_sld_norm2, directed = F)

v_inf2 = V(g_inf2)$name
v_sld2 = V(g_sld2)$name


v_missinf2 = subset(v_cominf, v_cominf %ni% v_inf2) 
v_misssld2 = subset(v_comsld, v_comsld %ni% v_sld2) 


g_inf_com2 = g_inf2 + add_vertices(g_inf2, length(v_missinf2), attr=list(name=v_missinf2))
g_sld_com2 = g_sld2 + add_vertices(g_sld2, length(v_misssld2), attr = list(name=v_misssld2))
```



```r
deg_inf2 = degree(g_inf_com2) %>% as.data.frame()
deg_sld2 = degree(g_sld_com2) %>% as.data.frame()



ggplot(deg_inf2) +
 aes(x = .) +
 geom_histogram(bins = 100L, fill = "#6EC3D2") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
ggplot(deg_sld2) +
 aes(x = .) +
 geom_histogram(bins = 100L, fill = "pink") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-13-2.png)<!-- -->



```r
#Rankeamos grado para infectados

V(g_inf_com)$degree = degree(g_inf_com)
df_deg_inf = cbind.data.frame(V(g_inf_com)$name, V(g_inf_com)$degree)
deg_inf_ord = df_deg_inf[order(-df_deg_inf$`V(g_inf_com)$degree`),]



#Rankeamos grado para saludables
V(g_sld_com)$degree = degree(g_sld_com)

df_deg_sld = cbind.data.frame(V(g_sld_com)$name, V(g_sld_com)$degree)
deg_sld_ord = df_deg_sld[order(-df_deg_sld$`V(g_sld_com)$degree`),]


head(deg_sld_ord)
```

```
##     V(g_sld_com)$name V(g_sld_com)$degree
## 202            295485                 324
## 283            324214                 324
## 183            528652                 315
## 347            531888                 315
## 211            369014                 303
## 144            287978                 299
```

```r
head(deg_inf_ord)
```

```
##     V(g_inf_com)$name V(g_inf_com)$degree
## 637            519398                 256
## 343            367889                 250
## 30             591285                 243
## 340            575844                 243
## 294            345801                 233
## 349            369379                 231
```

Definimos una función para calcular el índice de Jaccard


```r
jaccard = function(A,B) {
  
  int = length(intersect(A,B))
  uni = length(union(A,B))
  jac = int/uni
  return(jac)
}
```



```r
#JI para todos los nodos de la red más pequeña

jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:1408], 
        deg_inf_ord$`V(g_inf_com)$name`[1:1408])
```

```
## [1] 0.7490683
```

```r
#JI para top 10 nodos
jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:10], 
        deg_inf_ord$`V(g_inf_com)$name`[1:10])
```

```
## [1] 0.05263158
```

```r
#JI para top 20 nodos
jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:20], 
        deg_inf_ord$`V(g_inf_com)$name`[1:20])
```

```
## [1] 0.05263158
```

```r
#JI para top 50 nodos
jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:50], 
        deg_inf_ord$`V(g_inf_com)$name`[1:50])
```

```
## [1] 0.1904762
```

```r
#JI para top 100 nodos
jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:100], 
        deg_inf_ord$`V(g_inf_com)$name`[1:100])
```

```
## [1] 0.1904762
```

```r
#JI para top 500 nodos
jaccard(deg_sld_ord$`V(g_sld_com)$name`[1:500], 
        deg_inf_ord$`V(g_inf_com)$name`[1:500])
```

```
## [1] 0.3297872
```



## Netviz



```r
com_inf = cluster_louvain(g_inf_com, weights = abs(E(g_inf_com)$norm_1))
V(g_inf_com)$community = com_inf$membership
V(g_inf_com)$label = NA
V(g_inf_com)$size = degree(g_inf_com, v=V(g_inf_com)) %>% minmax() *7 %>% +.2
V(g_inf_com)$frame.color = "white"
V(g_inf_com)$color = V(g_inf_com)$community  #"lightblue3"

E(g_inf_com)$width = abs(E(g_inf_com)$norm_1) / 2.5
E(g_inf_com)$color = ifelse(E(g_inf_com)$norm_1 >0, "darkseagreen3","lightpink")
E(g_inf_com)$curved = 0

l <- layout_with_fr(g_inf_com, niter = 10000, weights = E(g_inf_com)$norm_1)



plot(g_inf_com, rescale = T, layout = l*1)
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-17-1.png)<!-- -->



```r
com_sld = cluster_louvain(g_sld_com, weights = abs(E(g_sld_com)$norm_1))
V(g_sld_com)$community = com_sld$membership
V(g_sld_com)$label = NA
V(g_sld_com)$size = degree(g_sld_com, v=V(g_sld_com)) %>% minmax() *7 %>% +.2
V(g_sld_com)$frame.color = "white"
V(g_sld_com)$color =  V(g_sld_com)$community #"lightblue3"
#
E(g_sld_com)$width = abs(E(g_sld_com)$norm_1) / 2.5
E(g_sld_com)$color = ifelse(E(g_sld_com)$norm_1 >0, "darkseagreen3","lightpink")
E(g_sld_com)$curved = 0

ls <- layout_with_fr(g_sld_com, niter = 50 ,weights = E(g_sld_com)$norm_1)



plot(g_sld_com, rescale = T, layout = ls)
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-18-1.png)<!-- -->




```r
com_inf2 = cluster_louvain(g_inf_com2, weights = abs(E(g_inf_com2)$norm_1))
V(g_inf_com2)$community = com_inf2$membership
V(g_inf_com2)$label = NA
V(g_inf_com2)$size = degree(g_inf_com2, v=V(g_inf_com2)) %>% minmax() *7 %>% +.2
V(g_inf_com2)$frame.color = "white"
V(g_inf_com2)$color = V(g_inf_com2)$community  #"lightblue3"

E(g_inf_com2)$width = abs(E(g_inf_com2)$norm_1) / 2.5
E(g_inf_com2)$color = ifelse(E(g_inf_com2)$norm_1 >0, "darkseagreen3","lightpink")
E(g_inf_com2)$curved = 0

l2 <- layout_with_fr(g_inf_com2, niter = 10000, weights = E(g_inf_com2)$norm_1)




plot(g_inf_com2, rescale = T, layout = l2)
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-19-1.png)<!-- -->



```r
com_sld2 = cluster_louvain(g_sld_com2, weights = abs(E(g_sld_com2)$norm_1))
V(g_sld_com2)$community = com_sld2$membership
V(g_sld_com2)$label = NA
V(g_sld_com2)$size = degree(g_sld_com2, v=V(g_sld_com2)) %>% minmax() *7 %>% +.2
V(g_sld_com2)$frame.color = "white"
V(g_sld_com2)$color =  V(g_sld_com2)$community #"lightblue3"
#
E(g_sld_com2)$width = abs(E(g_sld_com2)$norm_1) / 2.5
E(g_sld_com2)$color = ifelse(E(g_sld_com2)$norm_1 >0, "darkseagreen3","lightpink")
E(g_sld_com2)$curved = 0

ls2 <- layout_with_fr(g_sld_com2, niter =  50 ,weights = E(g_sld_com2)$norm_1)


plot(g_sld_com2, rescale = T, layout = ls2)
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-20-1.png)<!-- -->



# Redes separadas con aristas negativas y positivas

Observamos individualmente las subredes con interacciones solo positivas o negativas



```r
g_inf_neg <- graph_from_data_frame(inf_neg2, directed = F)
g_inf_pos <- graph_from_data_frame(inf_pos2, directed = F)

g_sld_neg <- graph_from_data_frame(sld_neg2, directed = F)
g_sld_pos <- graph_from_data_frame(sld_pos2, directed = F)
```


Calculamos la fuerza y el grado de cada nodo para cada red.



```r
#Calculamos grado
V(g_inf_neg)$degree = degree(g_inf_neg)
V(g_inf_pos)$degree = degree(g_inf_pos)

V(g_sld_neg)$degree = degree(g_sld_neg)
V(g_sld_pos)$degree = degree(g_sld_pos)

V(g_inf_com2)$degree = degree(g_inf_com2)
V(g_sld_com2)$degree = degree(g_sld_com2)

#Calculamos fuerza
V(g_inf_neg)$str = strength(g_inf_neg, weights = abs(E(g_inf_neg)$norm))
V(g_inf_pos)$str = strength(g_inf_pos, weights = E(g_inf_pos)$norm)

V(g_sld_neg)$str = strength(g_sld_neg, weights = abs(E(g_sld_neg)$norm))
V(g_sld_pos)$str = strength(g_sld_pos, weights = E(g_sld_pos)$norm)

V(g_inf_com2)$str = strength(g_inf_com2, weights = abs(E(g_inf_com2)$norm_1))
V(g_sld_com2)$str = strength(g_sld_com2, weights = abs(E(g_sld_com2)$norm_1))


#Calculamos fuerza promedio para cada nodo (s_i = str_i/deg_i)

V(g_inf_neg)$strnorm = V(g_inf_neg)$str/V(g_inf_neg)$degree
V(g_inf_pos)$strnorm = V(g_inf_pos)$str/V(g_inf_pos)$degree

V(g_sld_neg)$strnorm = V(g_sld_neg)$str/V(g_sld_neg)$degree
V(g_sld_pos)$strnorm = V(g_sld_pos)$str/V(g_sld_pos)$degree

V(g_inf_com2)$strnorm = V(g_inf_com2)$str/V(g_inf_com2)$degree
V(g_sld_com2)$strnorm = V(g_sld_com2)$str/V(g_sld_com2)$degree
```



Checamos plot de grado vs fuerza promedio



```r
deg_vs_str_infneg = cbind.data.frame(V(g_inf_neg)$degree, V(g_inf_neg)$strnorm)
deg_vs_str_infpos = cbind.data.frame(V(g_inf_pos)$degree, V(g_inf_pos)$strnorm)
deg_vs_str_infcom = cbind.data.frame(V(g_inf_com2)$degree, V(g_inf_com2)$strnorm)


deg_vs_str_sldneg = cbind.data.frame(V(g_sld_neg)$degree, V(g_sld_neg)$strnorm)
deg_vs_str_sldpos = cbind.data.frame(V(g_sld_pos)$degree, V(g_sld_pos)$strnorm)
deg_vs_str_sldcom = cbind.data.frame(V(g_sld_com2)$degree, V(g_sld_com2)$strnorm)
```




```r
ggplot(deg_vs_str_infneg) +
 aes(x = `V(g_inf_neg)$degree`, y = `V(g_inf_neg)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "#C9092C") + geom_vline(xintercept = mean( deg_vs_str_infneg$`V(g_inf_neg)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_infneg$`V(g_inf_neg)$strnorm`), col = "black") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

```r
ggplot(deg_vs_str_infpos) +
 aes(x = `V(g_inf_pos)$degree`, y = `V(g_inf_pos)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "green") + geom_vline(xintercept = mean( deg_vs_str_infpos$`V(g_inf_pos)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_infpos$`V(g_inf_pos)$strnorm`), col = "black") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-24-2.png)<!-- -->

```r
ggplot(deg_vs_str_infcom) +
 aes(x = `V(g_inf_com2)$degree`, y = `V(g_inf_com2)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "steelblue3") + geom_vline(xintercept = mean( deg_vs_str_infcom$`V(g_inf_com2)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_infcom$`V(g_inf_com2)$strnorm`, na.rm = T), col = "black") +
 theme_minimal()
```

```
## Warning: Removed 24 rows containing missing values (geom_point).
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-24-3.png)<!-- -->




```r
ggplot(deg_vs_str_sldneg) +
 aes(x = `V(g_sld_neg)$degree`, y = `V(g_sld_neg)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "#C9092C") + geom_vline(xintercept = mean( deg_vs_str_sldneg$`V(g_sld_neg)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_sldneg$`V(g_sld_neg)$strnorm`), col = "black") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
ggplot(deg_vs_str_sldpos) +
 aes(x = `V(g_sld_pos)$degree`, y = `V(g_sld_pos)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "green") + geom_vline(xintercept = mean( deg_vs_str_sldpos$`V(g_sld_pos)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_sldpos$`V(g_sld_pos)$strnorm`), col = "black") +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

```r
ggplot(deg_vs_str_sldcom) +
 aes(x = `V(g_sld_com2)$degree`, y = `V(g_sld_com2)$strnorm`) +
 geom_point(shape = "circle", 
 size = 1.5, colour = "plum2") + geom_vline(xintercept = mean( deg_vs_str_sldcom$`V(g_sld_com2)$degree`), col = "black") +
  geom_hline(yintercept = mean(deg_vs_str_sldcom$`V(g_sld_com2)$strnorm`, na.rm = T), col = "black") +
 theme_minimal()
```

```
## Warning: Removed 25 rows containing missing values (geom_point).
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-25-3.png)<!-- -->


Checamos grado de nodo vs. rank.



```r
deg_rank_inf = cbind.data.frame(V(g_inf_com2)$name, V(g_inf_com2)$degree)
deg_rank_inf =  deg_rank_inf[order(-deg_rank_inf$`V(g_inf_com2)$degree`),]
deg_rank_inf$rank = rank(-deg_rank_inf$`V(g_inf_com2)$degree`, ties.method = "min")



deg_rank_sld = cbind.data.frame(V(g_sld_com2)$name, V(g_sld_com2)$degree)
deg_rank_sld =  deg_rank_sld[order(-deg_rank_sld$`V(g_sld_com2)$degree`),]
deg_rank_sld$rank = rank(-deg_rank_sld$`V(g_sld_com2)$degree`, ties.method = "min")



ggplot() + 
geom_line(data=deg_rank_inf, aes(x=log10(`V(g_inf_com2)$degree`), y=log10(rank)), color='midnightblue') + 
geom_line(data=deg_rank_sld, aes(x=log10(`V(g_sld_com2)$degree`), y=log10(rank)), color='tomato') + theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-26-1.png)<!-- -->


Checamos fuerza norm v.s. rank



```r
str_rank_inf = cbind.data.frame(V(g_inf_com2)$name, V(g_inf_com2)$strnorm)
str_rank_inf =  str_rank_inf[order(-str_rank_inf$`V(g_inf_com2)$strnorm`),]
str_rank_inf$rank = rank(-str_rank_inf$`V(g_inf_com2)$strnorm`, ties.method = "min")



str_rank_sld = cbind.data.frame(V(g_sld_com2)$name, V(g_sld_com2)$strnorm)
str_rank_sld =  str_rank_sld[order(-str_rank_sld$`V(g_sld_com2)$strnorm`),]
str_rank_sld$rank = rank(-str_rank_sld$`V(g_sld_com2)$strnorm`, ties.method = "min")



ggplot() + 
geom_line(data=str_rank_inf, aes(x=log10(`V(g_inf_com2)$strnorm`), y=log10(rank)), color='midnightblue') + 
geom_line(data=str_rank_sld, aes(x=log10(`V(g_sld_com2)$strnorm`), y=log10(rank)), color='tomato') + theme_minimal()
```

```
## Warning: Removed 24 row(s) containing missing values (geom_path).
```

```
## Warning: Removed 25 row(s) containing missing values (geom_path).
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-27-1.png)<!-- -->


Checamos cómo cambiaron los rankings de los nodos comunes para ambas condiciones 



```r
itsct_sld = filter(deg_rank_sld, deg_rank_sld$`V(g_sld_com2)$name` %in% deg_rank_inf$`V(g_inf_com2)$name`) 
itsct_sld_ord = itsct_sld[order(itsct_sld$`V(g_sld_com2)$name`),]



itsct_inf = filter(deg_rank_inf, deg_rank_inf$`V(g_inf_com2)$name` %in% deg_rank_sld$`V(g_sld_com2)$name`)
itsct_inf_ord = itsct_inf[order(itsct_inf$`V(g_inf_com2)$name`),]


ranks_compare = itsct_sld_ord
ranks_compare$rank_inf = itsct_inf_ord$rank
ranks_compare$grado_inf = itsct_inf_ord$`V(g_inf_com2)$degree`

colnames(ranks_compare) = c("OTU", "grado_sld", "rank_salud", "rank_inf", "grado_inf")
```



```r
ggplot(ranks_compare) +
 aes(x =-log10(rank_inf), y = -log10(rank_salud)) +
 geom_point(shape = "circle", size = 1, 
 colour = "#FF269E") +
geom_abline(intercept = 0, slope = 1) +
 theme_minimal()
```

![](Ranking_norm_viz_files/figure-html/unnamed-chunk-29-1.png)<!-- -->



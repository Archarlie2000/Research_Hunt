---
title: "Validation"
date: '2023-08-12'
output:
  html_document:  
    keep_md: true
    toc: true
    code_folding: hide
    fig_height: 6
    fig_width: 12
    fig_align: 'center'
---
---





```r
library(tidyverse)
library(mosaic)
library(pander)
library(ggplot2)
library(ggfortify)
```

## Introduction

To discover the Differential Gen Expression (DIGS) before and after implantation arrest, we collected NGS results on 9 mink, 3 at before the diapause, 3 different the diapause, and 3 after diapause.


### Remove outliers from 9 mink



```r
E1 <- read.csv('13E-1_S26_L007.csv', header = TRUE)  %>% 
  mutate(phase = "E1")
E2 <- read.csv('13E-2_S29_L007.csv', header = TRUE)  %>% 
  mutate(phase = "E2")
E3 <- read.csv('13E-3_S32_L007.csv', header = TRUE)  %>% 
  mutate(phase = "E3")

D1 <- read.csv('13D-1_S27_L007.csv', header = TRUE) %>% 
  mutate(phase = "D1")
D2 <- read.csv('13D-2_S30_L007.csv', header = TRUE) %>% 
  mutate(phase = "D2")
D3 <- read.csv('13D-3_S33_L007.csv', header = TRUE) %>% 
  mutate(phase = "D3")

P1 <- read.csv('13P-1_S28_L007.csv', header = TRUE)  %>% 
  mutate(phase = "P1")
P2 <- read.csv('13P-2_S31_L007.csv', header = TRUE)  %>% 
  mutate(phase = "P2")
P3 <- read.csv('13P-3_S34_L007.csv', header = TRUE)  %>% 
  mutate(phase = "P3")


calculate_rank_sums <- function(A, B, C) {
  
  data <- cbind(A, B, C)

  # 1. Calculate the difference
  diff_A <- abs(A - (B + C) / 2)
  diff_B <- abs(B - (A + C) / 2)
  diff_C <- abs(C - (A + B) / 2)
  diff_data <- cbind(diff_A, diff_B, diff_C)

  # 2. Rank the differences
  rank_data <- t(apply(diff_data, 1, function(x) rank(-x)))

  return(colSums(rank_data))
}


results_TPM <- rbind(
  E = calculate_rank_sums(E1$TPM, E2$TPM, E3$TPM),
  D = calculate_rank_sums(D1$TPM, D2$TPM, D3$TPM),
  P = calculate_rank_sums(P1$TPM, P2$TPM, P3$TPM)
)

results_FPKM <- rbind(
  E = calculate_rank_sums(E1$FPKM, E2$FPKM, E3$FPKM),
  D = calculate_rank_sums(D1$FPKM, D2$FPKM, D3$FPKM),
  P = calculate_rank_sums(P1$FPKM, P2$FPKM, P3$FPKM)
)

pander(results_TPM, caption = "Results for TPM")
```


-----------------------------------
 &nbsp;   diff_A   diff_B   diff_C 
-------- -------- -------- --------
 **E**    51922    47426    47809  

 **D**    51803    47757    47596  

 **P**    50909    43437    52811  
-----------------------------------

Table: Results for TPM

```r
pander(results_FPKM, caption = "Results for FPKM")
```


-----------------------------------
 &nbsp;   diff_A   diff_B   diff_C 
-------- -------- -------- --------
 **E**    50876    47793    48488  

 **D**    52134    46987    48036  

 **P**    50812    41849    54496  
-----------------------------------

Table: Results for FPKM



## Including Plots



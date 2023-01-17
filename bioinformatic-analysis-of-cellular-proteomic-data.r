---
title: 'Bioinformatic Analysis of Cellular Proteomic Data'
subtitle: 'STA130 Course Project'
author: "Mahe Chen and Matthew Yu"
date: "December 8, 2022"
output:
  beamer_presentation:
    theme: AnnArbor
    colortheme: orchid
    fonttheme: structurebold
    slide_level: 2
    includes:
      in_header: preamble.tex
  slidy_presentation: default
classoption: aspectratio=169
fontsize: 10pt
urlcolor: blue
---

```{r, echo=TRUE, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE}
knitr::opts_chunk$set(eval=TRUE, include=TRUE, echo=TRUE, message=FALSE, warning=FALSE)
library(tidyverse)
knitr::opts_chunk$set(dev='pdf') #Default file format
#This chunk of code allows the size of printed code+output to be controlled using latex text size names
def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\n \\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})

color_block = function(color) {
  function(x, options) sprintf('\\color{%s}\\begin{verbatim}%s\\end{verbatim}',
                               color, x)
}
knitr::knit_hooks$set(message = color_block('red'))
knitr::knit_hooks$set(warning = color_block('red'))
knitr::knit_hooks$set(error = color_block('red'))
```

## Project Description:
- The purpose of this project is to apply statistics in bioinformatic research work to help fight cancer!

- Our data is based on advances of single cell analysis in the fields of Flow Cytometry and cellular proteomic processes in the fields of Mass Spectrometry.

- Based on the above analysis, we will measure the multivariate landscape of proteomic activity for a single cell in any experimental condition for any cell type at scale.

- We will explore how to direct deleterious cellular states to transition into non-deleterious states by analyzing the typical cellular homeostatis of healthy and deliterious cells and the phenotypical transformation of cellular proteomic homeostatsis over time in response to different experimental conditions.

## Our Plan

- By recognizing correlations and patterns between protein makeup(most specifically our outcome proteins) of a cell under some condition, we can classify the cell type and label the conditions as a possible causation

- Analyze protein levels of cells over time and under different treatment conditions

- Protein structure data to recognize patterns of healthy cell states turning deleterious, allow us to interfere with treatment as soon as possible

- Measure protein structure with Mass Spectrometry and Flow Cytometry
Each cell, under some condition and at some time, measures levels of 22 AP-1 transcription factors and 4 outcome phenotype proteins.

- Destroy cells to take measurements

- Split cells into groups and take measurements from each group

## Project Information
- The grid below illustrates what melanoma cells(which our skin cancer cells can spread to the rest of the body) can appear as during their differentiation states, meaning they are changing their function type from stem cell to specific cell.
```{r, echo=FALSE}
tab <- matrix(c("Low","Low","Low","High","Low","High","High","High","High","High","High","Low","High","High","Low","Low"), ncol=4, byrow=TRUE)
colnames(tab) <- c('MiTFg','NGFR','SOX10','AXL')
rownames(tab) <- c('Undifferentiated', 'Neural crest-like','Transitory','Melanocytic')
tab <- as.table(tab)
tab
```
- Some background information regarding our data set:
  - We are going to narrow our focus to the 4 phenotype protein levels
  - Measurements are taken at time intervals following administration of treatment
    - 0.5, 2, 6, 15, 24, 72, 120h
    - Drugs include Vem or Vem and Tram
    - Doses measured in micrometress

## Project Information
- Cellular Phenotype
  - Undifferentiated are like stem cells, these tissues do not have a specialised function yet (Not yet a “mature” cell). Melanoma cells that appear as this are challenging to diagnose and rapidly grow and spread
  - Neural crest-like melanoma cells behave very similarly to neural crest cells(which are normal healthy cells used in body regeneration) during their early stages of spreading and invasiveness, and thus can be challenging to detect. Derived from various tissue or stem cells, unlike neural crest cells.
  - Transitory cells - critical white blood cell (immune cell) that targets and kills antigens and cancer cells
  - Melanocytic cells give melanin, or tanning color, and is where melanoma cancer cells can develop.
- Phenotype Proteins:
  - MiTFg - involved with melanocyte protein and crucial for many cells’ function
  - NGFR - nerve growth factor receptor
  - SOX10 - important in the development of several things, most relevantly formation of melanocytes
  - AXL - pushes cells to divide, can get out of hand in case of cancer, while less of it promotes apoptosis(cells kill themselves)

## Project Objectives
### Our project will be able to answer the following questions:
1. Do phenotype protein levels in experimental condition ‘x’ change over time ‘t’?
  - Method: Two Sample Hypothesis Testing
2. Do phenotype protein levels at a time ‘t’ change between experimental conditions x1 and x2?
  - Method: Two Sample Hypothesis Testing
3. At a time of 0.5h with a dose of 0uM, what are the relationships between different proteins?
  - Method: Confidence Intervals
4. Can we predict cellular phenotype outcomes (Y) values/states from transcription factors?
  - Method: Regression and Correlation
5. Can we determine resulting cellular phenotype depending on levels of 4 phenotype proteins?
  - Method: Regression and Correlation

## Analysis
### Formula: $H_0: p_1 = p_2 \;\; \Longrightarrow \;\; H_0: p_1 - p_2 = 0$
Do proteins in experimental conditions Drug(Vem) and Dose(1uM) change over time(0.5h-120h)?
```{r, echo=FALSE, fig.width=6, fig.height=2.5}
data <- read_csv('data.csv')
data %>% 
  select(NF_kappaB, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 1) -> NF_KappaB_A
data %>%
  select(NF_kappaB, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 7) -> NF_KappaB_B
observed_test_statistic <- mean(NF_KappaB_A$NF_kappaB) - mean(NF_KappaB_B$NF_kappaB)
N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(NF_KappaB_A$NF_kappaB, NF_KappaB_B$NF_kappaB), size=4794+2879, replace=FALSE)
  tmp <- mean(shuffled_xs[1:4794])-mean(shuffled_xs[(4795):(4794+2879)])
  permutation_test_statistics[i] <- tmp
}
plot1 <- tibble("NF_KappaB"=permutation_test_statistics) %>%
  ggplot(aes(x=`NF_KappaB`)) +
  geom_histogram(bins=30, color = "black", fill = "light blue") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(MiTFg, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 1) -> MiTFg_A
data %>%
  select(MiTFg, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 7) -> MiTFg_B
observed_test_statistic <- mean(MiTFg_A$MiTFg) - mean(MiTFg_B$MiTFg)

N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(MiTFg_A$MiTFg, MiTFg_B$MiTFg), size=4794+2879, replace=FALSE)
  tmp <- mean(shuffled_xs[1:4794])-mean(shuffled_xs[(4795):(4794+2879)])
  permutation_test_statistics[i] <- tmp
}
plot2 <- tibble("MiTFg"=permutation_test_statistics) %>%
  ggplot(aes(x=`MiTFg`)) +
  geom_histogram(bins=30, color = "black", fill = "light green") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(AXL, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 1) -> AXL_A
data %>%
  select(AXL, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 7) -> AXL_B
observed_test_statistic <- mean(AXL_A$AXL) - mean(AXL_B$AXL)
N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(AXL_A$AXL, AXL_B$AXL), size=4794+2879, replace=FALSE)
  tmp <- mean(shuffled_xs[1:4794])-mean(shuffled_xs[(4795):(4794+2879)])
  permutation_test_statistics[i] <- tmp
}
plot3 <- tibble("AXL"=permutation_test_statistics) %>%
  ggplot(aes(x=`AXL`)) +
  geom_histogram(bins=30, color = "black", fill = "red") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 1) -> Sox10_A
data %>%
  select(Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 7) -> Sox10_B
observed_test_statistic <- mean(Sox10_A$Sox10) - mean(Sox10_B$Sox10)
N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(Sox10_A$Sox10, Sox10_B$Sox10), size=4794+2879, replace=FALSE)
  tmp <- mean(shuffled_xs[1:4794])-mean(shuffled_xs[(4795):(4794+2879)])
  permutation_test_statistics[i] <- tmp
}
plot4 <- tibble("Sox10"=permutation_test_statistics) %>%
  ggplot(aes(x=`Sox10`)) +
  geom_histogram(bins=30, color = "black", fill = "yellow") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))
gridExtra::grid.arrange(plot1, plot2, plot3, plot4, nrow=2)
```

## Explaination
### Formula: $H_0: p_1 = p_2 \;\; \Longrightarrow \;\; H_0: p_1 - p_2 = 0$
- Do proteins in experimental conditions Drug(Vem) and Dose(1uM) change over time(0.5h-120h)?
  1. We first assume the two populations are the same
  2. Then, use the permutation test in which we shuffle the two samples to determine if there's a difference in means based on the two populations.
  3. Use observed test statistics x1-x2 based on n1 and n2 sample, simulate the Sampling Distribution assuming the NULL Hypothesis is TRUE.
      - Note: the test statistics will be the difference between the statistics of the two groups.
      - Null hypothesis is a statement of equivalence.
      - Re-sample/shuffle both samples and calculate the test statistic each time.
- Conclusion: 
  - The Two Hypothesis Test based on NF_Kappa_B, MiTFg, AXL and Sox10 in experimental conditions Drug(Vem) and Dose(1uM) at timepoint 0.5h and 120h showed p-values of 0, observed statistics of 0.24(NF_KappaB), 0.48(MiTFg), 0.44(AXL), 0.14(Sox10). A p-value of 0 means there is no chance of observing results at least as extreme. We will set an a-significance level of 0.05(also the probability of a Type I error of rejecting a true H0) and we reject H0 at the a-significance level since 0 < 0.05. 

## Analysis
### Formula: Formula: $H_0: p_1 = p_2 \;\; \Longrightarrow \;\; H_0: p_1 - p_2 = 0$
Are protein levels at time 2h different between experimental conditions Vem and Vem+Tram?
```{r, echo=FALSE, fig.width=6, fig.height=2.5}
data <- read_csv('data.csv')
data %>% 
  select(NF_kappaB, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 1, timepoint_id == 2) -> NF_KappaB_A
data %>%
  select(NF_kappaB, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 2, dose_id == 1, timepoint_id == 2) -> NF_KappaB_B
observed_test_statistic <- mean(NF_KappaB_A$NF_kappaB) - mean(NF_KappaB_B$NF_kappaB)

N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(NF_KappaB_A$NF_kappaB, NF_KappaB_B$NF_kappaB), size=6224+6538, replace=FALSE)
  tmp <- mean(shuffled_xs[1:6224])-mean(shuffled_xs[(6225):(6224+6538)])
  permutation_test_statistics[i] <- tmp
}
plot1 <- tibble("NF_KappaB"=permutation_test_statistics) %>%
  ggplot(aes(x=`NF_KappaB`)) +
  geom_histogram(bins=30, color = "black", fill = "light blue") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(MiTFg, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 1, timepoint_id == 2) -> MiTFg_A
data %>%
  select(MiTFg, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 2, dose_id == 1, timepoint_id == 2) -> MiTFg_B
observed_test_statistic <- mean(MiTFg_A$MiTFg) - mean(MiTFg_B$MiTFg)

N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(MiTFg_A$MiTFg, MiTFg_B$MiTFg), size=6224+6538, replace=FALSE)
  tmp <- mean(shuffled_xs[1:6224])-mean(shuffled_xs[(6225):(6224+6538)])
  permutation_test_statistics[i] <- tmp
}
plot2 <- tibble("MiTFg"=permutation_test_statistics) %>%
  ggplot(aes(x=`MiTFg`)) +
  geom_histogram(bins=30, color = "black", fill = "light green") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(AXL, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 1, timepoint_id == 2) -> AXL_A
data %>%
  select(AXL, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 2, dose_id == 1, timepoint_id == 2) -> AXL_B
observed_test_statistic <- mean(AXL_A$AXL) - mean(AXL_B$AXL)

N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(AXL_A$AXL, AXL_B$AXL), size=6224+6538, replace=FALSE)
  tmp <- mean(shuffled_xs[1:6224])-mean(shuffled_xs[(6225):(6224+6538)])
  permutation_test_statistics[i] <- tmp
}
plot3 <- tibble("AXL"=permutation_test_statistics) %>%
  ggplot(aes(x=`AXL`)) +
  geom_histogram(bins=30, color = "black", fill = "red") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))

data %>% 
  select(Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 1, timepoint_id == 2) -> Sox10_A
data %>%
  select(Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 2, dose_id == 1, timepoint_id == 2) -> Sox10_B
observed_test_statistic <- mean(Sox10_A$Sox10) - mean(Sox10_B$Sox10)

N <- 1000
permutation_test_statistics <- 1:N
set.seed(1); for(i in 1:N){
  shuffled_xs <- 
    sample(c(Sox10_A$Sox10, Sox10_B$Sox10), size=6224+6538, replace=FALSE)
  tmp <- mean(shuffled_xs[1:6224])-mean(shuffled_xs[(6225):(6224+6538)])
  permutation_test_statistics[i] <- tmp
}
plot4 <- tibble("Sox10"=permutation_test_statistics) %>%
  ggplot(aes(x=`Sox10`)) +
  geom_histogram(bins=30, color = "black", fill = "yellow") + theme(axis.text=element_text(size=8),axis.title = element_text(size=8,face="bold"))
p4 <- mean(abs(permutation_test_statistics)>=abs(observed_test_statistic))
gridExtra::grid.arrange(plot1, plot2, plot3, plot4, nrow=2)
```

## Explaination
### Formula: Formula: $H_0: p_1 = p_2 \;\; \Longrightarrow \;\; H_0: p_1 - p_2 = 0$
- Are protein levels at time 2h different between experimental conditions Vem and Vem+Tram
  - The method we are using is the same as the previous one.
- Conclusion: 
  - The Two Hypothesis Test is based on NF_Kappa_B, MiTFg, AXL and Sox10 in experimental conditions Drug(Vem) and Dose(0uM) and Drug(Vem + Tram) and Dose(0uM) at timepoint 2h. It showed a p-value of 0 for NF_Kappa_B, AXL and Sox10, and a p-value of 0.74 for MiTFg. The observed statistics are -0.01(NF_Kappa_B), 0.001(MiTFg), -0.02(AXL), -0.1(Sox10). A p-value of 0 means there is no chance of observing results at least as extreme. A p-value of 0.72 means there is 72% chance of observing results at least as extreme. We will set an a-significance level of 0.05(probability of a Type I error of rejecting a true H0) and we reject H0 at the a-significance level since 0 < 0.05, failed to reject H0 at the a-significance level since 0.05 < 0.72. 
  
## Correlation Matrix
### Formula: $H_0: \mu=m_0$ $\;\;|\;\;\;$ $x_i \sim N\left(\mu, \; \sigma\right)$ $\;\;|\;\;\;$ $H_0:
At time 120h in experimental condition Vem and 1uM, what is the relationship between different proteins?
```{r, echo=FALSE, fig.width=12, fig.height=5}
data %>% select(NF_kappaB, MiTFg, AXL, Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 2) %>% select(NF_kappaB, MiTFg, AXL, Sox10) %>% 
  cor() %>% as_tibble(rownames="rowname") %>% pivot_longer(cols=!rowname, names_to="Variable 1",values_to="Correlation") %>%
  rename("Variable 2"=rowname) %>% ggplot(aes(x=`Variable 1`, y=`Variable 2`,
                                              fill=Correlation,
                                              label=round(Correlation,2))) + geom_tile() + geom_text(color="white")
```

## Linear Association
Formula: $$r = \frac{\sum_{i=1}^{n}(x_i-\bar x)(y_i-\bar y) }{\sqrt{\sum_{i=1}^{n}(x_i-\bar x)^2 \sum_{i=1}^{n}(y_i-\bar y)^2}} = \frac{\textcolor{gray}{(n-1)} s_{xy}}{\sqrt{\textcolor{gray}{(n-1)}s_x^2 \textcolor{gray}{(n-1)}s_y^2}} = \frac{s_{xy}}{s_x s_y} = \frac{\text{Cov}(x,y)}{\text{SD}(x) \text{SD}(y)}$$
```{r, echo=FALSE, fig.width=12, fig.height=4}
data %>% 
  select(NF_kappaB, MiTFg, AXL, Sox10, drug_id, dose_id, timepoint_id) %>% 
  filter(drug_id == 1, dose_id == 4, timepoint_id == 7) -> dara
fig1 <- dara %>% 
  ggplot(aes(x=dara$NF_kappaB, y=dara$MiTFg)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$NF_kappaB, y=dara$MiTFg), digits=2), sep=""), x="NF_kappaB", y="MiTFg") + theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
fig2 <- dara %>% 
  ggplot(aes(x=dara$NF_kappaB, y=dara$AXL)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$NF_kappaB, y=dara$AXL), digits=2), sep=""), x="NF_kappaB", y="AXL") + theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
fig3 <- dara %>% 
  ggplot(aes(x=dara$NF_kappaB, y=dara$Sox10)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$NF_kappaB, y=dara$Sox10), digits=2), sep=""), x="NF_kappaB", y="Sox10") + theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
fig4 <- dara %>% 
  ggplot(aes(x=dara$MiTFg, y=dara$AXL)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$MiTFg, y=dara$AXL), digits=2), sep=""), x="MiTFg", y="AXL") +
  theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
fig5 <- dara %>% 
  ggplot(aes(x=dara$MiTFg, y=dara$Sox10)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$MiTFg, y=dara$Sox10), digits=2), sep=""), x="MiTFg", y="Sox10") + theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
fig6 <- dara %>% 
  ggplot(aes(x=dara$AXL, y=dara$Sox10)) + geom_point(alpha=0.5) +
  labs(title=paste("Correlation: ", round(cor(x=dara$AXL, y=dara$Sox10), digits=2), sep=""), x="AXL", y="Sox10") + 
  theme(axis.text=element_text(size=5),axis.title = element_text(size=5))
library(gridExtra)
grid.arrange(fig1, fig2, fig3, fig4, fig5, fig6,ncol=3, nrow=2)
```

## Corelation Explained:
### At time 120h in experimental condition Vem and 1uM, what is the relationship between different proteins?
- Our ideas
  - The denominator scales the numerator so that the total $-1 \leq r \leq 1$ always
  - $r$ measures *linear association*, with $r>0$ *positive* and $r<0$ means *negative*
- As we can see from the graph:
  1. NF_kappa_B is somewhat associated with MiTFg. (Positive i.e. Increasing Correlation of 0.57)
  2. NF_kappa_B is strongly associated with AXL. (Positive i.e. Increasing Correlation of 0.71)
  3. NF_kappa_B is weakly associated with Sox10. (Positive i.e. Increasing Correlation of 0.14)
  4. MiTFg is somewhat associated with Sox10. (Positive i.e. Increasing Correlation of 0.54)
  5. MiTFg is very weakly associated with AXL. (Negative i.e. Decreasing Correlation of -0.01)
  6. AXL is very weakly associated with Sox10. (Negative i.e. Decreasing Correlation of -0.07)

## Bootstraping and Confidence Intervals
### Formula: $\bar{x} \pm z \frac{s}{\sqrt{n}}$
```{r, echo=FALSE, fig.width=12, fig.height=5}
test_stats2 <- data %>% select(MiTFg, Sox10, NF_kappaB, AXL, Timepoint, drug_id, dose_id) %>% filter(drug_id == 1, dose_id == 1, Timepoint == '0.5 h')
vector_test_stats2 <- pull(test_stats2, MiTFg)
vector_test_stats3 <- pull(test_stats2, Sox10)
set.seed(1)
num_trials <- 100
sample_size <- 1000
num_samples <- 1000

mu <- mean(vector_test_stats2) - mean(vector_test_stats3);
half_alpha = 0.1
percentile <- c(half_alpha, 1 - half_alpha)
plot <- ggplot() 

for (i in 1:num_trials){  
  bootstrap_trial <- sample(vector_test_stats2, size = sample_size, replace = TRUE)  #sample to sample from
  bootstrap_trial2 <- sample(vector_test_stats3, size = sample_size, replace = TRUE)
  bootstrap_samples <- 1:num_samples
  for (j in 1:num_samples){  #take a bunch of bootstrap samples
    tmp <- (sample(bootstrap_trial, size = sample_size, replace = TRUE))
    tmp2 <- (sample(bootstrap_trial2, size = sample_size, replace = TRUE))
    bootstrap_samples[j] <- mean(tmp) - mean(tmp2)
  }
  #calculate confidence interval
  ConfidenceInterval <- quantile(bootstrap_samples, percentile)
  if(all(ConfidenceInterval < mu) | all(ConfidenceInterval > mu)){
  col="red"}else{col="black"}
  plot <- plot + geom_line(color = col, data = tibble(x=ConfidenceInterval, y=c(i,i)), aes(x=x, y=y))
}

a <- plot+xlab("Mean Difference Between MiTFg and Sox10 Levels")+labs(title= paste((1-2*half_alpha)*100, "% Confidence Intervals", sep="")) + theme(axis.title.y=element_blank())


test_stats2 <- data %>% select(MiTFg, Sox10, NF_kappaB, AXL, Timepoint, drug_id, dose_id) %>% filter(drug_id == 1, dose_id == 1, Timepoint == '0.5 h')
vector_test_stats2 <- pull(test_stats2, MiTFg)
vector_test_stats3 <- pull(test_stats2, Sox10)
set.seed(1)
num_trials <- 100
sample_size <- 1000
num_samples <- 1000

mu <- mean(vector_test_stats2) - mean(vector_test_stats3);
half_alpha = 0.025
percentile <- c(half_alpha, 1 - half_alpha)
plot <- ggplot() 

for (i in 1:num_trials){  
  bootstrap_trial <- sample(vector_test_stats2, size = sample_size, replace = TRUE)  #sample to sample from
  bootstrap_trial2 <- sample(vector_test_stats3, size = sample_size, replace = TRUE)
  bootstrap_samples <- 1:num_samples
  for (j in 1:num_samples){  #take a bunch of bootstrap samples
    tmp <- (sample(bootstrap_trial, size = sample_size, replace = TRUE))
    tmp2 <- (sample(bootstrap_trial2, size = sample_size, replace = TRUE))
    bootstrap_samples[j] <- mean(tmp) - mean(tmp2)
  }
  #calculate confidence interval
  ConfidenceInterval1 <- quantile(bootstrap_samples, percentile)
  if(all(ConfidenceInterval1 < mu) | all(ConfidenceInterval1 > mu)){
  col="red"}else{col="black"}
  plot <- plot + geom_line(color = col, data = tibble(x=ConfidenceInterval1, y=c(i,i)), aes(x=x, y=y))
}

b <- plot+xlab("Mean Difference Between MiTFg and Sox10 Levels")+labs(title= paste((1-2*half_alpha)*100, "% Confidence Intervals", sep="")) + theme(axis.title.y=element_blank())
gridExtra::grid.arrange(a,b, ncol=2)
```


```{r, echo=FALSE, message=FALSE, warning=FALSE}
#Data management (split based on drug id, dose id, timepoint, and Rep)
# Separate Rep
drug1_dose5_t4_r1 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 1)
drug1_dose5_t4_r2 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 2)
drug1_dose5_t4_r3 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 3)

# All Rep
drug1_dose5_t4 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5)
drug1_dose5_t4 <- drug1_dose5_t4 %>% relocate(NGFR, .after = NF_kappaB)

```

## Confidence Intervals Explained
### Relationships Between Different Proteins at time of 0.5h with dose 0uM
- With confidence intervals, we can estimate many parameter values between our phenotype proteins. For example, if we were to keep all other factors constant(so in this case we are choosing a time point of 0.5h, drug dose of 0uM), we could estimate the difference between means of phenotype proteins MiTFg and Sox10.\
- We performed bootstrapping by taking 2 samples, one for each phenotype protein level, from     the population with our condition/time restrictions.\ 
- Then we resample(with replacement) from each sample, and calculate the difference between the  means of the 2 phenotype protein levels.\
- As our example we used the last bar for our 95% confidence interval and we can say we are 95%  confident that the interval (0.06317026, 0.08802395) will contain the true difference           between the means. This is a pretty tight interval with high confidence, so this is an         accurate and appropriate estimation for the true parameter.\
- I’ve simulated this bootstrapping 100 times to demonstrate the idea of 95% confidence(which are signified by the black intervals, making up the majority) and we could theoretically repeat this process between different proteins, conditions, and confidence intervals to get a better understanding of the relationships between protein levels

## Predicting Cellular Phenotype Outcomes From Transcription Factors
- Here we got R to calculate the pearson correlation coefficient between all 26 proteins(4 phenotype and 22 transcriptions) with conditions of Vem drug, 15 hours after administration, and dose of 3.16 uM.\ 
- Keeping in concise, the correlation between 2 proteins is calculated by Cov(x,y) ÷ SD(x)SD(y), where we divide the covariance(squared difference in values) by the standard deviation in order to standardize the value \
- We also decided to classify the phenotype values into “High” and “Low” categories (continuous data  → discrete data) to fit a classification model on it \
- The graphs below depict the distribution of the value of phenotype proteins in the cell (which follows normal distribution) \
- The distribution for phenotypes across repetitions is generally similar, therefore we can use the data across repetitions to classify any values for the phenotypes above the median as “high” and below the median as “low” 

## An Overview of Means of the indicator proteins
```{r, echo=FALSE, fig.width=12, fig.height=5}
#Data management (split based on drug id, dose id, timepoint, and Rep)
# Separate Rep
drug1_dose5_t4_r1 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 1)
drug1_dose5_t4_r2 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 2)
drug1_dose5_t4_r3 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5, Rep == 3)

# All Rep
drug1_dose5_t4 <- data %>% filter(drug_id == 1, timepoint_id == 4, dose_id == 5)
drug1_dose5_t4 <- drug1_dose5_t4 %>% relocate(NGFR, .after = NF_kappaB)

```

```{r, echo=FALSE, fig.width=12, fig.height=5}
h1 <- drug1_dose5_t4 %>% ggplot(aes(x = NGFR)) + geom_boxplot(fill="light blue") + xlim(2.5, 4.5)
h2 <- drug1_dose5_t4 %>% ggplot(aes(x = MiTFg)) + geom_boxplot(fill="light green") + xlim(2.5, 4.5)
h3 <- drug1_dose5_t4 %>% ggplot(aes(x = AXL)) + geom_boxplot(fill="red") + xlim(2.5, 4.5)
h4 <- drug1_dose5_t4 %>% ggplot(aes(x = Sox10)) + geom_boxplot(fill="yellow") + xlim(2.5, 4.5)
medians <- drug1_dose5_t4 %>% select(c(NGFR, MiTFg, AXL, Sox10)) %>% summarise(m_NGFR = median(NGFR), m_MiTFg = median(MiTFg), m_AXL = median(AXL), m_Sox10 = median(Sox10))
m_NGFR <- medians$m_NGFR
m_MiTFg <- medians$m_MiTFg
m_AXL <- medians$m_AXL
m_Sox10 <- medians$m_Sox10
grid.arrange(h1,h2,h3,h4)
```

```{r, echo=FALSE, message=FALSE, warning=FALSE}
drug1_dose5_t4 <- drug1_dose5_t4 %>% mutate("NGFR_level" = case_when(NGFR > m_NGFR ~ "HIGH", TRUE ~ "LOW"),
                                                  "MiTFg_level" = case_when(MiTFg > m_MiTFg ~ "HIGH", TRUE ~ "LOW"),
                                                  "AXL_level" = case_when(AXL > m_AXL ~ "HIGH", TRUE ~ "LOW"),
                                                  "Sox10_level" = case_when(Sox10 > m_Sox10 ~ "HIGH", TRUE ~ "LOW"))
```

## Continuing with Multivariate Linear Regression and Classification
- Utilizing the same conditions as above, we can create multivariate linear regression models using all available proteins as our explanatory variables. In real life, it may be costly to have to record and calculate data for all 25 protein levels just to predict 1 protein, but we are already provided with full data on all 26 protein levels. Since we are only predicting the phenotype values we don’t need to worry about multicollinearity. \
- Another concern that may arise with this many explanatory variables is overfitting, but we perform an 80-20 train split test(save 80% of data to build our models on, and test them out on the new 20% of the rest of data), then calculate the root mean square error(RMSE) to estimate the spread of our error. And the results show us the RMSE is quite low for the linear models of all 4 phenotype proteins. NGFR has the best model at about 0.0795 RMSE and AXL was the predictive model at about 0.137 RMSE \

```{r, echo=FALSE, message=FALSE, warning=FALSE}
# 80-20 split
n <- nrow(drug1_dose5_t4)
n_train <- as.integer(n*0.8)
n_test <- n - n_train
training_indices <- 
  sample(1:n,size=n_train,replace=FALSE)
 drug1_dose5_t4 <- drug1_dose5_t4 %>% rowid_to_column()
train <- drug1_dose5_t4 %>% 
  filter(rowid %in% training_indices)
test <- drug1_dose5_t4 %>% 
  filter(!(rowid %in% training_indices))

NGFR_test <- test$NGFR
MiTFg_test <- test$MiTFg
AXL_test <- test$AXL
Sox10_test <- test$Sox10
```

```{r, echo=FALSE, message = FALSE, warning=FALSE}
# NGFR
model_NGFR = lm(NGFR ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)


# MiTFg
model_MiTFg = lm(MiTFg ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

# AXL
model_AXL = lm(AXL ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

# Sox10
model_Sox10 = lm(Sox10 ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

```

```{r, echo=FALSE, message = FALSE, warning = FALSE}
model_NGFR_test <- predict(model_NGFR, newdata=test)
model_MiTFg_test <- predict(model_MiTFg, newdata=test)
model_AXL_test <- predict(model_AXL, newdata=test)
model_Sox10_test <- predict(model_Sox10, newdata=test)

model_NGFR_test_RMSE <- sqrt(mean((NGFR_test-model_NGFR_test)^2))
model_MiTFg_test_RMSE <- sqrt(mean((MiTFg_test-model_MiTFg_test)^2))
model_AXL_test_RMSE <- sqrt(mean((AXL_test-model_AXL_test)^2))
model_Sox10_test_RMSE <- sqrt(mean((Sox10_test-model_Sox10_test)^2))

library(glue)
glue("NGFR RMSE: ", round(model_NGFR_test_RMSE, digits=2), " | MiTFg RMSE: ", round(model_MiTFg_test_RMSE, digits=2), " | AXL RMSE: ", round(model_AXL_test_RMSE, digits=2), " | Sox10 RMSE: ", round(model_Sox10_test_RMSE, digits=2))
```

## Classification
- Below we built a classification tree based on all 25 explanatory variables as the covariates for MiTFG levels, we have decided to use all 25 seeing as we have access to all 25, however, this may have diminishing returns. It is not too complicated due to the stopping rules that don’t allow trees to split if there is not a significant increase in classification. Similarly, we can repeat this for all the other proteins. \
- The calculated accuracy(proportion we got correct), precision(proportion of ones we identified as high that are correct), sensitivity(proportion of actually high levels we identified as high), and specificity(proportion of actually low we identified as low) are all in the 0.61 to 0.77 range, indicating our tree is a solid predictor of these values \
- Finally, we used to trained classification model to predict the levels of the phenotypes (High or Low) and used those values to classify the types of cells (Undifferentiated, Transitory, etc)\

## Decision Tree for MiTFg
```{r, echo=FALSE, fig.width=12, fig.height=5}
library(rpart)
library(partykit)
set.seed(1)
tree <- rpart(MiTFg_level ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)
tree %>% as.party() %>%
plot(type="simple",gp=gpar(cex=0.7))

tree_test <- predict(tree, type = "class", newdata = test)
tree_test_confusion_matrix <- table(`y-hat`=tree_test,'observed y'=test$MiTFg_level)
n <- sum(tree_test_confusion_matrix)
n_TN <- tree_test_confusion_matrix[1,1]
n_FN <- tree_test_confusion_matrix[1,2]
n_FP <- tree_test_confusion_matrix[2,1]
n_TP <- tree_test_confusion_matrix[2,2]

glue("Accuracy: ", round((n_TP+n_TN)/n, digits=2), " | Precision: ", round(n_TP / (n_TP + n_FP), digits=2), " | Sensitivity: ", round(n_TP / (n_TP + n_FN), digits=2), " | Specificity: ", round(n_TN / (n_TN + n_FP), digits=2))
```
```{r, echo=FALSE, message = FALSE, warning = FALSE}
# NGFR
tree_NGFR <- rpart(NGFR_level ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

# MiTFg
tree_MiTFg <- rpart(MiTFg_level ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

# AXL
tree_AXL <- rpart(AXL_level ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)

# Sox10
tree_Sox10 <- rpart(Sox10_level ~ Phospho_c_Fos + Phospho_c_Jun + Phospho_ATF2 + Phospho_Fra1 + c_Fos + c_Jun + Fra1 + JunD + ATF2 + JunB + Fra2 + ATF4 + Phospho_ATF4 + Phospho_Erk1 + Phospho_ATF1 + ATF6 + Phospho_S6 + ATF3 + ATF5 + Phospho_p38 + Ki_67 + NF_kappaB, data = train)
```

## Conclusions Visualized
- The answer to Question 4 is yes. We can see that the low RMSE values for each phenotype and the fairly high accuracy of the classification model implies that our model fits the data well. For example, classifying a cell as melanocytic may suggest cells in early stages of cancer.
```{r, echo=FALSE, fig.width=12, fig.height=3}
tree_NGFR_pred <- predict(tree_NGFR, type = 'class', newdata = test)
tree_MiTFg_pred <- predict(tree_MiTFg, type = 'class', newdata = test)
tree_AXL_pred <- predict(tree_AXL, type = 'class', newdata = test)
tree_Sox10_pred <- predict(tree_Sox10, type = 'class', newdata = test)

predicted_ph_levels <- tibble(NGFR_level_hat = tree_NGFR_pred, MiTFg_level_hat = tree_MiTFg_pred, AXL_level_hat = tree_AXL_pred, Sox10_level_hat = tree_Sox10_pred)
predicted_ph_levels <- predicted_ph_levels %>% mutate("Type of cell" = case_when(NGFR_level_hat == "LOW" &
                                                                                   MiTFg_level_hat == "LOW" &
                                                                                   AXL_level_hat == "HIGH" &
                                                                                   Sox10_level_hat == "LOW" ~ "Undifferentiated",
                                                                                 NGFR_level_hat == "HIGH" &
                                                                                   MiTFg_level_hat == "LOW" &
                                                                                   AXL_level_hat == "HIGH" &
                                                                                   Sox10_level_hat == "HIGH" ~ "Neural crest-like",
                                                                                 NGFR_level_hat == "HIGH" &
                                                                                   MiTFg_level_hat == "HIGH" &
                                                                                   AXL_level_hat == "LOW" &
                                                                                   Sox10_level_hat == "HIGH" ~ "Transitory",
                                                                                 NGFR_level_hat == "HIGH" &
                                                                                   MiTFg_level_hat == "HIGH" &
                                                                                   AXL_level_hat == "LOW" &
                                                                                   Sox10_level_hat == "LOW" ~ "Melanocytic"))
predicted_ph_levels
```

## Project Acknowledgements
- Our ideas come from Dr. Scott Schwartz who is a seasoned professional at the University of Toronto worked in Integrative Biology, Nutrition and Complex Disease, and Next Generation Sequencing labs.
Link: https://github.com/pointOfive

- Our data come from the article (Refer to the below link) finding that the "AP-1 transcription factor network" (i.e., the relative distributions and dependency relationships of transcription factors) are predictive of "cellular plasticity in melanoma" (i.e., how easily changeable the phenotype are melanoma cell lines)
Link: https://www.biorxiv.org/content/10.1101/2021.12.06.471514v1.full

## Project References
- NCI Dictionaries, National Cancer Institute
  - https://www.cancer.gov/publications/dictionaries
- The Neural Crest and Cancer: A Developmental Spin on Melanoma, National Library of Medicine
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3809092
- What Is Melanoma Skin Cancer?
  - https://www.cancer.org/cancer/melanoma-skin-cancer/about/what-is-melanoma.html
- AP-1 transcription factor network explains diverse patterns of cellular plasticity in melanoma
  - https://www.biorxiv.org/content/10.1101/2021.12.06.471514v1.full
- MITF gene
  - https://medlineplus.gov/genetics/gene/mitf
- SOX10 gene
  - https://medlineplus.gov/genetics/gene/sox10/
- Nerve Growth Factor Receptor
  - https://www.sciencedirect.com/topics/medicine-and-dentistry/nerve-growth-factor-receptor
- AXL receptor tyrosine kinase as a promising anti-cancer approach
  - https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-019-1090-3

# ALTRA Microbiome Analysis

<br>

## Introduction

## The Data

``` r
# Import data
data <- paste(getwd(), "/../data/ALTRA_clinical&16S.merged.19-Dec-2022.analysis_samples.csv", sep="") %>% 
  read_delim(delim = ",") %>%
  mutate(rowname = ifelse(sample_type_16S == "Stool", paste(rowname, "F", sep = ""), rowname)) %>%
  column_to_rownames("rowname")

# OTU table
otu_table <- data %>% select(contains("Bacteria")) %>%
  t() %>% otu_table(taxa_are_rows = T)

# Sample metadata
metadata <- data %>% select(-contains("Bacteria")) %>%
  select(sample_id, ccp3_group, ccp3, everything())

# Taxonomy matrix
otu_names <- otu_table %>% as.data.frame() %>% rownames()
tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_tax_split <- sapply(otu_names, FUN=function(x) strsplit(x, "/"))
tax_table <- plyr::ldply(otu_tax_split, rbind)[-1] %>%
  set_rownames(otu_names) %>%
  set_colnames(tax_levels) %>%
  as.matrix()

# Phyloseq object
physeq <- phyloseq(otu_table(otu_table), tax_table(tax_table), sample_data(metadata))
```

### Missing Data

``` r
# Visualize missingness
data %>% select(-contains("Bacteria/")) %>% vis_miss() +
  scale_fill_manual(values = c("gray15", "darkgoldenrod")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0)) +
  guides(fill = guide_legend(title = "Missing"))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-3-1.png)<!-- --> <br>
<br>

### Subject Counts

``` r
# Subject counts by grouping
# All samples
subject_counts <- data %>% select(sample_id, ccp3_group) %>% 
  unique() %>% select(ccp3_group) %>% table() %>% as.data.frame() %>%
  column_to_rownames(".") %>% dplyr::rename("All Subjects" = "Freq") %>% t()

# Stool samples
stool_counts <- data %>% filter(sample_type_16S == "Stool") %>% select(sample_id, ccp3_group) %>% 
  select(ccp3_group) %>% table() %>% as.data.frame() %>%
  column_to_rownames(".") %>% dplyr::rename("Stool" = "Freq") %>% t()

# Sputum samples
sputum_counts <- data %>% filter(sample_type_16S == "Sputum") %>% select(sample_id, ccp3_group) %>% 
  select(ccp3_group) %>% table() %>% as.data.frame() %>%
  column_to_rownames(".") %>% dplyr::rename("Sputum" = "Freq") %>% t()

# Output table
rbind(subject_counts, stool_counts, sputum_counts) %>% knitr::kable()
```

|              | NegControl | PosConverter | PosNonconverter | PosRA |
|:-------------|-----------:|-------------:|----------------:|------:|
| All Subjects |         32 |           14 |              40 |     5 |
| Stool        |         32 |           14 |              40 |     5 |
| Sputum       |         19 |            8 |              38 |     2 |

<br> <br>

## EDA

### Alpha Diversity by CCP-Group

The violin/boxplots below show the distributions of various alpha
diversity measurements, stratified by CCP-group and sample type.  
<br>

#### Shannon Diversity

``` r
plot_alpha_div(metadata, y="ShannonH.Median", title="Shannon Diversity")
```

![](altra-inference_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### Species Richness

``` r
plot_alpha_div(metadata, y="Sobs.Median", title="Species Richness")
```

![](altra-inference_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### Species Evenness

``` r
plot_alpha_div(metadata, y="ShannonE.Median", title="Species Evenness")
```

![](altra-inference_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Community Compositions

Below are stacked bar charts showing community compositions. Bars are
ordered by community similarity, as determined through hierarchical
clustering (using Bray-Curtis distance). Squares at the top of each bar
represent CCP group (see legend below). The top 15 taxa observed among
all samples are included - rarer taxa are grouped into the “Other”
group. Taxa that were did meet detection criteria (1e-4% relative
abundance in at least 10% of samples) are not included, which is why not
all stacked bar charts reach 100%.  
<br>

<img src="/Users/bgraham/Github/altra/inference/ccp_label.png" id="id"
class="class" style="width:20.0%;height:5.0%" />

#### Stool, Species-Level

``` r
sample_type = "Stool"
ccp3_groups = c("NegControl", "PosNonconverter", "PosConverter")

# Select 
ra.physeq <- physeq %>%
    microbiome::transform("compositional") %>%         # Transform to RA
    subset_samples(sample_type_16S == sample_type) %>% # Pick sample type
    subset_samples(ccp3_group %in% ccp3_groups) %>%    # Choose groups
    core(detection=1e-6, prevalence=0.10)              # Pick the core
```

``` r
# Display compositional barchart
composition_barchart(ra.physeq=ra.physeq, marker_size=2)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# Set cluster and pheatmap data
cluster.df <- ra.physeq %>% 
  otu_table(taxa_are_rows = TRUE) %>% t() %>% as.data.frame() %>% 
  set_colnames(abbrev_taxa(ra.physeq))

# Create pheatmap
core.physeq <- ra.physeq %>% core(detection=1/100, prevalence=25/100)
pheatmap.df <- core.physeq %>%
  otu_table() %>% t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(core.physeq))

# Pheatmap annotations
annotation <- ra.physeq %>% sample_data() %>% as_tibble() %>% as.data.frame() %>%
  select(CCP = ccp3, contains("95_pos")) %>%
  mutate(
    CCP = ifelse(CCP == 0, "-", "+"),
    sp_rf_ig_m_95_pos = ifelse(sp_rf_ig_m_95_pos == 0, "-", "+"),
    sp_rf_ig_a_95_pos = ifelse(sp_rf_ig_a_95_pos == 0, "-", "+"),
    sp_ccp_ig_a_95_pos = ifelse(sp_ccp_ig_a_95_pos == 0, "-", "+"),
    sp_ccp_ig_g_95_pos = ifelse(sp_ccp_ig_g_95_pos == 0, "-", "+")) %>%
  set_rownames(colnames(t(pheatmap.df)))

# Set colors
annotation_color = list(ccp3 = c(CCPminus = "gray10", CCPplus = "darkgoldenrod"))
```

``` r
# Display pheatmap
pheatmap(t(pheatmap.df), legend = F, color = viridis(100),
         cluster_cols = cluster.df %>% vegdist(method = "bray") %>% hclust(method = "ward.D2"),
         annotation = annotation, annotation_colors = annotation_color,
         border_color = NA, cluster_rows = F, show_colnames = F)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Sputum, Species-Level

``` r
sample_type = "Sputum"
ccp3_groups = c("NegControl", "PosNonconverter", "PosConverter", "PosRA")

ra.physeq <- physeq %>%
    microbiome::transform("compositional") %>%         # Transform to RA
    subset_samples(sample_type_16S == sample_type) %>% # Pick sample type
    subset_samples(ccp3_group %in% ccp3_groups) %>%    # Choose groups
    core(detection=1e-6, prevalence=0.10)           # Pick the core
```

``` r
# Display compositional barchart
composition_barchart(ra.physeq=ra.physeq, marker_size=2.8)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Set cluster and pheatmap data
cluster.df <- ra.physeq %>% 
  otu_table(taxa_are_rows = TRUE) %>% t() %>% as.data.frame() %>% 
  set_colnames(abbrev_taxa(ra.physeq))

# Create pheatmap
core.physeq <- ra.physeq %>% core(detection=1/100, prevalence=25/100)
pheatmap.df <- core.physeq %>%
  otu_table() %>% t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(core.physeq))

# # Pheatmap annotations
annotation <- ra.physeq %>% sample_data() %>% as_tibble() %>% as.data.frame() %>%
  select(CCP = ccp3, contains("95_pos")) %>%
  mutate(
    CCP = ifelse(CCP == 0, "-", "+"),
    sp_rf_ig_m_95_pos = ifelse(sp_rf_ig_m_95_pos == 0, "-", "+"),
    sp_rf_ig_a_95_pos = ifelse(sp_rf_ig_a_95_pos == 0, "-", "+"),
    sp_ccp_ig_a_95_pos = ifelse(sp_ccp_ig_a_95_pos == 0, "-", "+"),
    sp_ccp_ig_g_95_pos = ifelse(sp_ccp_ig_g_95_pos == 0, "-", "+")) %>%
  set_rownames(colnames(t(pheatmap.df)))

# Set colors
annotation_color = list(ccp3 = c(CCPminus = "gray10", CCPplus = "darkgoldenrod"))
```

``` r
# Display pheatmap
pheatmap(
  t(pheatmap.df), legend = F, color = viridis(100),
  cluster_cols = cluster.df %>% vegdist(method = "bray") %>% hclust(method = "ward.D2"),
  annotation = annotation, annotation_colors = annotation_color,
  border_color = NA, cluster_rows = F, show_colnames = F)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

<br> <br>

## Inference

### Hypothesis 1

#### CCP(+) vs CCP(-)

*There are microbiome differences between groups suggesting
relationships between ‘autoimmune states’.*  
<br> For this section, we will consider CCP(+) individuals as those
individuals who are CCP(+) but do not have RA (both CCPPosNonconverters
and CCPPosConverters). We will compare this group to the CCP(-)
(NegControl) group.  
<br>

#### Stool

``` r
# Pick relative abundances (compositional) and sample metadata
sample_type = "Stool"
pseq <- physeq %>% 
  subset_samples(sample_type_16S == sample_type) %>%
  subset_samples(ccp3_group != "PosRA") #%>%
  #tax_glom(taxrank="Genus") %>%
  
pseq.rel <- pseq %>% 
  microbiome::transform("compositional") #%>%
  #core(detection = 0.01, prevalence = 0.50)
otu <- abundances(pseq.rel) %>% t() %>% as.data.frame()
meta <- meta(pseq.rel) %>% mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+")))
```

##### Alpha Diversity

``` r
# Shannon diversity plot
shannon <- meta %>%
  ggplot(aes(x=ccp3, y=ShannonH.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Shannon Diversity") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Sobs diversity plot
sobs <- meta %>%
  ggplot(aes(x=ccp3, y=Sobs.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Richness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Species Evenness plot
evenness <- meta %>%
  ggplot(aes(x=ccp3, y=ShannonE.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Evenness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Display three on one plot
ggarrange(shannon, sobs, evenness, ncol = 3, nrow = 1)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# Perform tests
shan.stat.test <- wilcox.exact(ShannonH.Median ~ ccp3, data=meta, paired=F)
rich.stat.test <- wilcox.exact(Sobs.Median ~ ccp3, data=meta, paired=F)
even.stat.test <- wilcox.exact(ShannonE.Median ~ ccp3, data=meta, paired=F)

# Output p-value table
data.frame(
  "Measurement" = c("Shannon Diversity", "Species Richness", "Species Evenness"),
  "p.val" = c(
    paste("p = ", round(shan.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(rich.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(even.stat.test$p.value, 2), sep = ""))) %>% 
  knitr::kable()
```

| Measurement       | p.val    |
|:------------------|:---------|
| Shannon Diversity | p = 0.98 |
| Species Richness  | p = 0.15 |
| Species Evenness  | p = 0.66 |

##### PCoA

``` r
# Need core taxa to save time
pcoa_otu <- pseq %>% 
  microbiome::transform("compositional") %>%
  core(detection = 0.01, prevalence = 0.10)
core_taxa <- abbrev_taxa(pcoa_otu) # Extract shortened taxa names
pcoa_otu <- pcoa_otu %>%
  abundances() %>% 
  t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(pcoa_otu))

# Determine coordinates for samples
PCoA <- vegdist(pcoa_otu, method="bray") %>%
  cmdscale() %>%
  as.data.frame() %>%
  select(Dim1=`V1`, Dim2=`V2`)

# Get vectors for taxa
taxa_vectors <- envfit(ord = PCoA, env = pcoa_otu)
taxa_vector_coords <- taxa_vectors$vectors$arrows * sqrt(taxa_vectors$vectors$r)
taxa_vector_p.vals <- taxa_vectors$vectors$pvals
vector_df <- data.frame(p_val = taxa_vector_p.vals) %>%
  bind_cols(taxa_vector_coords) %>%
  rownames_to_column("Taxa") %>%
  filter(p_val <= 0.05) %>%
  arrange(p_val) %>% head(5)

# Add metadata to ordination values
pcoa_plot_df <- PCoA %>% 
  merge(pseq %>% sample_data() %>% as.data.frame(), by = 'row.names') %>%
  mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+"))) %>%
  column_to_rownames('Row.names')

# Ordination bi-plot
pcoa_plot_df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = ccp3)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_segment(data = vector_df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               #arrow = arrow(length = unit(0.2, "cm")), 
               colour = "black", stat = "identity", alpha = 0.7, inherit.aes = FALSE) +
  geom_text_repel(data = vector_df, #vjust = "inward", hjust = "inward",
                  aes(x = Dim1, y = Dim2, label = Taxa), 
                  inherit.aes = FALSE, size=3) +
  theme_bw() +
  scale_color_manual(values = c("gray10", "darkgoldenrod")) +
  ggtitle("PCoA, Beta Diversity", "Bray-Curtis") +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 12))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# Ordination bi-plot - color differently
pcoa_plot_df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = ccp3_group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_segment(data = vector_df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               #arrow = arrow(length = unit(0.2, "cm")), 
               colour = "black", stat = "identity", alpha = 0.7, inherit.aes = FALSE) +
  geom_text_repel(data = vector_df, #vjust = "inward", hjust = "inward",
                  aes(x = Dim1, y = Dim2, label = Taxa), 
                  inherit.aes = FALSE, size=3) +
  theme_bw() +
  scale_color_manual(values = c("gray10", "darkgoldenrod1", "darkred")) +
  ggtitle("PCoA, Beta Diversity", "Bray-Curtis") +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 12))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

##### PERMANOVA

``` r
# PERMANOVA test using Bray-Curtis distance
set.seed(007) # Set seed for reproducibility - permutation-based test
permanova <- adonis2(
  otu ~ ccp3 + age + gender + race, data = meta, by="margin", permutations = 999, method = "bray")

permanova %>% as.data.frame() %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  dplyr::rename("p.val" = "Pr(>F)") %>%
  mutate(p.val = ifelse(p.val <= 0.05, paste("**", p.val, "****", sep=""), p.val)) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>%
  knitr::kable(align = 'ccccc')
```

|          | Df  | SumOfSqs |  R2  |  F   |     p.val     |
|:---------|:---:|:--------:|:----:|:----:|:-------------:|
| ccp3     |  1  |   0.37   | 0.06 | 4.88 | **0.001**\*\* |
| age      |  1  |   0.08   | 0.01 | 1.04 |     0.375     |
| gender   |  1  |   0.12   | 0.02 | 1.55 |     0.104     |
| race     |  6  |   0.34   | 0.05 | 0.75 |     0.833     |
| Residual | 76  |   5.74   | 0.86 |  NA  |      NA       |
| Total    | 85  |   6.70   | 1.00 |  NA  |      NA       |

**Dispersions Plot**

``` r
# Plot dispersion distances for each "group"
beta_dispersion <- otu %>% vegdist(method = "bray") %>% betadisper(meta$ccp3)
plot(beta_dispersion, hull=FALSE, ellipse=TRUE)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

**Homogeneity of Dispersons**

``` r
# Hypothesis test
set.seed(007)
otu %>% vegdist() %>% betadisper(meta$ccp3) %>% permutest()
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00083 0.0008345 0.1086    999  0.764
    ## Residuals 84 0.64537 0.0076830

##### Differential Abundance Analysis

The model reports “ccp3_CCPplus_vs_CCPminus”, which should indicate that
the **CCPminus group is the reference**. An $\alpha$-level of 0.05 is
used as the threshold for selecting taxa.

``` r
# Convert physeq object to deseq and fit model
deseq <- phyloseq_to_deseq2(pseq, ~ ccp3) # Convert physeq object to deseq
fit <- DESeq2::DESeq(deseq, test="Wald", fitType="parametric") # Fit model

# Taxonomy info to get short names
tax_info <- pseq %>% tax_table() %>% as.data.frame() %>%
  mutate(short_name = abbrev_taxa(pseq))

# Extract, filter and sort results
options(digits = 3)
results <- DESeq2::results(fit, cooksCutoff = F, tidy = TRUE) %>%
  filter(padj < 0.05) %>% # Select significant p-vals
  column_to_rownames("row") %>%
  merge(tax_info, by = "row.names") %>%
  arrange(log2FoldChange) # Sort by log2FoldChange

# Plot
order <- results$short_name
results$short_name <- factor(results$short_name, levels = order)

# Plot taxa that fit core criteria and are also differentially abundant
core_taxa <- pseq.rel %>% core_members(detection = 1/100, prevalence = 10/100) # Filter out rare stuff
significant_taxa <- results$Row.names
filtered_sig_taxa <- intersect(core_taxa, significant_taxa) # Taxa both significant and not rare
sig_taxa.df <- pseq.rel %>% 
  otu_table() %>% t() %>% as.data.frame() %>%
  select(filtered_sig_taxa) %>%
  merge(meta %>% select(ccp3), by = "row.names") %>%
  column_to_rownames("Row.names")

sig_taxa.df %>%
  melt(idvars = ccp3) %>%
  dplyr::rename(Row.names = variable, rel_ab = value) %>%
  merge(results %>% select(Row.names, short_name), by = "Row.names") %>%
  group_by(ccp3, short_name) %>%
  dplyr::summarize(Median.RA = median(rel_ab) * 100, IQR = IQR(rel_ab) * 100) %>%
  ungroup() %>%
  dplyr::rename(Taxa = short_name) %>%
  arrange(Taxa, ccp3) %>% knitr::kable()
```

| ccp3 | Taxa                  | Median.RA |    IQR |
|:-----|:----------------------|----------:|-------:|
| \-   | Peptostreptococcaceae |     3.889 |  5.686 |
| \+   | Peptostreptococcaceae |     1.016 |  2.908 |
| \-   | Coriobacteriaceae     |     0.474 |  0.784 |
| \+   | Coriobacteriaceae     |     0.469 |  0.475 |
| \-   | Blautia               |    17.429 | 15.685 |
| \+   | Blautia               |    18.102 | 10.125 |
| \-   | Bacteroidales         |     0.028 |  0.115 |
| \+   | Bacteroidales         |     0.140 |  0.830 |
| \-   | Bacteroides           |     0.245 |  0.843 |
| \+   | Bacteroides           |     2.797 |  7.502 |

``` r
results %>%
  filter(Row.names %in% filtered_sig_taxa) %>%
  ggplot(aes(x = short_name, y = log2FoldChange, fill = Phylum)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c("#EBCC2A", "#F21A00", "#3B9AB2")) +
  ggtitle("Changes in Relative Abundance for Significant Taxa", "CCP(+) vs CCP(-)") + 
  theme_bw() +
  theme(axis.title.y = element_blank())
```

![](altra-inference_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
# maaslin2.ccp <- Maaslin2(
#   input_data = otu, 
#   input_metadata = meta, 
#   output = "maaslin2_output.stool.ccpplus_vs_ccpneg", 
#   fixed_effects = c("ccp3"))
```

#### 

<br>

#### Sputum

``` r
# Pick relative abundances (compositional) and sample metadata
sample_type = "Sputum"
pseq <- physeq %>% 
  subset_samples(sample_type_16S == sample_type) %>%
  subset_samples(ccp3_group != "PosRA") #%>%
  #tax_glom(taxrank="Genus") %>%
  
pseq.rel <- pseq %>% 
  microbiome::transform("compositional") #%>%
  #core(detection = 0.01, prevalence = 0.50)

otu <- abundances(pseq.rel) %>% t() %>% as.data.frame()
meta <- meta(pseq.rel) %>% mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+")))
```

##### Alpha Diversity

``` r
# Shannon diversity plot
shannon <- meta %>%
  ggplot(aes(x=ccp3, y=ShannonH.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Shannon Diversity") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Sobs diversity plot
sobs <- meta %>%
  ggplot(aes(x=ccp3, y=Sobs.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Richness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Species Evenness plot
evenness <- meta %>%
  ggplot(aes(x=ccp3, y=ShannonE.Median)) +
  geom_violin(aes(fill = ccp3), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Evenness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Display three on one plot
ggarrange(shannon, sobs, evenness, ncol = 3, nrow = 1)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
# Perform tests
shan.stat.test <- wilcox.exact(ShannonH.Median ~ ccp3, data=meta, paired=F)
rich.stat.test <- wilcox.exact(Sobs.Median ~ ccp3, data=meta, paired=F)
even.stat.test <- wilcox.exact(ShannonE.Median ~ ccp3, data=meta, paired=F)

# Output p-value table
data.frame(
  "Measurement" = c("Shannon Diversity", "Species Richness", "Species Evenness"),
  "p.val" = c(
    paste("p = ", round(shan.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(rich.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(even.stat.test$p.value, 2), sep = ""))) %>% 
  knitr::kable()
```

| Measurement       | p.val    |
|:------------------|:---------|
| Shannon Diversity | p = 0.05 |
| Species Richness  | p = 0.12 |
| Species Evenness  | p = 0.08 |

##### PCoA

``` r
# Need core taxa to save time
pcoa_otu <- pseq %>% 
  microbiome::transform("compositional") %>%
  core(detection = 0.01, prevalence = 0.10)
core_taxa <- abbrev_taxa(pcoa_otu) # Extract shortened taxa names
pcoa_otu <- pcoa_otu %>%
  abundances() %>% 
  t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(pcoa_otu))

# Determine coordinates for samples
PCoA <- vegdist(pcoa_otu, method="bray") %>%
  cmdscale() %>%
  as.data.frame() %>%
  select(Dim1=`V1`, Dim2=`V2`)

# Get vectors for taxa
taxa_vectors <- envfit(ord = PCoA, env = pcoa_otu)
taxa_vector_coords <- taxa_vectors$vectors$arrows * sqrt(taxa_vectors$vectors$r)
taxa_vector_p.vals <- taxa_vectors$vectors$pvals
vector_df <- data.frame(p_val = taxa_vector_p.vals) %>%
  bind_cols(taxa_vector_coords) %>%
  rownames_to_column("Taxa") %>%
  filter(p_val <= 0.05) %>%
  arrange(p_val) %>% head(5)

# Add metadata to ordination values
pcoa_plot_df <- PCoA %>% 
  merge(pseq %>% sample_data() %>% as.data.frame(), by = 'row.names') %>%
  mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+"))) %>%
  column_to_rownames('Row.names')

# Ordination bi-plot
pcoa_plot_df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = ccp3)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_segment(data = vector_df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               colour = "black", stat = "identity", alpha = 0.7, inherit.aes = FALSE) +
  geom_text_repel(data = vector_df, aes(x = Dim1, y = Dim2, label = Taxa), inherit.aes = FALSE, size=3) +
  theme_bw() +
  scale_color_manual(values = c("gray10", "darkgoldenrod")) +
  ggtitle("PCoA, Beta Diversity", "Bray-Curtis") +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 12))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

##### PERMANOVA

**Homogeneity of dispersions assumption not met!** We should therefore
not report the p-value here and should opt for another test. Looking at
the dispersion plots, nonetheless, shows us that there is not a
significance in community composition between groups (even if we were to
use another test).

``` r
# PERMANOVA test using Bray-Curtis distance
set.seed(007) # Set seed for reproducibility - permutation-based test
permanova <- adonis2(otu ~ ccp3 + age + gender + race, data = meta, 
                     by="margin", permutations = 999, method = "bray")

permanova %>% as.data.frame() %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  dplyr::rename("p.val" = "Pr(>F)") %>%
  mutate(p.val = ifelse(p.val <= 0.05, paste("**", p.val, "****", sep=""), p.val)) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>%
  knitr::kable(align = 'ccccc')
```

|          | Df  | SumOfSqs |  R2  |  F   | p.val |
|:---------|:---:|:--------:|:----:|:----:|:-----:|
| ccp3     |  1  |   0.16   | 0.02 | 1.25 | 0.22  |
| age      |  1  |   0.19   | 0.02 | 1.49 | 0.14  |
| gender   |  1  |   0.20   | 0.02 | 1.54 | 0.11  |
| race     |  5  |   0.62   | 0.07 | 0.97 | 0.53  |
| Residual | 56  |   7.20   | 0.87 |  NA  |  NA   |
| Total    | 64  |   8.24   | 1.00 |  NA  |  NA   |

**Dispersions Plot**

``` r
# Plot dispersion distances for each "group"
beta_dispersion <- otu %>% vegdist(method = "bray") %>% betadisper(meta$ccp3)
plot(beta_dispersion, hull=FALSE, ellipse=TRUE)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

**Homogeneity of Dispersons**

``` r
# Hypothesis test
set.seed(007)
otu %>% vegdist() %>% betadisper(meta$ccp3) %>% permutest()
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq    F N.Perm Pr(>F)  
    ## Groups     1  0.043  0.0433 4.18    999  0.047 *
    ## Residuals 63  0.653  0.0104                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### 

<br> <br>

### Hypothesis 2

#### CCP(+) Nonconverters vs Converters

*There are microbiome differences between CCP+ subjects who do and do
not go on to develop clinical RA.*  
<br> For this section, we will compare CCP(+)-Nonconverters and
CCP(+)-Converters. We are testing the hypothesis that certain microbiota
associated with ‘real’ development of future RA.

<br>

#### Stool

``` r
# Pick relative abundances (compositional) and sample metadata
sample_type = "Stool"
pseq <- physeq %>% 
  subset_samples(sample_type_16S == sample_type) %>%
  subset_samples(ccp3_group != "PosRA") %>%
  subset_samples(ccp3_group != "NegControl")
pseq.rel <- pseq %>% 
  microbiome::transform("compositional") #%>%
otu <- abundances(pseq.rel) %>% t() %>% as.data.frame()
meta <- meta(pseq.rel) %>% mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+")))
```

##### Alpha Diversity

``` r
# Shannon diversity plot
shannon <- meta %>%
  ggplot(aes(x=ccp3_group, y=ShannonH.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Shannon Diversity") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Sobs diversity plot
sobs <- meta %>%
  ggplot(aes(x=ccp3_group, y=Sobs.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Richness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Species Evenness plot
evenness <- meta %>%
  ggplot(aes(x=ccp3_group, y=ShannonE.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Evenness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Display three on one plot
ggarrange(shannon, sobs, evenness, ncol = 3, nrow = 1)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
# Perform tests
shan.stat.test <- wilcox.exact(ShannonH.Median ~ ccp3_group, data=meta, paired=F)
rich.stat.test <- wilcox.exact(Sobs.Median ~ ccp3_group, data=meta, paired=F)
even.stat.test <- wilcox.exact(ShannonE.Median ~ ccp3_group, data=meta, paired=F)

# Output p-value table
data.frame(
  "Measurement" = c("Shannon Diversity", "Species Richness", "Species Evenness"),
  "p.val" = c(
    paste("p = ", round(shan.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(rich.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(even.stat.test$p.value, 2), sep = "")
  )
) %>% knitr::kable()
```

| Measurement       | p.val    |
|:------------------|:---------|
| Shannon Diversity | p = 0.55 |
| Species Richness  | p = 0.38 |
| Species Evenness  | p = 0.9  |

##### PCoA

``` r
# Need core taxa to save time
pcoa_otu <- pseq %>% 
  microbiome::transform("compositional") %>%
  core(detection = 0.01, prevalence = 0.10)
core_taxa <- abbrev_taxa(pcoa_otu) # Extract shortened taxa names
pcoa_otu <- pcoa_otu %>%
  abundances() %>% 
  t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(pcoa_otu))

# Determine coordinates for samples
PCoA <- vegdist(pcoa_otu, method="bray") %>%
  # Morisita requires integer data; morisita-horn can handle abundance
  cmdscale() %>%
  as.data.frame() %>%
  select(Dim1=`V1`, Dim2=`V2`)

# Get vectors for taxa
taxa_vectors <- envfit(ord = PCoA, env = pcoa_otu)
taxa_vector_coords <- taxa_vectors$vectors$arrows * sqrt(taxa_vectors$vectors$r)
taxa_vector_p.vals <- taxa_vectors$vectors$pvals
vector_df <- data.frame(p_val = taxa_vector_p.vals) %>%
  bind_cols(taxa_vector_coords) %>%
  rownames_to_column("Taxa") %>%
  filter(p_val <= 0.05) %>%
  arrange(p_val) %>% head(5)

# Add metadata to ordination values
pcoa_plot_df <- PCoA %>% 
  merge(pseq %>% sample_data() %>% as.data.frame(), by = 'row.names') %>%
  column_to_rownames('Row.names')

# Ordination bi-plot
pcoa_plot_df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = ccp3_group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_segment(data = vector_df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               #arrow = arrow(length = unit(0.2, "cm")), 
               colour = "black", stat = "identity", alpha = 0.7, inherit.aes = FALSE) +
  geom_text_repel(data = vector_df, #vjust = "inward", hjust = "inward",
                  aes(x = Dim1, y = Dim2, label = Taxa), 
                  inherit.aes = FALSE, size=3) +
  #coord_fixed() +
  #xlim(-1, 0.8) +
  #ylim(-0.6, 1) +
  theme_bw() +
  scale_color_manual(values = c("gray10", "darkgoldenrod")) +
  ggtitle("PCoA, Beta Diversity", "Bray-Curtis") +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 12))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

##### PERMANOVA

``` r
# PERMANOVA test using Bray-Curtis distance
set.seed(007) # Set seed for reproducibility - permutation-based test
permanova <- adonis2(otu ~ ccp3_group + age + gender + race, data = meta, 
                     by = "margin", permutations = 999, method = "bray")

permanova %>% as.data.frame() %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  dplyr::rename("p.val" = "Pr(>F)") %>%
  mutate(p.val = ifelse(p.val <= 0.05, paste("**", p.val, "****", sep=""), p.val)) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>%
  knitr::kable(align = 'ccccc')
```

|            | Df  | SumOfSqs |  R2  |  F   | p.val |
|:-----------|:---:|:--------:|:----:|:----:|:-----:|
| ccp3_group |  1  |   0.05   | 0.01 | 0.63 | 0.76  |
| age        |  1  |   0.10   | 0.02 | 1.24 | 0.25  |
| gender     |  1  |   0.09   | 0.02 | 1.11 | 0.36  |
| race       |  4  |   0.32   | 0.08 | 1.05 | 0.41  |
| Residual   | 46  |   3.54   | 0.86 |  NA  |  NA   |
| Total      | 53  |   4.14   | 1.00 |  NA  |  NA   |

**Dispersions Plot**

``` r
# Plot dispersion distances for each "group"
beta_dispersion <- otu %>% vegdist(method = "bray") %>% betadisper(meta$ccp3_group)
plot(beta_dispersion, hull=FALSE, ellipse=TRUE)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

**Homogeneity of Dispersons**

``` r
# Hypothesis test
set.seed(007)
otu %>% vegdist() %>% betadisper(meta$ccp3_group) %>% permutest()
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq    F N.Perm Pr(>F)
    ## Groups     1  0.014 0.01395 1.41    999   0.24
    ## Residuals 52  0.514 0.00988

#### 

<br>

#### Sputum

``` r
# Pick relative abundances (compositional) and sample metadata
sample_type = "Sputum"
pseq <- physeq %>% 
  subset_samples(sample_type_16S == sample_type) %>%
  subset_samples(ccp3_group != "PosRA") %>%
  subset_samples(ccp3_group != "NegControl")
  #tax_glom(taxrank="Genus") %>%

pseq.rel <- pseq %>% 
  #tax_glom(taxrank="Genus") %>%
  microbiome::transform("compositional") #%>%
  #core(detection = 0.01, prevalence = 0.50)
otu <- abundances(pseq.rel) %>% t() %>% as.data.frame()
meta <- meta(pseq.rel) %>% mutate(ccp3 = factor(ccp3, levels = c(0, 1), labels = c("-", "+")))
```

##### Alpha Diversity

``` r
# Shannon diversity plot
shannon <- meta %>%
  ggplot(aes(x=ccp3_group, y=ShannonH.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Shannon Diversity") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Sobs diversity plot
sobs <- meta %>%
  ggplot(aes(x=ccp3_group, y=Sobs.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Richness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Species Evenness plot
evenness <- meta %>%
  ggplot(aes(x=ccp3_group, y=ShannonE.Median)) +
  geom_violin(aes(fill = ccp3_group), color = "gray5", lwd = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", color = "gray5") +
  geom_jitter(aes(color = ccp3_group), width = 0.03, alpha = 0.8) +
  geom_line(aes(group = sample_id), alpha = 0.5) +
  ggtitle("Species Evenness") +
  theme_bw() +
  scale_fill_manual(values = c("gray10", "darkgoldenrod")) +
  scale_color_manual(values = c("gray60", "black")) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Display three on one plot
ggarrange(shannon, sobs, evenness, ncol = 3, nrow = 1)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
# Perform tests
shan.stat.test <- wilcox.exact(ShannonH.Median ~ ccp3_group, data=meta, paired=F)
rich.stat.test <- wilcox.exact(Sobs.Median ~ ccp3_group, data=meta, paired=F)
even.stat.test <- wilcox.exact(ShannonE.Median ~ ccp3_group, data=meta, paired=F)

# Output p-value table
data.frame(
  "Measurement" = c("Shannon Diversity", "Species Richness", "Species Evenness"),
  "p.val" = c(
    paste("p = ", round(shan.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(rich.stat.test$p.value, 2), sep = ""),
    paste("p = ", round(even.stat.test$p.value, 2), sep = "")
  )
) %>% knitr::kable()
```

| Measurement       | p.val    |
|:------------------|:---------|
| Shannon Diversity | p = 0.49 |
| Species Richness  | p = 0.53 |
| Species Evenness  | p = 0.77 |

##### PCoA

``` r
# Need core taxa to save time
pcoa_otu <- pseq %>% 
  microbiome::transform("compositional") %>%
  core(detection = 0.01, prevalence = 0.10)
core_taxa <- abbrev_taxa(pcoa_otu) # Extract shortened taxa names
pcoa_otu <- pcoa_otu %>%
  abundances() %>% 
  t() %>% as.data.frame() %>%
  set_colnames(abbrev_taxa(pcoa_otu))

# Determine coordinates for samples
PCoA <- vegdist(pcoa_otu, method="bray") %>%
  # Morisita requires integer data; morisita-horn can handle abundance
  cmdscale() %>%
  as.data.frame() %>%
  select(Dim1=`V1`, Dim2=`V2`)

# Get vectors for taxa
taxa_vectors <- envfit(ord = PCoA, env = pcoa_otu)
taxa_vector_coords <- taxa_vectors$vectors$arrows * sqrt(taxa_vectors$vectors$r)
taxa_vector_p.vals <- taxa_vectors$vectors$pvals
vector_df <- data.frame(p_val = taxa_vector_p.vals) %>%
  bind_cols(taxa_vector_coords) %>%
  rownames_to_column("Taxa") %>%
  filter(p_val <= 0.05) %>%
  arrange(p_val) %>% head(5)

# Add metadata to ordination values
pcoa_plot_df <- PCoA %>% 
  merge(pseq %>% sample_data() %>% as.data.frame(), by = 'row.names') %>%
  column_to_rownames('Row.names')

# Ordination bi-plot
pcoa_plot_df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = ccp3_group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_segment(data = vector_df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), 
               #arrow = arrow(length = unit(0.2, "cm")), 
               colour = "black", stat = "identity", alpha = 0.7, inherit.aes = FALSE) +
  geom_text_repel(data = vector_df, #vjust = "inward", hjust = "inward",
                  aes(x = Dim1, y = Dim2, label = Taxa), 
                  inherit.aes = FALSE, size=3) +
  #coord_fixed() +
  #xlim(-1, 0.8) +
  #ylim(-0.6, 1) +
  theme_bw() +
  scale_color_manual(values = c("gray10", "darkgoldenrod")) +
  ggtitle("PCoA, Beta Diversity", "Bray-Curtis") +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 12))
```

![](altra-inference_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

##### PERMANOVA

``` r
# PERMANOVA test using Bray-Curtis distance
set.seed(007) # Set seed for reproducibility - permutation-based test
permanova <- adonis2(otu ~ ccp3_group + age + gender + race, data = meta, 
                     by="margin", permutations = 999, method = "bray")

permanova %>% as.data.frame() %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  dplyr::rename("p.val" = "Pr(>F)") %>%
  mutate(p.val = ifelse(p.val <= 0.05, paste("**", p.val, "****", sep=""), p.val)) %>% 
  mutate_if(is.numeric, ~round(., 2)) %>%
  knitr::kable(align = 'ccccc')
```

|            | Df  | SumOfSqs |  R2  |  F   | p.val |
|:-----------|:---:|:--------:|:----:|:----:|:-----:|
| ccp3_group |  1  |   0.03   | 0.00 | 0.20 | 1.00  |
| age        |  1  |   0.14   | 0.02 | 0.98 | 0.43  |
| gender     |  1  |   0.17   | 0.03 | 1.18 | 0.27  |
| race       |  4  |   0.54   | 0.09 | 0.94 | 0.55  |
| Residual   | 38  |   5.43   | 0.86 |  NA  |  NA   |
| Total      | 45  |   6.30   | 1.00 |  NA  |  NA   |

**Dispersions Plot**

``` r
# Plot dispersion distances for each "group"
beta_dispersion <- otu %>% vegdist(method = "bray") %>% betadisper(meta$ccp3_group)
plot(beta_dispersion, hull=FALSE, ellipse=TRUE)
```

![](altra-inference_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

**Homogeneity of Dispersons**

``` r
# Hypothesis test
set.seed(007)
otu %>% vegdist() %>% betadisper(meta$ccp3_group) %>% permutest()
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq    F N.Perm Pr(>F)
    ## Groups     1  0.003 0.00253 0.21    999   0.66
    ## Residuals 44  0.529 0.01203

#### 

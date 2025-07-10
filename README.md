# Tutorial for metabarcoding analyses with 16S rRNA gene amplicons

## 1. Getting ready to start with R

### 1.1 Overall workflow
![alt text](Workflow.PNG)

### 1.2 Package installation
First, some packages need to be installed.

For some of the packages, the installation needs to go through BiocManager
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
```

For the whole workflow, three main packages are necessary: [phyloseq](https://joey711.github.io/phyloseq/) , [microbiome](https://microbiome.github.io/tutorials/) and [vegan](https://www.rdocumentation.org/packages/vegan/versions/2.6-4)

```
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
install.packages(“vegan”)
```

For the whole decontamination of the dataset: 
```
BiocManager::install("decontam") #for decontamination of the dataset
```

For discriminant analyses
```
BiocManager::install("microbiomeMarker")
BiocManager::install("DESeq2")
install.packages("metacoder")
```

For the graphical outputs
```
install.packages("ggplot2")
install.packages("randomcoloR")
install.packages("ggh4x")
install.packages("hrbrthemes")
install.packages("cowplot")
install.packages("RColorBrewer")
install.packages("ggpubr")
install.packages("ggrepel")
```

For data frame manipulations
```
install.packages("stringi")
install.packages("rlang")
install.packages("Rcpp")
install.packages(“dplyr”)
install.packages("reshape2")
```

For additional statistical tests 
```
install.packages("agricolae")
install.packages(“RVAideMemoire”)
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages(“ape”)
install.packages("rcompanion")
install.packages("multcompView")
```

### 1.3 Packages loading
Now, let's load all the packages we need for the moment
```
library(stringi)
library(vegan)
library(Rcpp)
library(ggplot2)
library(randomcoloR)
library(rlang)
library(phyloseq)
library(ape)
library(dplyr)
library(agricolae)
library(reshape2)
library(RVAideMemoire) 
library(microbiome)
library(hrbrthemes)
library(cowplot) 
#library(ggthemes)
library(RColorBrewer) 
library(decontam)
library(ggrepel) 
library(rcompanion)
library(pairwiseAdonis) 
library(ggpubr)
library(ggh4x)
library(rcompanion)
library(multcompView)
```

### 1.4 Setting and preparing your working directory
```
setwd("your_path_to_insert_here")
```
We will now make subfolders where all the output from our analysis will be saved
```
dir.create("1_Data_prep_results")
dir.create("2_Alpha_div_results")
dir.create("3_Beta_div_results")
dir.create("4_Compositional_results")
dir.create("5_Core_community_results")
```


## 2. Preparation of the 16S metabarcoding dataset for a ready-to-go analysis

### 2.1. Import the dataset
First, we will import 3 ".csv" data tables 
- The "ASV_table" which includes all the numbers of reads of each ASVs in every samples
- The "Taxonomy_table" which includes all taxonomic information of each ASVs (usually levels such as Phylum, Class, Order, Family, Genus)
- The "Metadata_table" which includes all supplementary information on the sample (e.g. sampling date, location, etc)
The ".csv" files should be delimited with ;

```
ASV_table = read.csv(file = "ASV_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(ASV_table)

TAX_table = read.csv(file = "Taxonomy_table.csv" , sep = ";" , header = T , row.names = 1)
TAX_table = as.matrix(TAX_table)
dim(TAX_table)

META_table = read.csv(file = "Metadata_table.csv" , header = TRUE , sep = ";" , row.names = 1)
dim(META_table)
```

Using the dim() make sure that the dimension of each dataset is correct. 
- The number of rows should be the same for the ASV_table and Taxonomy_table
- The number of columns should be the same for the ASV_table and Metadata_table

Optional: if you have a tree file you can also import it for phyloseq

```
TREE_file = read_tree("Tree-file.nwk.nhx")
TREE_file
```

### 2.2. Create your main phyloseq object (output: physeq_raw)
Phyloseq is a package that allows you to generate a single object combining the ASV_table, the TAX_table, and the META_table (and optionally the TREE_file).

Using a phyloseq object makes your life more simple, as every manipulation of the dataset will be done on each tables. 

For example, removing a sample from the phyloseq object will remove it simulteanously from the ASV_table and META_table. Same idea for removing an ASV.
```
ASV_phylo = otu_table(ASV_table, taxa_are_rows = TRUE)
dim(ASV_phylo)

TAX_phylo = tax_table(TAX_table)
dim(TAX_phylo)

META_phylo = sample_data(META_table)
dim(META_phylo)

physeq_raw = phyloseq(ASV_phylo, TAX_phylo, META_phylo)
physeq_raw

```
Optional: your TREE_file can also be added to the phyloseq object
```
TREE_phylo = phy_tree(TREE_file)
physeq_raw_tree = merge_phyloseq(physeq_raw, TREE_phylo)
physeq_raw_tree
```
Your raw phyloseq object is now ready. 

Try to have a closer look at what's inside! 
Check how many ASVs and samples you have in total


<details>
  <summary>See the answer</summary>
  
  ```
ntaxa(physeq_raw)
nsamples(physeq_raw)
  ```
We have 1181 ASVs in total and 81 samples
</details>

Check how many ASVs per sample you have 


<details>
  <summary>See the answer</summary>
  
  ```
taxa_sums(physeq_raw)
min(taxa_sums(physeq_raw))
max(taxa_sums(physeq_raw))
mean(taxa_sums(physeq_raw))
  ```
We have around 652 ASVs on average. 53107 in max, and 81 in min (probably an extraction blank or negative control) 
</details>

### 2.2. Decontaminate your dataset using the decontam package (output: physeq_decontam)
You will need first to create a data frame "Libsize_table" to inspect your library size
 ```
Libsize_table <- as.data.frame(sample_data(physeq_raw)) 
Libsize_table$LibrarySize <- sample_sums(physeq_raw)
Libsize_table <- Libsize_table[order(Libsize_table$LibrarySize),]
Libsize_table$Index <- seq(nrow(Libsize_table))
Libsize_table

plot_libsize = ggplot(data=Libsize_table, aes(x=Index, y=LibrarySize, color=Sample_or_Control))
plot_libsize = plot_libsize + geom_point()
plot_libsize

ggsave(filename = "Plot_libsize.pdf", 
       plot = plot_libsize, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1_Data_prep_results")
 ```


<details>
  <summary>See figure</summary>
  
![alt text](1_Data_prep_results/Plot_libsize.png)

</details>

Based on the figure, we can observe that the control samples (negative controls and extraction blanks) have a library size significantly lower compared to most of our true samples.

Let's now identify the contaminants using the prevalence method
 
 ```
sample_data(physeq_raw)$is.neg <- sample_data(physeq_raw)$Sample_or_Control == "Control"

contam_prev05 <- isContaminant(physeq_raw, method="prevalence", neg="is.neg", threshold=0.5)
contam_prev05 = cbind(as.data.frame(tax_table(physeq_raw)) , contam_prev05)
contam_prev05

 ```
Based on these results, how many contaminants were identified? Are they abundant in the dataset? 

<details>
  <summary>See the answer</summary>
  
  ```
table(contam_prev05$contaminant)
  ```
We have identified 5 ASVs as contaminants. 
  ```
subset(contam_prev05, contaminant == "TRUE")
  ```
Most of them have a low prevalence in our dataset (which is good)
</details>

Now, let's save this information in a csv file, if we need to go back later on these contamination results 
  ```
write.csv(contam_prev05, file.path("./1_Data_prep_results" , "Contamination_table_prev05.csv"))
  ```

We can now make a plot of the prevalence of the contaminants ASV in the controls and true samples

We have to make a phyloseq object of presence-absence in negative controls and true samples ...
  ```
ps.pa <- transform_sample_counts(physeq_raw, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True sample", ps.pa)
 ```
... and then make a data.frame of prevalence in positive and negative samples
 ```
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam_prev05$contaminant)
 ```
The prevalence can be plotted as follows
 ```
plot_prevalence = ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) 
plot_prevalence = plot_prevalence + geom_point() 
plot_prevalence = plot_prevalence + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
plot_prevalence
 ```
Conclusion: the contaminants are mostly prevalent in the negative controls and not in the true samples! The results are ok!

We can save this plot as a figure to confirm the quality of our decontamination
 ```
ggsave(filename = "Plot_prevalence.pdf", 
       plot = plot_prevalence, 
       device = "pdf" , 
       width = 15 , height = 10, units = "cm", 
       path = "./1_Data_prep_results")
 ```

<details>
  <summary>See figure</summary>
  
![alt text](1_Data_prep_results/Plot_prevalence.png)

</details>

We can finally make our phyloseq object without these contaminants
 ```
physeq_decontam = prune_taxa(!contam_prev05$contaminant, physeq_raw)
physeq_decontam
 ```
Do you confirm that the right number of contaminating ASVs has been removed? 
<details>
  <summary>See the answer</summary>
  
  ```
ntaxa(physeq_raw)
ntaxa(physeq_decontam)
  ```
Yes, 5 ASVs were removed!
</details>



### 2.3. Filter additional non-prokaryotic taxa by removing ASVs with 16S from Eukaryotes, chloroplasts, and mitochondria

  ```
physeq_filtered = subset_taxa(physeq_decontam, 
                              (Kingdom != "Eukaryota" | is.na(Kingdom)) &
                              (Order!= "Bacteria | Cyanobacteria | Cyanobacteriia | Chloroplast"| is.na(Order)) &
                              (Family != "Bacteria | Proteobacteria | Alphaproteobacteria | Rickettsiales | Mitochondria"| is.na(Family)))

physeq_filtered
  ```
How many non-prokaryotic ASVs were filtered in total? 

<details>
  <summary>See the answer</summary>
  
  ```
ntaxa(physeq_decontam) - ntaxa(physeq_filtered)
  ```
14 ASVs in total
</details>

### 2.4. Checking the rarefaction curves and removing samples with not enough reads 
Rarefaction curves are a good representation to verify if all the diversity is covered in each sample. 

On the x-axis, we have the number of reads, on the y-axis, the number of new ASVs discovered. 

If we reach a plateau, it means that no more new ASVs are identified, even if we increase or sample size (number of reads).

  ```
as.data.frame(t(otu_table(physeq_filtered)))

rarecurve(as.data.frame(t(otu_table(physeq_filtered))), step = 20, cex = 0.5)
# Save the plot manually to ./1. Data_prep results
  ```

<details>
  <summary>See figure</summary>
  
![alt text](1_Data_prep_results/Rarefaction_curves.png)

</details>



Future update soon with the package [inext](https://johnsonhsieh.github.io/iNEXT/)

Based on these rarefaction curves, we can now set a minimum reads threshold. 

The samples below this threshold will be excluded from our analysis, as the sequencing depth is expected to be not enough. 


  ```
min_reads_threshold = 2000
physeq_subsampled = prune_samples(sample_sums(physeq_filtered)>=min_reads_threshold, physeq_filtered)
physeq_subsampled
  ```
With this code, we also remove the negative controls and extraction blanks. 

How many samples were removed in total? 
<details>
  <summary>See the answer</summary>
  
  ```
nsamples(physeq_filtered) - nsamples(physeq_subsampled)
  ```
3 samples in total (only the negative controls and extraction blanks, in this case)

</details>

physeq_subsampled is your phyloseq object ready to be transformed for the alpha-diversity, beta-diversity, and compositional analyses! 
Save it as a RDS object to avoid all these steps in a future session! 

  ```
saveRDS(physeq_subsampled, "Physeq_subsampled.RDS")
  ```
When you start a new session (a few days later for example), here is how to load the phyloseq object again

  ```
physeq_subsampled = readRDS("Physeq_subsampled.RDS")
physeq_subsampled
  ```

## 3. Before starting starting the analyses, let's see why the dataset should be transformed

Rarefying your datasets means that in every sample, you will randomly remove ASVs counts until they reach the same number of ASVs of the sample with the lower size. 

As a result, all your samples will have the same read numbers, which correspond to the minimum library size

Why do we remove such information in our dataset? To make diversity comparison possible between samples, as some samples end up with more reads than others after the sequencing. 

How does it impact my dataset? Rarefaction usually removes very rare ASVs from your datasets, as their is a high chance that their original counts (e.g. 1 or 2) goes to 0.

This approach is currently criticized and debated. Some groups argue that rare ASVs are essential to consider, while others find that the effect of the rarefaction is almost non-existent

For more details, see these papers:

McMurdie, P.J. and Holmes, S. (2014) Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Comput Biol 10: e1003531.

Gloor, G.B., Macklaim, J.M., Pawlowsky-Glahn, V., and Egozcue, J.J. (2017) Microbiome datasets are compositional: and this is not optional. Front Microbiol 8: 2224.

Schloss, P.D. (2023) Waste not, want not: revisiting the analysis that called into question the practice of rarefaction. mSphere 9: e00355-23.
 

<details>
  <summary>See my opinion here on the question </summary>
  
I think that the choice of the rarefaction transformation depends on the look of the rarefaction curves.

If most of your samples almost reach the plateau, but not entirely, then rarefaction is needed for the alpha-diversity analyses.

If all samples clearly showed a full plateau, it means that nothing is missing in your samples, and then there is no need to worry about the difference in library size.

In this rare situation, there is no need to rarefy your dataset, and the alpha-diversity (especially the richness) can be plotted with the real counts from your dataset, with all the rare ASVs.

For all the rest of the analyses (beta-diversity and composition), I prefer to keep all the information and just work with normalization to the sum (i.e., transforming your counts into percentages).

</details>


## 4. Alpha diversity analyses

### 4.1. Create a phyloseq rarefied object specifically for alpha-div analyses

As explained just before, we will rarefy our dataset, so all samples will end up with the same number of reads in total.
  ```
min(sample_sums(physeq_subsampled))

physeq_rarefied = rarefy_even_depth(physeq_subsampled)
physeq_rarefied
  ```
How many ASVs per sample do we end up with? 

<details>
  <summary>See the answer</summary>
  
  ```
sample_sums(physeq_rarefied)
  ```
2674 ASVs per sample

</details>


### 4.2. Consider the factors of comparison for this analysis and prepare your color vectors

Before plotting the alpha diversity, let's go back to the current metadata

  ```
metadata_subsampled = sample_data(physeq_subsampled)
metadata_subsampled
  ```

What factor of interest is important to compare in this dataset? How many levels (groups) there is in these factors?
<details>
  <summary>See the answer</summary>

- For the temporal changes: the sampling date

```
factor(metadata_subsampled$Month)
```
We have 6 distinct sampling dates

- For the spatial comparison: the sampling site

```
factor(metadata_subsampled$Site)
```
We have 5 distinct sampling sites

</details>

Let's create the color vector for the sampling dates. 

For temporal changes, I like to have my colors changing gradually from a cold color for the cold months, to hot colors for the warmer months

```
color_months = c("Violet","SkyBlue", "LightGreen", "Yellow", "Orange", "Red")

```

### 4.3. Create a data frame with the results of the alpha-diversity indices and the metadata

The function estimate_richness allows you to calculate many different alpha-diversity indices. 

Here we will keep only the observed richness, the Shannon index, and the Chao1 (estimated richness)

```
data_alpha = estimate_richness(physeq_rarefied , measures = c("Observed", "Shannon", "Chao1"))
data_alpha
```

Unfortunately, phyloseq does not provide the Pielou index, so we have to add it to calculate it with this formula
```
Pielou = data_alpha$Shannon / log(data_alpha$Observed)
Pielou
```

We can now add the Pielou index and the factors of interest to the data frame
```
data_alpha_all = cbind(metadata_subsampled[, c("Site","Month","Site_code","Month_code")], data_alpha , Pielou)
data_alpha_all
```

What if I don't want to include all the samples in my analysis? (just only the site S1, for example)

<details>
  <summary>See the answer</summary>
To plot your alpha diversity results if you want to include just a specific selection of samples, use the subset function

```
data_alpha_all_subset = subset(data_alpha_all, Site_code == "S1"  )
data_alpha_all_subset 
```

</details>

### 4.4. Generate a single figure for the 3 alpha-diversity indices together

To have the three diversity indices (Shannon, Chao1, and Pielou) in different facets, we need to rearrange the data.frame with a single column for all values

```
data_long = melt(data_alpha_all, id.var=c("Site","Month","Site_code","Month_code"))
data_long
```

We can now make boxplots using ggplot2 with facet_grid for the figure

```
plot_alpha_all = ggplot(data_long, aes(x = Month_code, y = value))
plot_alpha_all = plot_alpha_all + geom_point(aes(pch = Site, fill = Month), size = 2, alpha = 0.8, pch = 21)
plot_alpha_all = plot_alpha_all + geom_boxplot(aes(fill = Month), alpha = 0.8, size = 0.8)
plot_alpha_all = plot_alpha_all + scale_fill_manual(values = color_months)
plot_alpha_all = plot_alpha_all + theme_bw(base_size = 15)  + theme(legend.position="left")
plot_alpha_all = plot_alpha_all + facet_grid(variable~Site, scales = "free")
plot_alpha_all = plot_alpha_all + scale_shape_manual(values = c(21,22,23,24,25))
plot_alpha_all = plot_alpha_all + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_alpha_all = plot_alpha_all +theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_alpha_all


ggsave(filename = "Plot_alpha_all.pdf", 
       plot = plot_alpha_all, 
       device = "pdf" , 
       width = 35 , height = 25, units = "cm", 
       path = "./2_Alpha_div_results")
```
<details>
  <summary>See figure</summary>
  
![alt text](2_Alpha_div_results/Plot_alpha_all.png)

</details>



### 4.5. Test the normality of your distribution for each index

```
shapiro_data_shannon = shapiro.test(data_alpha_all$Shannon)
shapiro_data_shannon

shapiro_data_chao1 = shapiro.test(data_alpha_all$Chao1)
shapiro_data_chao1

shapiro_data_pielou = shapiro.test(data_alpha_all$Pielou)
shapiro_data_pielou
```

Let's make a table summarizing these results

```
shapiro_data_alpha <- matrix(nrow = 3 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_alpha) = c("W","p-value")
rownames(shapiro_data_alpha) = c("Shannon","Chao1","Pielou")

shapiro_data_alpha[,1] <- c(shapiro_data_shannon$statistic,
                               shapiro_data_chao1$statistic,
                               shapiro_data_pielou$statistic)

shapiro_data_alpha[,2]<- c(shapiro_data_shannon$p.value,
                           shapiro_data_chao1$p.value,
                           shapiro_data_pielou$p.value)

shapiro_data_alpha

write.csv(shapiro_data_alpha, file.path("./2_Alpha_div_results" , "Shapiro_data_alpha.csv"))

```
Based on these results, all p-values are < 0.05, implying that the distribution of the data is significantly different from normal distribution. 

In other words, we cannot assume normality. The Kruskal-Wallis test needs to be used (instead of a parametric ANOVA)


### 4.6. Test the differences of alpha-diversity according to your factors with analyses of variances

Let's test the differences between months and sites, with the Kruskal-Wallis test

```
data_kruskal_Shannon_month = kruskal.test(Shannon ~ Month, data_alpha_all)
data_kruskal_Shannon_month 

data_kruskal_Chao1_month = kruskal.test(Chao1 ~ Month, data_alpha_all)
data_kruskal_Chao1_month 

data_kruskal_Pielou_month = kruskal.test(Pielou ~ Month, data_alpha_all)
data_kruskal_Pielou_month 
```
We can now create a table summarizing this information all together

```
data_kruskal_alpha_month  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(data_kruskal_alpha_month) = c("Chi-square","Df","p-value")
rownames(data_kruskal_alpha_month) = c("Shannon","Chao1","Pielou")

data_kruskal_alpha_month

data_kruskal_alpha_month[,1] <- c(data_kruskal_Shannon_month$statistic,
                                  data_kruskal_Chao1_month$statistic,
                                  data_kruskal_Pielou_month$statistic)

data_kruskal_alpha_month[,2]<- c(data_kruskal_Shannon_month$parameter,
                                 data_kruskal_Chao1_month$parameter,
                                 data_kruskal_Pielou_month$parameter)

data_kruskal_alpha_month[,3]<- c(data_kruskal_Shannon_month$p.value,
                                 data_kruskal_Chao1_month$p.value,
                                 data_kruskal_Pielou_month$p.value)

data_kruskal_alpha_month

write.csv(data_kruskal_alpha_month, file.path("./2_Alpha_div_results" , "Data_kruskal_alpha_month.csv"))
```
Conclusion: For each alpha-diversity index, there are significant differences across time.

What about the spatial differences?

<details>
  <summary>See the answer</summary>

```
data_kruskal_Shannon_site = kruskal.test(Shannon ~ Site, data_alpha_all)
data_kruskal_Shannon_site 

data_kruskal_Chao1_site = kruskal.test(Chao1 ~ Site, data_alpha_all)
data_kruskal_Chao1_site 

data_kruskal_Pielou_site = kruskal.test(Pielou ~ Site, data_alpha_all)
data_kruskal_Pielou_site 



data_kruskal_alpha_site  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(data_kruskal_alpha_site) = c("Chi-square","Df","p-value")
rownames(data_kruskal_alpha_site) = c("Shannon","Chao1","Pielou")

data_kruskal_alpha_site

data_kruskal_alpha_site[,1] <- c(data_kruskal_Shannon_site$statistic,
                                  data_kruskal_Chao1_site$statistic,
                                  data_kruskal_Pielou_site$statistic)

data_kruskal_alpha_site[,2]<- c(data_kruskal_Shannon_site$parameter,
                                 data_kruskal_Chao1_site$parameter,
                                 data_kruskal_Pielou_site$parameter)

data_kruskal_alpha_site[,3]<- c(data_kruskal_Shannon_site$p.value,
                                 data_kruskal_Chao1_site$p.value,
                                 data_kruskal_Pielou_site$p.value)

data_kruskal_alpha_site

write.csv(data_kruskal_alpha_site, file.path("./2_Alpha_div_results" , "Data_kruskal_alpha_site.csv"))

```
Conclusion: for the three indices there are no differences between the sites

</details>

### 4.7. Test the differences between groups with pairwise comparisons

Since we performed a Kruskal-Wallis test, we will now conduct a Wilcoxon test for the pairwise comparison between each month. 

Let's try with the Shannon index
```
data_wilcox_shannon_month = pairwise.wilcox.test(data_alpha_all$Shannon, data_alpha_all$Month, p.adjust.method = "bonferroni")
data_wilcox_shannon_month
data_wilcox_shannon_month$p.value
```

We can now create a table summarizing these results. 
This table should include the lowercase indices ("a", "b", "c", ...) allowing us to easily identify the significant groups with the boxplots.

```
data_wilcox_shannon_month_full = fullPTable(data_wilcox_shannon_month$p.value)
data_wilcox_shannon_month_full

indices_wilcox_shannon_month = multcompLetters(data_wilcox_shannon_month_full,
                                           compare="<",
                                           threshold=0.05,
                                           Letters=letters,
                                           reversed = FALSE)
indices_wilcox_shannon_month$Letters



data_wilcox_shannon_month_full2 = cbind(data_wilcox_shannon_month_full, indices_wilcox_shannon_month$Letters )
data_wilcox_shannon_month_full2

write.csv(data_wilcox_shannon_month_full2, file.path("./2_Alpha_div_results" , "Data_wilcox_shannon_month.csv"))
```

Now let's do the same with the Pielou and the Chao1 indices!


<details>
  <summary>See the answer</summary>

Just copy-paste the same script and replace the index name using ctrl+F !!

</details>



## 5. Beta diversity analyses

### 5.1. Create a phyloseq compositional objects specifically for beta-div and compositional analyses

We need now to create a phyloseq object without rarefaction. Instead, the dataset will be transformed into compositional values. 

```
physeq_compo <- transform(physeq_subsampled, "compositional")   
physeq_compo
```

How does the data look like in the ASV_table now? What is the sum of all ASVs in a sample?

<details>
  <summary>See the answer</summary>
  
```
otu_table(physeq_compo)
sample_sums(physeq_compo)
```
All values are now in percentage (100% --> 1)
</details>


### 5.2. Consider the factors of comparison for this analysis and prepare your color and shape vectors

We want to keep the same color code used in the alpha-diversity analysis

```
color_months = c("Violet","SkyBlue", "LightGreen", "Yellow", "Orange", "Red")
```

But for the site, we want to distinguish them based on another  graphical option. Let's use different pch shapes for the geom_point. 

```
shape_sites = c(21,22,23,24,25)
```

### 5.3. NMDS analysis

We can now make our NMDS analysis with all samples, and calculate the 2D stress.

```
nmds <- ordinate(physeq_compo, "NMDS", "bray")
nmds$stress
```

We can create a data frame with the coordinates of the two NMDS axes. 

```
data_nmds <- plot_ordination(physeq_compo, nmds, type = "Samples", justDF = TRUE)
data_nmds
```

Now, let's plot the NMDS analysis

```
plot_nmds = ggplot(data_nmds, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = Month, pch = Site), size = 5, alpha = 0.8)
plot_nmds = plot_nmds + scale_shape_manual(values = shape_sites)
plot_nmds = plot_nmds + scale_fill_manual(values = color_months)
plot_nmds = plot_nmds + theme_bw(base_size = 15)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7)))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.16", x = -1, y = 1)
plot_nmds

ggsave(filename = "Plot_NMDS.pdf", 
       plot = plot_nmds, 
       device = "pdf" , 
       width = 30 , height = 20, units = "cm", 
       path = "./3_Beta_div_results")
```
<details>
  <summary>See figure</summary>
  
![alt text](3_Beta_div_results/Plot_NMDS.png)

</details>

### 5.5 Multivariate test to test the effect of each factor on the overall variance of the beta-diversity

Let's calculate the dissimilarity matrix between every sample

```
data_distbeta = as.matrix(distance(physeq_compo, method="bray"))
data_distbeta
```

Let's check again the metadata

```
metadata_subsampled <- as(sample_data(physeq_compo), "data.frame")
metadata_subsampled
```
We have two factors of comparison that can interact with each other (Month and Site). So a two-way PERMANOVA test (aka two-way Adonis test) is required. 

```
data_permnova = adonis2(data_distbeta ~ Month*Site, data = metadata_subsampled)
data_permnova

write.csv(as.data.frame(data_permnova), 
          file.path("./3_Beta_div_results" , "Data_twoway_permnova.csv"))
 
```
Conclusion: there are significant differences according to the sampling time and site. Both factors interact with each other in a significant way. 

We can now perform multivariate pairwise comparisons to check which months are different from each other. Same for the sites. 

```
data_pairwiseadonis_month = pairwise.adonis(data_distbeta, metadata_subsampled$Month)
data_pairwiseadonis_month

write.csv(data_pairwiseadonis_month, 
          file.path("./3_Beta_div_results" , "Pairwise_adonis_month.csv"))

data_pairwiseadonis_site = pairwise.adonis(data_distbeta, metadata_subsampled$Site)
data_pairwiseadonis_site

write.csv(data_pairwiseadonis_month, 
          file.path("./3_Beta_div_results" , "Pairwise_adonis_site.csv"))

```

### 5.4. Beta-diversity dispersion analysis 

The beta-diversity dispersion analysis allows you to calculate the dispersion within groups of samples. It can help to confirm visual observation based on the NMDS plots, when you observe that some groups have replicates more dispersed than others. 

Let's calculate the beta-dispersion based on the distance to centroids within each group (1 group of replicate is one sampling at a specific time on the specific site) 

```
betadisper_result <- betadisper(distance(physeq_compo, method="bray"), data_nmds$Month_site)
betadisper_result

dist_centroid = betadisper_result$distances
dist_centroid
```

We can now combine these results with the metadata and plot these distances to centroids for 

```
data_betadisper = cbind(data_nmds, dist_centroid)
data_betadisper

plot_betadisper = ggplot(data_betadisper, aes(x = Month_code, y = dist_centroid))
plot_betadisper = plot_betadisper + geom_point(aes(pch = Site, fill = Month), size = 2, alpha = 0.8, pch = 21)
plot_betadisper = plot_betadisper + geom_boxplot(aes(fill = Month), alpha = 0.8, size = 0.8)
plot_betadisper = plot_betadisper + facet_grid(Site~., scales = "free")
plot_betadisper = plot_betadisper + scale_fill_manual(values = color_months)
plot_betadisper = plot_betadisper + theme_bw(base_size = 15)  + theme(legend.position="left")
plot_betadisper = plot_betadisper + scale_shape_manual(values = c(21,22,23,24,25))
plot_betadisper = plot_betadisper + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_betadisper = plot_betadisper + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_betadisper

ggsave(filename = "Plot_betadisper.pdf", 
       plot = plot_betadisper, 
       device = "pdf" , 
       width = 15 , height = 30, units = "cm", 
       path = "./3_Beta_div_results")

```
<details>
  <summary>See figure</summary>
  
![alt text](3_Beta_div_results/Plot_betadisper.png)

</details>

Ok, what about the statistical differences? 

<details>
  <summary>See the answer</summary>

You can use the same script as the stats for the alpha-diversity !

</details>


### 5.5. db-RDA analysis to identify environmental parameters as explanatory factors of the beta-diversity

Let's perform a db-RDA with the vegan package. For this, we will need the distance matrix of the compositional dataset and the metadata which include the environmental parameters

Good news, we already have these two data frames! Can you find them back from above? 

<details>
  <summary>See the answer</summary>

```
data_distbeta
metadata_subsampled
```
</details>

Now let's consider the important environmental parameters for our analysis. In this case, our choice goes to  the temperature, salinity, silicates, phosphate, and the filtered fraction of the copper. We can perform the db-RDA analysis based on a model that include all of these parameters

```
dbrda_result <- capscale(data_distbeta~ Temperature + Salinity + SIOH4 + PO4 + Cu_F, metadata_subsampled)
dbrda_result
```

We can now test the significance of the model used

```
anova.cca(dbrda_result, step = 999)
```

Ok, the model is significant, now let's plot the db-RDA. For a good graphical output (using ggplot2) some data need to be extracted first. The coordinates of the two axes, but also the coordinates for the vectors related to each environmental parameter. 

```
dbrda_result_scores <- scores(dbrda_result, display = "sites")
dbrda_result_scores <- as.data.frame(dbrda_result_scores)
dbrda_result_scores

dbrda_result_scores_metadata <- cbind(dbrda_result_scores, metadata_subsampled)
dbrda_result_scores_metadata

dbrda_result_arrow <- scores(dbrda_result, display = "bp")
dbrda_result_arrow <- as.data.frame(dbrda_result_arrow)
dbrda_result_arrow
```

Based on the data extracted from the db-RDA analysis, we can now plot using ggplot2

```
plot_rda <- ggplot(dbrda_result_scores_metadata, aes(x = CAP1 , y = CAP2))
plot_rda <- plot_rda + geom_point(aes(pch = Site,  fill = Month), size = 5, alpha = 0.8)
plot_rda <- plot_rda + scale_shape_manual(values = shape_sites)
plot_rda <- plot_rda + scale_fill_manual(values = color_months)
plot_rda <- plot_rda + geom_segment(data = dbrda_result_arrow, aes(x=0, y=0, xend= CAP1, yend = CAP2), color = "black", arrow = arrow(length = unit(0.03, "npc")), linewidth = 1 , alpha = 0.7)
plot_rda <- plot_rda + theme_bw(base_size = 15)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_rda <- plot_rda + geom_text(data = dbrda_result_arrow, aes(x=1.1*CAP1, y=1.1*CAP2), label = row.names(dbrda_result_arrow), color = "black", cex = 6)
plot_rda

ggsave(filename = "Plot_rda.pdf", 
       plot = plot_rda, 
       device = "pdf" , 
       width = 25 , height = 20, units = "cm", 
       path = "./3_Beta_div_results")
```

<details>
  <summary>See figure</summary>
  
![alt text](3_Beta_div_results/Plot_rda.png)

</details>

### 5.6. Variance partitioning

The variance partitioning analysis allows you to determine the relative contribution of each environmental parameter (or group of parameters) as explanatory factors of the beta-diversity variance. 

As for the db-RDA, we will perform the analysis first with vegan using the distance matrix and the metadata

```
varpart_result <- varpart(data_distbeta, ~ Temperature + Salinity, ~Cu_F, ~PO4 + NO3, data=metadata_subsampled)
varpart_result
showvarparts(3, bg=2:4)
plot(varpart_result, bg=2:4)
```

<details>
  <summary>See figure and interpretation</summary>
  
![alt text](3_Beta_div_results/Plot_varpart.png)

X1 corresponds to temperature and salinity together

X2 corresponds to the copper-filtered

X3 corresponds to phosphate and silicates

The results on the Venn diagrams indicate the percentage explained by each group, independently or combined. 

The residual variance unexplained reach 70%

</details>


## 6. Compositional and differential analyses

### 6.1. Barplots at the family level

For compositional analyses, I have a personal preference for barplots instead of heatmaps or bubble plots, as I find it difficult to compare different values just based on colors (for the heatmaps) are circle area (bubble plots). 

One of the difficult parts of the compositional representation is related to the high number of variables. Of course, making barplots at the ASVs or genus level will make the figure very messy, with too many variables. One solution is to make barplots at the family levels, and combine all the rare families into a single group named "others". 

For this aggregation of the dataset at the family level we can use the microbiome package, to make a new phyloseq object with families instead of ASVs. 

We can combine all rare families below 3% together with the detection threshold

```
physeq_aggreg <- aggregate_rare(physeq_compo, level = "Family", detection = 3/100, prevalence = 0/100)
physeq_aggreg
```

Let's have a closer look at this phyloseq object, especially with the different families to plot

```
ntaxa(physeq_aggreg)
tax_table_family = as.data.frame(tax_table(physeq_aggreg))
tax_table_family$Family

```

We can create a color vector, based on the affiliation of the families. For example different orange colors for Bacteroidetes, green colors for Gammaproteobacteria, and blue colors for Alphaproteobacteria


```
color_families =  c( "Bacteria | Actinobacteria | Acidimicrobiia | Microtrichales | Microtrichaceae"="darkslateblue",                   
                   "Bacteria | Actinobacteria | Actinobacteria | Micrococcales | Intrasporangiaceae",                 
                   "Bacteria | Bacteroidetes | Bacteroidia | Chitinophagales | Saprospiraceae"="#FF914D",                       
                   "Bacteria | Bacteroidetes | Bacteroidia | Chitinophagales | unknown family"="#FFBD59",                       
                   "Bacteria | Bacteroidetes | Bacteroidia | Cytophagales | Amoebophilaceae"="#FF8200",                         
                   "Bacteria | Bacteroidetes | Bacteroidia | Cytophagales | Cyclobacteriaceae"="#FFC78E",                       
                   "Bacteria | Bacteroidetes | Bacteroidia | Flavobacteriales | Flavobacteriaceae"="#FFE38E",                   
                   "Bacteria | Bacteroidetes | Bacteroidia | Sphingobacteriales | NS11-12 marine group"="#E8B374",              
                   "Bacteria | Chloroflexi | Anaerolineae | Anaerolineales | Anaerolineaceae"="#FF66C4",                        
                   "Bacteria | Chloroflexi | Anaerolineae | Ardenticatenales | unknown family"="#FFAEE0",                       
                   "Bacteria | Chloroflexi | Anaerolineae | Caldilineales | Caldilineaceae"="#E4299C",                    
                   "Bacteria | Cyanobacteria | Oxyphotobacteria | Nostocales | Nostocaceae"="#FF3131",                          
                   "Bacteria | Cyanobacteria | Oxyphotobacteria | Nostocales | unknown family"="#C20F0F",                       
                   "Bacteria | Cyanobacteria | Oxyphotobacteria | Nostocales | Xenococcaceae"="#FD8181",                        
                   "Bacteria | Cyanobacteria | Oxyphotobacteria | Phormidesmiales | Phormidesmiaceae"="#F60000",                
                   "Bacteria | Cyanobacteria | Oxyphotobacteria | Synechococcales | Synechococcales Incertae Sedis"="#B40202", 
                   "Bacteria | Firmicutes | Clostridia | Clostridiales | Ruminococcaceae"="#D9D9D9", 
                   "Bacteria | Planctomycetes | Planctomycetacia | Pirellulales | Pirellulaceae"="#FFF847",                     
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Caulobacterales | Hyphomonadaceae"="#0097B2",             
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Micavibrionales | unknown family"="#0CC0DF",              
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales | Hyphomicrobiaceae"="#5CE1E6",   
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales | Rhizobiaceae"="#38B6FF",     
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae"="#5271FF",         
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodovibrionales | Kiloniellaceae"="#004AAD",     
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Sphingomonadales | Sphingomonadaceae"="#00FFE8",          
                   "Bacteria | Proteobacteria | Alphaproteobacteria | unknown order | unknown family"="#0D00FF", 
                   "Bacteria | Proteobacteria | Deltaproteobacteria | Bdellovibrionales | Bacteriovoracaceae"="yellow",        
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Alteromonadales | Alteromonadaceae"="#C1FF72",            
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Alteromonadales | Shewanellaceae"="#7ED957", 
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Cellvibrionales | Cellvibrionaceae"="#00BF63",     
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Cellvibrionales | Spongiibacteraceae"="#217F4D",          
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Ectothiorhodospirales | Ectothiorhodospiraceae"="#5EFFA8",
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Pseudomonadaceae"="#B5E8A6",            
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Thiohalorhabdales | Thiohalorhabdaceae"="#48FFAD",        
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Thiotrichales | Thiotrichaceae"="#57FF48",                
                   "Bacteria | Verrucomicrobia | Verrucomicrobiae | Verrucomicrobiales | DEV007"="#9D784D",     
                   "Bacteria | Verrucomicrobia | Verrucomicrobiae | Verrucomicrobiales | Rubritaleaceae"="#724C1F",             
                   "None | unknown | unknown | unknown | unknown"="#A6A6A6",                                                    
                   "Other"="gray" )
```

Now, let's make a modified version of the plot_bar function from phyloseq to remove the borders of each barplot

See: https://stackoverflow.com/questions/52747802/how-to-remove-very-thin-bar-plot-outline-border

```
plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

```

We can now make the barplots at the family level


```
barplot_family = plot_bar_2(physeq_aggreg, "Month_site_replicate", fill = "Family")
barplot_family = barplot_family + facet_grid(~Site_code + Month_code, scales = "free", space="free") 
barplot_family = barplot_family + scale_fill_manual(values = color_families)
barplot_family = barplot_family + theme(legend.text = element_text(size=9),
                                    axis.text.x = element_blank(),
                                    legend.title = element_blank(),
                                    axis.ticks.x=element_blank(), 
                                    axis.title = element_blank())
barplot_family = barplot_family + guides (fill = guide_legend(ncol = 1)) 
barplot_family
barplot_family = barplot_family + theme(panel.spacing = unit(0, "cm", data = NULL),panel.border = element_rect(color = "black", fill = NA, size = 0.9))
barplot_family = barplot_family +scale_y_continuous(expand = c(0,0))
barplot_family

ggsave(filename = "Barplot_family.pdf", 
       plot = barplot_family, 
       device = "pdf" , 
       width = 40 , height = 20, units = "cm", 
       path = "./4_Compositional_results")

```

<details>
  <summary>See figure</summary>
  
![alt text](4_Compositional_results/Barplot_family.png)

</details>



### 6.2. Differential analyses

Coming soon. See DeSeq2, LEfSe or Metacoder


## 7. Core community analyses

The core community is defined as the part of the microbial community always present in a group of samples, whatever the changes. 

In this study, we will focus on the spatial-temporal core community, always present whatever the site or the month. 


### 7.1. Choosing the prevalence threshold and preparing a new phyloseq object

We can use the microbiome package to make a new phyloseq object specifically subsetting this core community. 

A prevalence threshold of 90% will allow us to keep ASVs present at least in 90% of the samples. All other ASVs (not core) will be excluded from the dataset. 

The choice of this threshold depends on your dataset (the relation between your samples) and the number of samples

```
physeq_core <- core(physeq_compo, detection = 0.01/100, prevalence = 90/100)  ##phyloseq object of the core microbiota is obtained 
physeq_core
```

How many core ASVs do we have in total? In which class are they mainly affiliated?

<details>
  <summary>See the answer</summary>

```
ntaxa(physeq_core)
tax_table(physeq_core)
```
63 core ASVs in total, mostly affiliated to Alpha- and Gammaproteobacteria, and Bacteroidia

</details>


### 7.2. Plotting the relative abundance of the core community using barplots

The relative abundance of sequences you core ASV can be plotted using barplots at the family level. 

As there is only 63 core ASVs, the plots could also be plotted at the genus level and might not be too dense. 


```
physeq_core_family = tax_glom(physeq_core, taxrank = "Family")
physeq_core_family

barplot_core_family = plot_bar_2(physeq_core_family, "Month_site_replicate", fill = "Family") + theme_bw()
barplot_core_family = barplot_core_family + geom_bar(stat = "identity" , position="stack")
barplot_core_family = barplot_core_family + theme(legend.text = element_text(size=9),
                          axis.ticks.y = element_blank(),
                          legend.position="right",
                          axis.text.x = element_blank(),
                          legend.title = element_blank(),
                          axis.ticks.x=element_blank(), 
                          axis.title = element_blank())
barplot_core_family = barplot_core_family + scale_fill_manual(values = color_families)
barplot_core_family = barplot_core_family + facet_grid(~Site_code + Month_code, scales = "free", space = "free")
barplot_core_family =  barplot_core_family + guides(fill = guide_legend(ncol = 1))
barplot_core_family =  barplot_core_family + theme(axis.title = element_blank())
barplot_core_family = barplot_core_family  + theme(panel.spacing = unit(0, "cm", data = NULL),panel.border = element_rect(color = "black", fill = NA, size = 1))
barplot_core_family = barplot_core_family + scale_y_continuous(expand = c(0,0))
barplot_core_family


ggsave(filename = "Barplot_coreASVs_familylevel.pdf", 
       plot = barplot_core_family, 
       device = "pdf" , 
       width = 40 , height = 20, units = "cm", 
       path = "./5_Core_community_results")

```

<details>
  <summary>See figure</summary>
  
![alt text](5_Core_community_results/Barplot_coreASVs_familylevel.png)

</details>


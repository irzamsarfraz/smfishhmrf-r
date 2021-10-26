## Transitioning guide for HMRF

We recently updated the HMRF module in Giotto to make it easier to use and at the same time making it more powerful.

To take advantage of these changes, users should adjust their workflow according to the following:

### Change 1
As a new addition, after running binSpect/silhouetteRank we automatically record spatial gene information into the gene metadata table in the Giotto object. 
As a result, the return should be gobject instead of a data.table, and return_gobject should be set to TRUE:
```R
# old way
kmtest = binSpect(visium_brain, calc_hub = T, hub_min_int = 5,spatial_network_name = 'spatial_network')
```
should be changed to:
```R
# new way
visium_brain = binSpect(visium_brain, calc_hub = T, hub_min_int = 5,spatial_network_name = 'spatial_network', return_gobject=TRUE)
```

Spatial gene information can be accessed via:
```R
fDataDT(visium_brain)
# see new columns in the gene metadata table
names(fDataDT(visium_brain))
```
This change is required for running HMRF, since HMRF reviews the entire list of genes to set spatial gene filtering and sampling criteria. So it is essential to have spatial gene results in the Giotto object.

### Change 2

Usage of HMRF. HMRF now consists of two steps to make it more flexible. Step 1 is initialization (the initHMRF() function). Step 2 is HMRF (the doHMRF() function). 

#### Old way
In the old way, running HMRF requires several more steps of filtering spatial genes, spatial coexpression modules, and sampling spatial genes to determine a smaller subset, and finally doHMRF() function. Running these steps can require a lot of user input.
```R
#old way
# plot the spatial scores
plot(x=seq(1, 14414), y=-log10(spatial_genes$pval), xlab="Rank of genes by spatial score", ylab="-log10Pvalue")
abline(v=c(1500))

# take top 1500 genes to do spatial coexpression modules
ext_spatial_genes = spatial_genes[1:1500,]$gene
spat_cor_netw_DT = detectSpatialCorGenes(visium_brain, method = 'network', spatial_network_name = 'spatial_network', subset_genes = ext_spatial_genes, network_smoothing=0)
# cluster spatial genes
spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 20)

# sampling spatial genes
sample_rate=2; target=500; tot=0; num_cluster=20
gene_list = list()
clust = spat_cor_netw_DT$cor_clusters$spat_netw_clus
for(i in seq(1, num_cluster)){
    gene_list[[i]] = colnames(t(clust[which(clust==i)]))
}
for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    tot = tot+num_g/(num_g^(1/sample_rate))
}
factor=target/tot; num_sample=c()
for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    num_sample[i] = round(num_g/(num_g^(1/sample_rate)) * factor)
}
set.seed(10); samples=list(); union_genes = c()
for(i in seq(1, num_cluster)){
    if(length(gene_list[[i]])<num_sample[i]){
        samples[[i]] = gene_list[[i]]
    }else{
        samples[[i]] = sample(gene_list[[i]], num_sample[i])
    }
    union_genes = union(union_genes, samples[[i]])
}
union_genes = unique(union_genes)

# do HMRF with different betas on 500 spatial genes
my_spatial_genes <- union_genes
hmrf_folder = fs::path("11_HMRF")
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
HMRF_spatial_genes = doHMRF(gobject = visium_brain, expression_values = 'scaled', spatial_genes = my_spatial_genes, k = 20,   spatial_network_name="spatial_network", betas = c(0, 10, 5),  output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k20_scaled'))
```

#### New way

In the new way, running HMRF is much simplified, as we package several of the preprocessing steps together. Namely, there are now two separate functions, initHMRF(), and doHMRF() should be run in sequence. The initHMRF() function will perform spatial gene filtering, spatial coexpression module determination, sampling spatial genes, Kmeans initialization. In the old version of the code, all of these require separate commands. Now they are automated into the initHMRF() function. The doHMRF() will then perform HMRF given the initialization.

```R
# new way
hmrf=initHMRF(visium_brain, expression_values="scaled",
     spatial_network_name="spatial_network",
     use_spatial_genes="silhouetteRank", 
     gene_samples = 500, gene_sampling_rate = 2, gene_sampling_from_top = 2500,
     k=10, nstart=100, filter_method="none")
res=doHMRF(hmrf, betas=c(0, 10, 5))
```

- `use_spatial_genes`: `silhouetteRank` or `binSpect`. This specifies which version of Giotto's spatial genes to use.
- `gene_samples`, `gene_sampling_rate`, `gene_sampling_from_top`: Gene sampling parameters. The gene sampling rate is 1-50. When `gene_sampling_rate` is 1, it will sample equal number of genes from different spatial coexpression gene modules. When `gene_sampling_rate` is large, it will sample a number that is proportional to spatial coexpression module size. 
- `filter_method`: `none`: filter spatial genes based on top N genes set by `gene_sampling_from_top` sorted by spatial gene pvalue. `elbow`: filter genes based on elbow point of -log10Pvalue vs gene rank plot.
- `user_gene_list`: user-specified gene list, alternative to using spatial genes determined by Giotto (default NULL).

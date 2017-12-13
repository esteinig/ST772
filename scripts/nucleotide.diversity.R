# This script calculates pairwise nucleotide diversity by categorical data (e.g. geographic regions)
# as outlined by Stucki et al. (doi:10.1038/ng.3704) and pairwise SNP distances; 
# it generates plots for Figure 2b
#       
#     snippy:           output from snippy-core (.tab)
#     data:             file with two columns: sample id, categorical data (.tab)
#     groups:           vector of categories to include in analysis
#     bootstrap:        bootstrap replicates
#     verbose:          print results
#     seed:             random number seed for bootstrap
#     colours:          vector of colours (length groups) or palette
#     background:       include ggplot background, otherwise blank
#     outputs:          file names for statistics, histogram, boxplot and diversity plot

if (!require("pacman")) install.packages("pacman", repos="http://cran.us.r-project.org")
pacman::p_load("ape", "cowplot", "ggplot2", "RColorBrewer")

# Nucleotide Diversity Function

nucleotide.diversity <- function(snippy, data, groups, bootstrap=1000, verbose=T, save=TRUE,
                                 seed=23061989, colours="Dark2", background=FALSE, outputs=NULL) {
  
  if (length(colours) == 1){
    colours <- brewer.pal(n=length(groups), name=colours)
  }
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  if (is.null(outputs)){
    outputs = c("diversity.statistics.txt", "pairwise.histogram.pdf", "pairwise.boxplot.pdf", "diversity.pdf")
  }
  
  
  
  snps = read.csv(snippy, header=T, sep="\t")
  snps = snps[, -which(names(snps) %in% c("LOCUS_TAG", "GENE", "PRODUCT", "EFFECT", "CHR", "POS", "X"))]
  
  alignment.length = dim(snps)[1]
  
  data.matrix = t(as.matrix(snps)) 
  
  meta = read.csv(data, header=F, sep="\t")
  names(meta) <- c("key", "value")
  assignments <- as.character(meta$value)
  names(assignments) <- as.character(meta$key)
  
  
  rownames(data.matrix) <- as.vector(sapply(rownames(data.matrix), function(x){ assignments[x] }))
  names(rownames(data.matrix)) <- NULL
  
  dist.dfs <- lapply(groups, function(group){
    
    # Convert data matrix sub-setted for group to DNAbin:
    dna.bin <- as.DNAbin(data.matrix[rownames(data.matrix) == group,])
    
    # Calculate SNP difference counts deleting pairwise sites with at least one missing:
    dna.dist <- as.vector(dist.dna(dna.bin, model="N", pairwise.deletion = F))
    
    # Convert to data frame:
    dna.dist.df <- data.frame(group=rep(group,length(dna.dist)), dist=dna.dist)
    
  })
  
  dist.df <- do.call("rbind", dist.dfs)
  
  dist.df$group <- factor(dist.df$group, levels = groups)
  
  dist.boxplot <- ggplot(dist.df, aes(x=group, y=dist, fill=group, lower = 0.25, middle = 0.5, upper = 0.75)) +
    geom_boxplot() + xlab("\nRegion") + ylab("Pairwise SNP Distance\n") + guides(fill=F) +
    scale_fill_manual(values=colours)
  
  dist.df$group <- factor(dist.df$group, levels = groups)
  
  dist.histogram <- ggplot(dist.df, aes(x=dist, fill=group)) + geom_histogram(binwidth=5) + facet_wrap(~group, ncol=1, scales="free_y") +
    xlab("\nPairwise SNP Distance") + ylab("") + guides(fill=F) +
    scale_fill_manual(values=colours)
  
  test.comb <- combn(groups, 2, simplify=F)
  
  wilcox.results <- lapply(test.comb, function(x){
    
    wilcox.test(x=dist.df$dist[dist.df$group == x[1]], y=dist.df$dist[dist.df$group == x[2]])
    
  })
  
  wilcox.names <- lapply(test.comb, function(x){
    paste0(x[1], " vs. ", x[2])
  })
  
  names(wilcox.results) <- wilcox.names
  
  kruskal.list <- lapply(groups, function(x){
    dist.df$dist[dist.df$group == x]
  })
  
  kruskal.result <- kruskal.test(kruskal.list)
  
  dunn.result <- dunn.test::dunn.test(kruskal.list, method="bonferroni", kw=T)
  
  nuc.df <- dist.df
  
  div.results <- lapply(groups, function(group){
    
    g.dist <- dist.dna(as.DNAbin(data.matrix[rownames(data.matrix) == group,], model="N"))/alignment.length
    g.mean <- mean(as.vector(g.dist))
    
    values <- c()
    for (i in 1:bootstrap){
      bs.matrix <- data.matrix[,sample(1:ncol(data.matrix), replace=T)]
      d.dist <- dist.dna(as.DNAbin(bs.matrix[rownames(bs.matrix) == group,], model="N"))/alignment.length
      d.mean <- mean(as.vector(d.dist))
      values[i] <- d.mean
    }
    
    ci <- quantile(values, probs=c(0.025, 0.975))
    results <- list("group" = group,
                    "true.mean" = g.mean,
                    "low" = ci[1],
                    "high" = ci[2])
    
    return(results)
    
  })
  
  div.mat <- do.call("rbind", div.results)
  
  div.df <- data.frame(group=as.character(div.mat[,1]),
                       mean=as.numeric(div.mat[,2]),
                       low=as.numeric(div.mat[,3]),
                       high=as.numeric(div.mat[,4]))
  
  div.df$group <- factor(div.df$group, levels = groups)
  
  div.plot <- ggplot(data=div.df, aes(x=group, y=mean)) + geom_point(colour=colours, size=3) + 
    geom_errorbar(aes(ymin=low, ymax=high, width=0.1), colour=colours) + theme_cowplot() +
    ylab(expression(paste("Nucleotide diversity ( ", pi," )"))) + xlab("") +
    theme(legend.title=element_blank()) #+ scale_x_discrete(limits=groups)
  
  if (background == FALSE){
    bck <- theme_cowplot()
    dist.histogram <- dist.histogram + bck
    dist.boxplot <- dist.boxplot + bck
    div.plot <- div.plot + bck
  }
  
  
  if (verbose) {
    print(dist.histogram)
    print(dist.boxplot)
    print(wilcox.results)
    print(kruskal.result)
    print(div.plot)
  }
  
  if (save) {
    
    print(outputs)
    
    sink(outputs[1])
    print(wilcox.results)
    print(kruskal.result)
    print(dunn.result)
    sink()
    
    ggsave(outputs[2], plot=dist.histogram, dpi=600)
    ggsave(outputs[3], plot=dist.boxplot, dpi=600)
    ggsave(outputs[4], plot=div.plot, dpi=600)
    
  }
  
  final <- list("pairwise.histogram" = dist.histogram,
                "pairwise.boxplot" = dist.boxplot,
                "pairwise.data" = dist.df,
                "diversity.plot" = div.plot,
                "diversity.data" = div.df,
                "pairwise.wilcox" = wilcox.results,
                "pairwise.kruskal" = kruskal.result)
  
  return(final)
  
}

# The function uses a renamed version of the Snippy core alignment due to reading sample
# names starting with numbers (10V3510, 11ANR_12420) incorrectly. These have been renamed to: 
#
# SA_10V3510 and SA_11ANR_12420
#
# in both input files: snippy.core.renamed.tab and regions.tab

#nucleotide.diversity(snippy="diversity.snippy.tab", data="diversity.regions.tab", groups=c("Australasia", "Europe", "South Asia"),
#                     colours=c("#fd8d3c", "#4292c6", "#810f7c"), bootstrap=1000, seed=23061989, save=TRUE)
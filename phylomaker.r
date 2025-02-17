library(V.PhyloMaker2)
library(phytools)
library(geiger)
library(tidyverse)
spp <- read.csv("sample_species_list_asteraceae.csv")
tree <- phylo.maker(sp.list = spp, scenarios="S3")
#write.tree(tree$scenario.3, "sample.tre")
plot(tree$scenario.3, type="fan",lwd=70)


par1 = par()
#pdf('tree1.pdf')
par(mar = c(0.1,0.1,0.1,0.1))
plot(tree$scenario.3,lwd = 2, cex = c(0.65,0.65), type="fan")
dev.off()

write.tree(tree, 'sample_tree2.tre')

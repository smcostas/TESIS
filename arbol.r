library(tidyverse);library(phytools);library(corHMM);library(geiger)
library(viridis); library(nlme); library(ape)

## tree #####

tree <- read.tree('sample_tree2.tre')
plot(tree)
## DATA #########
df_sum <- read.csv('final_df.csv', row.names = 1)

## check tree data ####
nchk = name.check(tree, df_sum)

nchk$data_not_tree

nchk$tree_not_data

fake = nchk$tree_not_data

true = c('Angelphytum', 'Schkuria_pinnata', 'Chaptalia', 'Chaptalia_3', 'Lorentzianthus_viscidus', 'Sanvitalia_sp', 'Flourensia_oolepis',
         'Flourensia_riparia', 'Kauna_saltensis', 'Heliantheae_monte', 'Burkartia_lanigera', 'Verbesina_encelioides', 'Austrobrickellia_patens') ### estan en orden necesario ###


for (i in seq_along(fake)) {
  # Busca el índice del nombre incorrecto en tree$tip.label
  idx <- which(tree$tip.label == fake[i])  
  
  if (length(idx) > 0) {
    tree$tip.label[idx] <- true[i]
  }
}

plot(tree)

nchk = name.check(tree, df_sum)


df_sum = df_sum[tree$tip.label,] ## ya ordenado en el sentido del arbol ## necesario
write.tree(tree,'final_tree.tre')
pdf('final_tree.pdf', width = 10, height = 10)
png('final_tree.png', width = 21, height = 20, units = 'cm', res = 600)
par(mar = c(0.1,0.1,0.1,0.1))
plot(tree,lwd = 2, cex = c(0.65,0.65), type="fan")
dev.off()



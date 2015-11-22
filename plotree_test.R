library(ape)

tree <- read.tree("ph_dmc.nre")
tree2 <- read.tree("rev_dmc.nre")  # rev_dmc and ph1_DNA_F can be read now, although they are just 2 nodes.
tree3 <- read.tree("ph1_DNA_F.nre")

plot(tree)
plot(tree2)
plot(tree3)

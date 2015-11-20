library(ape)

tree <- read.tree("ph_dmc.nre")
# tree2 <- read.tree("rev_dmc.nre")  # this tree still just has "singleton nodes" and can't be read

plot(tree)
# plot(tree2)

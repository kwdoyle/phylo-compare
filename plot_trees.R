library(ape)

#labels <- c("Taxon1", "Taxon2", "Taxon3", "Taxon4", "Taxon5", "Taxon6",
#            "Taxon7", "Taxon8", "Taxon9", "Taxon10")

# make taxon labels just be their numbers instead
labels <- seq(1,10)

# 50 nodes
T_H_50 <- read.tree("ph50.nre")
T_DMC_50 <- read.tree("ph50_dmc.nre")
T_seq_50 <- read.tree("ph50_seq.nre")
T_F_50 <- read.tree("ph50_DNA_F.nre")

plot(T_H_50, main="True Tree?", cex=0.6)
plot(T_DMC_50, main="DMC Tree", cex=0.6)
plot(T_seq_50, main="Seq Tree", cex=0.6)
plot(T_F_50, main="(DMC + Seq) Tree", cex=0.6)

### 50 nodes is very messy to look at, and the symmetric difference proportions between
# each of the trees are the same for 50 nodes and 10 nodes, so it is probably better
# to use 10 nodes instead.

# 10 nodes
T_H_10 <- read.tree("ph10.nre")
T_DMC_10 <- read.tree("ph10_dmc.nre")
T_seq_10 <- read.tree("ph10_seq.nre")
T_F_10 <- read.tree("ph10_DNA_F.nre")

# change taxon names to just be numbers on the trees
T_H_10$tip.label <- sub("Taxon", "", T_H_10$tip.label, fixed = T)
T_DMC_10$tip.label <- sub("Taxon", "", T_DMC_10$tip.label, fixed = T)
T_seq_10$tip.label <- sub("Taxon", "", T_seq_10$tip.label, fixed = T)
T_F_10$tip.label <- sub("Taxon", "", T_F_10$tip.label, fixed = T)

# plot
plot(T_H_10, main="True Tree?")
plot(T_DMC_10, main="DMC Tree")
plot(T_seq_10, main="Seq Tree")
plot(T_F_10, main="(DMC + Seq) Tree")

# make matrix of labels to be used in cophyloplot
#mat <- matrix(data=c(labels, labels), ncol=2)
mat <- matrix(c(labels, labels), ncol=2)

# True tree(?) with sequence tree
# looks as if the sequence tree gives the same relationships as in the true relationship tree
cophyloplot(T_H_10, T_seq_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# True tree(?) with DMC tree
cophyloplot(T_H_10, T_DMC_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# True tree(?) with (sequence + DMC) tree
cophyloplot(T_H_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# DMC tree with sequence tree
cophyloplot(T_DMC_10, T_seq_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# DMC tree with (sequence + DMC) tree
cophyloplot(T_DMC_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# sequence tree with (sequence + DMC) tree
cophyloplot(T_seq_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3)




# calculates symmetric difference; same as in calc_distance.py
dist.topo(T_H_10, T_DMC_10)
dist.topo(T_DMC_10, T_seq_10)

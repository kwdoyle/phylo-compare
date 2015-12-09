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

T_F_10_diff_qmod_and_qcon <- read.tree("ph10_2.1_DNA_F.nre")

# 10 nodes, alpha = 1
T_H_alpha1 <- read.tree("alpha-1.nre")
T_F_alpha1 <- read.tree("alpha-1_DNA_F.nre")
rev_dmc_alpha1 <- read.tree("rev_dmc_alpha-1.nre")

# change taxon names to just be numbers on the trees
T_H_10$tip.label <- sub("Taxon", "", T_H_10$tip.label, fixed = T)
T_DMC_10$tip.label <- sub("Taxon", "", T_DMC_10$tip.label, fixed = T)
T_seq_10$tip.label <- sub("Taxon", "", T_seq_10$tip.label, fixed = T)
T_F_10$tip.label <- sub("Taxon", "", T_F_10$tip.label, fixed = T)
T_F_10_diff_qmod_and_qcon$tip.label <- sub("Taxon", "", T_F_10$tip.label, fixed = T)
T_H_alpha1$tip.label <- sub("Taxon", "", T_H_alpha1$tip.label, fixed = T)
T_F_alpha1$tip.label <- sub("Taxon", "", T_F_alpha1$tip.label, fixed = T)
rev_dmc_alpha1$tip.label <- sub("Taxon", "", rev_dmc_alpha1$tip.label, fixed = T)

# plot
plot(T_H_10, main="True Tree")
plot(T_DMC_10, main="DMC Tree")
plot(T_seq_10, main="Seq Tree")
plot(T_F_10, main="(DMC + Seq) Tree")
#nodelabels()  # this will put labels on the nodes of the tree. they're kind of big and useless though
plot(T_F_10_diff_qmod_and_qcon, main="(DMC + Seq) Tree", sub="qmod=0.9, qcon=0.1")

# make matrix of labels to be used in cophyloplot
#mat <- matrix(data=c(labels, labels), ncol=2)
mat <- matrix(c(labels, labels), ncol=2)

# True tree(?) with sequence tree
# looks as if the sequence tree gives the same relationships as in the true relationship tree
cophyloplot(T_H_10, T_seq_10, assoc=mat, length.line = 4, space = 48, gap = 3, rotate=T)

# True tree(?) with DMC tree
cophyloplot(T_H_10, T_DMC_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# True tree(?) with (sequence + DMC) tree
cophyloplot(T_H_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# DMC tree with sequence tree
cophyloplot(T_DMC_10, T_seq_10, assoc=mat, length.line = 4, space = 48, gap = 3, rotate=T)

# DMC tree with (sequence + DMC) tree
cophyloplot(T_DMC_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3)

# sequence tree with (sequence + DMC) tree
cophyloplot(T_seq_10, T_F_10, assoc=mat, length.line = 4, space = 48, gap = 3, rotate=T)

# Rotating these trees reveals that they all have the same relationships..
# the branch lengths are different between trees though.. I guess that's what's important?

# sequence tree with (sequence + DMC) tree which used qmod=0.9 and qcon=0.1
cophyloplot(T_seq_10, T_F_10_diff_qmod_and_qcon, assoc=mat, length.line = 4, space = 48, gap = 3, rotate=T)


# true tree and combo tree with alpha = 1
cophyloplot(T_H_alpha1, T_F_alpha1, assoc=mat, length.line = 4, space = 48, gap = 3)

# true tree and rev_dmc tree w/ alpha=1?
cophyloplot(T_H_alpha1, rev_dmc_alpha1, assoc=mat, length.line = 4, space = 48, gap = 3)







# calculates symmetric difference; same as in calc_distance.py
dist.topo(T_H_10, T_DMC_10)
dist.topo(T_DMC_10, T_seq_10)



# These are the branch lengths
TH_branch <- T_H_10$edge.length
DMC_branch <- T_DMC_10$edge.length
seq_branch <- T_seq_10$edge.length
combo_branch <- T_F_10$edge.length

branch_lengths_mat <- matrix(c(TH_branch, DMC_branch, seq_branch, combo_branch), nrow=length(TH_branch), ncol=4, byrow=F)
colnames(branch_lengths_mat) <- c("T_H", "T_DMC", "T_seq", "T_comb")
branch_lengths <- data.frame(branch_lengths_mat)


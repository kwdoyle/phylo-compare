library(ape)

T_H <- read.tree("ph50.nre")
T_DMC <- read.tree("ph_dmc.nre")
T_seq <- read.tree("ph_seq.nre")
T_F <- read.tree("ph50_DNA_F.nre")

plot(T_H, main="True Tree?", cex=0.6)
plot(T_DMC, main="DMC Tree", cex=0.6)
plot(T_seq, main="Seq Tree", cex=0.6)
plot(T_F, main="(DMC + Seq) Tree", cex=0.6)
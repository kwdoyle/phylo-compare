import dendropy
from dendropy.calculate import treecompare

T_H = dendropy.Tree.get(file=open("ph6.nre", "r"), schema="newick")
T_DMC = dendropy.Tree.get(file=open("ph_dmc.nre", "r"), schema="newick")
T_seq = dendropy.Tree.get(file=open("ph_seq.nre", "r"), schema="newick")
T_F = dendropy.Tree.get(file=open("ph6_DNA_F.nre", "r"), schema="newick")


print "Distance between T_H with T_DMC " + str(treecompare.weighted_robinson_foulds_distance(T_H, T_DMC))#T_H.robinson_foulds_distance(T_DMC))
print "Distance between T_H with T_seq " + str(T_H.robinson_foulds_distance(T_seq))
print "Distance between T_H with T_F " + str(T_H.robinson_foulds_distance(T_F))

# Calculating these distances gives an error about a non-identical taxon namespace between the two trees.
# This most likely means that the taxons assigned to each tree literally need to be the same object.

print "Distance between T_H with T_DMC " + str(T_H.symmetric_difference(T_DMC))
print "Distance between T_H with T_seq " + str(T_H.symmetric_difference(T_seq))
print "Distance between T_H with T_F " + str(T_H.symmetric_difference(T_F))

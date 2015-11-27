import dendropy
from dendropy.calculate import treecompare

# Using tree with largest possible amount of nodes (6)
T_H = dendropy.Tree.get(file=open("ph50.nre", "r"), schema="newick")
T_DMC = dendropy.Tree.get(file=open("ph_dmc.nre", "r"), schema="newick")
T_seq = dendropy.Tree.get(file=open("ph_seq.nre", "r"), schema="newick")
T_F = dendropy.Tree.get(file=open("ph50_DNA_F.nre", "r"), schema="newick")

# Turn tree files into strings
T_H_string = str(T_H) + ";"
T_DMC_string = str(T_DMC) + ";"
T_seq_string = str(T_seq) + ";"
T_F_string = str(T_F) + ";"

# Make taxon namespaces using those strings
T_H_list = dendropy.TreeList()
T_H_list.read(data=T_H_string, schema="newick")

T_DMC_list = dendropy.TreeList(taxon_namespace=T_H_list.taxon_namespace)
T_DMC_list.read(data=T_DMC_string, schema="newick")

T_seq_list = dendropy.TreeList(taxon_namespace=T_H_list.taxon_namespace)
T_seq_list.read(data=T_seq_string, schema="newick")

T_F_list = dendropy.TreeList(taxon_namespace=T_H_list.taxon_namespace)
T_F_list.read(data=T_F_string, schema="newick")


# so this creation of the namespaces works, as these two lists below will print out the contents of
# the original tree files, AND the statistical functions below also successfully run.
#print T_H_list[0]
#print T_DMC_list[0]

# Calculating symmetric differences
# symmetric difference is the number of splits found in one of the trees but not the other
print "Symmetric difference between T_H and T_DMC: " + str(treecompare.symmetric_difference(T_H_list[0], T_DMC_list[0]))
print "Symmetric difference between T_H and T_seq: " + str(treecompare.symmetric_difference(T_H_list[0], T_seq_list[0]))
print "Symmetric difference between T_H with T_F: " + str(treecompare.symmetric_difference(T_H_list[0], T_F_list[0]))


# Calculating the robinson foulds distances
# this is the weighted symmetric difference, which is the sum of the square of differences in branch lengths for equivalent splits between two trees
print "Robinson-Foulds distance between T_H and T_DMC: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_DMC_list[0]))
print "Robinson-Foulds distance between T_H and T_seq: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_seq_list[0]))
print "Robinson-Foulds distance between T_H and T_F: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_F_list[0]))

import sys
import dendropy
from dendropy.calculate import treecompare

# Use files to analyze as command line arguments
try:
    ph_tree = sys.argv[1]
    dmc_tree = sys.argv[2]
    seq_tree = sys.argv[3]
    combo_tree = sys.argv[4]
except:
    print "Usage: python calc_distance.py phX.nre ph_dmc.nre ph_seq.nre phX_DNA_F.nre"
    exit()

T_H = dendropy.Tree.get(file=open(ph_tree, "r"), schema="newick") # "ph50.nre"
T_DMC = dendropy.Tree.get(file=open(dmc_tree, "r"), schema="newick") # "ph_dmc.nre"
T_seq = dendropy.Tree.get(file=open(seq_tree, "r"), schema="newick") # "ph_seq.nre"
T_F = dendropy.Tree.get(file=open(combo_tree, "r"), schema="newick") # "ph50_DNA_F.nre"

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

# Calculating symmetric differences (unweighted robinson foulds).
# symmetric difference is the number of splits found in one of the trees but not the other.
# it is defined as the number of transformations needed to turn one tree into the other.
print "Symmetric difference between T_H and T_DMC: " + str(treecompare.symmetric_difference(T_H_list[0], T_DMC_list[0]))
print "Symmetric difference between T_H and T_seq: " + str(treecompare.symmetric_difference(T_H_list[0], T_seq_list[0]))
print "Symmetric difference between T_H with T_F: " + str(treecompare.symmetric_difference(T_H_list[0], T_F_list[0]))
print "Symmetric difference between T_DMC with T_seq: " + str(treecompare.symmetric_difference(T_DMC_list[0], T_seq_list[0]))
print "Symmetric difference between T_DMC with T_F: " + str(treecompare.symmetric_difference(T_DMC_list[0], T_F_list[0]))
print "Symmetric difference between T_seq with T_F: " + str(treecompare.symmetric_difference(T_seq_list[0], T_F_list[0]))

# Calculating the robinson foulds distances
# This is the weighted symmetric difference, which is the sum of the square of differences in branch lengths for equivalent splits between two trees.
# It takes edge lengths into account, and therefore will yield a non-zero answer for trees with identical relationships, but have different branch lengths.
# This explains why the unweighted distance between T_H and T_seq is 0, but is >0 for the weighted distance.
print "Robinson-Foulds distance between T_H and T_DMC: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_DMC_list[0]))
print "Robinson-Foulds distance between T_H and T_seq: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_seq_list[0]))
print "Robinson-Foulds distance between T_H and T_F: " + str(treecompare.weighted_robinson_foulds_distance(T_H_list[0], T_F_list[0]))
print "Robinson-Foulds distance between T_DMC and T_seq: " + str(treecompare.weighted_robinson_foulds_distance(T_DMC_list[0], T_seq_list[0]))
print "Robinson-Foulds distance between T_DMC and T_F: " + str(treecompare.weighted_robinson_foulds_distance(T_DMC_list[0], T_F_list[0]))
print "Robinson-Foulds distance between T_seq and T_F: " + str(treecompare.weighted_robinson_foulds_distance(T_seq_list[0], T_F_list[0]))

### Note: running the framework twice for the same nodes, edges, qmod, etc, will give the same exact trees each time
# Maybe we can try changing qmod and qcon? Does it even matter?

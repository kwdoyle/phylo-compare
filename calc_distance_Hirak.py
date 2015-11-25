import dendropy
from dendropy.calculate import treecompare

T_H = dendropy.Tree.get(file=open("ph6.nre", "r"), schema="newick")
T_DMC = dendropy.Tree.get(file=open("ph_dmc.nre", "r"), schema="newick")
T_seq = dendropy.Tree.get(file=open("ph_seq.nre", "r"), schema="newick")
T_F = dendropy.Tree.get(file=open("ph6_DNA_F.nre", "r"), schema="newick")

# Make taxon namespaces
T_DMC_string = str(T_DMC) + ";"
T_seq_string = str(T_seq) + ";"

# testing
tree_list1 = dendropy.TreeList()
tree_list1.read(data=T_DMC_string, schema="newick")
tree_list2 = dendropy.TreeList(taxon_namespace=tree_list1.taxon_namespace)
tree_list2.read(data=T_seq_string, schema="newick")

# Results in: 6
print(treecompare.symmetric_difference(tree_list1[0], tree_list2[0]))

# so this creation of the namespaces works, as these two tree_lists below will print out the contents of
# the original tree files, AND the symmetric_difference function above also successfully runs.
#print tree_list1[0]
#print tree_list2[0]

# Calculating the robinson foulds distance will now work with this method too!!!
print(treecompare.weighted_robinson_foulds_distance(tree_list1[0], tree_list2[0]))


#print "Distance between T_H with T_DMC " + str(treecompare.weighted_robinson_foulds_distance(T_H, T_DMC))#T_H.robinson_foulds_distance(T_DMC))
#print "Distance between T_H with T_seq " + str(T_H.robinson_foulds_distance(T_seq))
#print "Distance between T_H with T_F " + str(T_H.robinson_foulds_distance(T_F))

# Calculating these distances gives an error about a non-identical taxon namespace between the two trees.
# This most likely means that the taxons assigned to each tree literally need to be the same object.

#print "Distance between T_H with T_DMC " + str(T_H.symmetric_difference(T_DMC))
#print "Distance between T_H with T_seq " + str(T_H.symmetric_difference(T_seq))
#print "Distance between T_H with T_F " + str(T_H.symmetric_difference(T_F))

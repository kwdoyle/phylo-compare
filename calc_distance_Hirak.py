import dendropy

T_H = dendropy.Tree(stream=open("ph3.nre"), schema="newick")
T_DMC = dendropy.Tree(stream=open("ph_dmc.nre"), schema="newick")
T_seq = dendropy.Tree(stream=open("ph_seq.nre"), schema="newick")
T_F = dendropy.Tree(stream=open("ph3_DNA_F.nre"), schema="newick")


print "Distance between T_H with T_DMC " + str(T_H.robinson_foulds_distance(T_DMC))
print "Distance between T_H with T_seq " + str(T_H.robinson_foulds_distance(T_seq))
print "Distance between T_H with T_F " + str(T_H.robinson_foulds_distance(T_F))


print "Distance between T_H with T_DMC " + str(T_H.symmetric_difference(T_DMC))
print "Distance between T_H with T_seq " + str(T_H.symmetric_difference(T_seq))
print "Distance between T_H with T_F " + str(T_H.symmetric_difference(T_F))
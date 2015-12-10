from __future__ import division
from optparse import OptionParser
from collections import deque
from math import log
import sys,time,logging,random
import os
import dendropy
import StringIO
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
from Bio import AlignIO
import copy

random.seed(10301949)

H = []

#==============================================================================
#                               NODE OBJECT
#==============================================================================
class Node:
    def __init__(self,id):
        self.id = id
        self.nodes = [id]
        self.neighbors = set()

    def __cmp__(self,other):
        if self.id > other.id: return 1
        elif self.id < other.id: return -1
        else: return 0

    def __str__(self):
        s =  "%s --> " %(self.id)
        for index,neighbor in enumerate(self.neighbors):
            s += "(%i) " %(neighbor) # self.weights[index]
        s += "\n"
        return s


#==============================================================================
#                               GRAPH OBJECT
#==============================================================================
class Graph:
    def __init__(self,name):
        self.node_map = {} # Mapping from node id#'s to Node objects.
        self.num_nodes = 0
        self.num_edges = 0
        self.name = name

    def __iter__(self):
        return iter(self.node_map)

    def __str__(self):
        s = "%s: #nodes: %s, #edges: %s\n\n" %(self.name,self.num_nodes,self.num_edges)
        for i in self.node_map.iterkeys():
            s += "%s" %(self.node_map[i])
        return s

    @staticmethod
    def read_graph(filename):
        """ Creates a graph object from an edgelist file. """
        G = Graph("")
        input = open(filename)
        for line in input:
            [u,v] = line.strip().split()
            if u == v: continue
            if not G.has_node(u): G.add_node(u)
            if not G.has_node(v): G.add_node(v)
            G.add_edge(u,v)

        logging.info("Graph: #nodes=%i, #edges=%i" %(G.num_nodes,G.num_edges))
        return G
    @staticmethod
    def rand_graph(nodes,edges,outfile):
        """ Make a random graph according to arguments """
        G = Graph("")
        target = open(outfile,'w')

        for i in range(edges):
            u = random.randint(1,nodes)
            v = random.randint(1,nodes)
            if u == v: continue
            if not G.has_node(u): G.add_node(u)
            if not G.has_node(v): G.add_node(v)
            if not (G.has_edge(u,v) or (G.has_edge(v,u))):
                G.add_edge(u,v)
                line = str(u) + "\t" + str(v)
                target.write(line)
                target.write("\n")

        target.close()
        return G


    def add_node(self,i):
        """ Adds a given node to the graph, as an isolate. """
        self.num_nodes += 1
        self.node_map[i] = Node(i)

    def remove_node(self,i):
        """ Deletes the given node from the graph. """
        I = self.node_map[i]

        # Remove i from the neighbor lists of all of i's neighbors.
        for neighbor in I.neighbors:
            if i != neighbor:
                self.node_map[neighbor].neighbors.remove(i)
            self.num_edges -= 1

        del self.node_map[i]
        self.num_nodes -= 1

    def has_node(self,u):
        """ Returns True if the given node is in the graph; False otherwise. """
        return u in self.node_map

    def add_edge(self,u,v):
        """ Adds edge to the graph. """
        self.node_map[u].neighbors.add(v)
        if u != v: self.node_map[v].neighbors.add(u)
        self.num_edges += 1

    def delete_edge(self,u,v):
        """ Deletes edge from the graph. """
        self.node_map[u].neighbors.remove(v)
        if u != v: self.node_map[v].neighbors.remove(u)
        self.num_edges -= 1

    def has_edge(self,u,v):
        """ Returns True if the edge exists in the graph; False otherwise. """
        return v in self.node_map[u].neighbors

    def nodes(self):
        """ Return all keys in the graph object. """
        return self.node_map.keys()

    def random_node(self):
        """ Returns a random node id. """
        # TODO: something better than populating keys() each time?
        return random.choice(self.node_map.keys())

    def neighbors(self,u):
        """ Returns set of neighbors of u. """
        return self.node_map[u].neighbors

    def dmc_likelihood(self,u,v,qmod,qcon):
        """ Computes the likelihood of u merging with v according to the DMC model. """
        # Ignore edges to each other because those will be taken care of by qcon.
        U = self.node_map[u].neighbors - set([v])
        V = self.node_map[v].neighbors - set([u])

        Cuv = len(U & V)
        Suv = len(U) + len(V) - 2*Cuv
        gamma = qcon if self.has_edge(u,v) else 1-qcon

        return (1-qmod)**(Cuv) * (qmod/2)**(Suv) * gamma


#==============================================================================
# We start with two simple node and generate a graph. Each step we randomly
# select a node that would be expanded. So a node $u$ will spawn of another
# node say $v$. Now in this way we will go on adding edges. The following
# function has the following functionality
# Input: number of nodes, number of edges, q_{mod} and q_{con}
# Output: a graph G, during the process generated phylogeny tree
#==============================================================================

def forward_dmc(nodes,edges,outfile,qmod,qcon,treefile):
    """ Make a DMC random graph according to arguments """
    G = Graph("")
    target = open(outfile,'w')
    anchor_node_dic = {}


    new_nodes = set(range(1,nodes+1))


    # Add two nodes in the graph
    # They can have an edge according to
    # q_con probabbility
    print len(new_nodes)


    while True:
        print 'in loop'
        u = random.randint(1,nodes+1)
        v = random.randint(1,nodes+1)
        if not G.has_node(u):
            G.add_node(u)
            new_nodes.remove(u)
        if not G.has_node(v):
            G.add_node(v)
            new_nodes.remove(v)  # is this causing the KeyError when 7 nodes is attempted? I always get KeyError: 3, so maybe somehow node 3 is removed, and then is checked for again but now isn't there?
        if u != v:
            if random.random() < qcon:  # this isn't actually adding the distances from each node, which I think is what qcon and qmod should essentially be...
                G.add_edge(u,v)
                # Create newick tree
                # tree set-up should be like:
                # ((Taxon[number] : [distance from node], other Taxon[number] : [distance from node])
                # [useless number?] : [distance of this clade from node],
                # repeat for other taxons);

                # yes, that useless long string of numbers that gets added to the tree file is messing it up
                # indel-seq-gen is able to read the tree file if that sequence of numbers isn't there

                ### the useless string of numbers is actually the object ID obtained from id()

                string_root = "(" + "Taxon" + str(u) + ":" + str(1) + "," + "Taxon" + str(v) + ":" + str(1) + ")" + ";"
                tree = dendropy.Tree.get(data=string_root,schema="newick")
            break

    print len(new_nodes)


    while len(new_nodes) > 0:
        # Chooding anchor node
        new_parent = G.random_node()
        #print new_parent

        # Create new node
        child_node = random.sample(new_nodes,1)
        child_node = child_node[0]
        G.add_node(child_node)
        anchor_node_dic[child_node] = new_parent
        new_nodes.remove(child_node)
        # Create newick tree
        string_temp = "(" + "Taxon" + str(new_parent) + ":" + str(0.1) + "," + "Taxon" + str(child_node) + ":" + str(0.1) + ")" + ";"
        temp_tree = dendropy.Tree.get(data=string_temp,schema="newick")
        current_parent = filter(lambda x: x.taxon.label == "Taxon" + str(new_parent), [y for y in tree.leaf_nodes()])
        temp_tree_nodes = [x for x in temp_tree.nodes()]
        current_parent[0].add_child(temp_tree_nodes[1])
        current_parent[0].add_child(temp_tree_nodes[2])
        current_parent[0].taxon.label = id(current_parent[0]) #.oid  # used to be .oid at the end of current_parent[0]. Removing oid causes a maximum recursion depth error
        # I think we need the pygit2 module to use .oid, which is a git attribute
        # but this actually makes sense considering git uses a tree structure to handle all
        # commits and such

        # maybe the function id() can be used instead of .oid
        # nope, that didn't fix the namespace error.
        # not having id() or .oid there doesn't give a max recursion overflow error anymore, but it still doesn't fix the namespace error.

        # attemping to fix the .oid identifier deal by just assigning each label as the current_parent in string form

        # Add nodes according to q_mod
        #
        for x in G.neighbors(new_parent):
            G.add_edge(x,child_node)


        # shallow copy of neighbors
        adjacent_nodes = copy.copy(G.neighbors(child_node))

        for x in adjacent_nodes:
            if random.random() < qmod:
                if random.random() < 0.5:  # should these 0.5s be qmod and qcon instead?
                    if G.has_edge(x,child_node):  # or maybe not, because this is for the forward dmc step.
                        G.delete_edge(x,child_node)
                elif (G.has_edge(x,new_parent)):
                        G.delete_edge(x,new_parent)
        if random.random() < 0.5:
            if not G.has_edge(new_parent,child_node):
                G.add_edge(new_parent,child_node)

        print len(new_nodes)


    con_nodes = set()

    for u in G:
        for v in G:
            if u >= v: continue
            elif (G.has_edge(u,v) or G.has_edge(v,u)):
                con_nodes.add(u)
                con_nodes.add(v)
                line = str(u) + "\t" + str(v)
                target.write(line)
                target.write("\n")

    """ Without this part the resulting graph was disconnected"""
    """ This step makes sure that no node remains isolated """

    for u in G:
        if not u in con_nodes:
            G.add_edge(u,anchor_node_dic[u])
            line = str(u) + "\t" + str(anchor_node_dic[u])
            target.write(line)
            target.write("\n")



    target.close()


    out = tree.as_string('newick')
    out =  out.replace('[&U]','[50]') ## 50 unit long sequences
    open(treefile,'w').write(out)    # The first T produced by the history
    print(tree.as_ascii_plot())



    return G

#==============================================================================
# This function creates a tree made in reverse direction. In each step
# it merges two nodes that has maximum likelihood to be emerged from the
# same
# In this version of framework we introduce a sampling while merging
# that is we would merge two nodes with probability proportional to
# their likelihood.
#==============================================================================

def dmc_delorean(G,qmod,qcon):
    """ Reconstructs the network using the dmc model and delorean algorithm. """

    # Make initial pairwise likelihoods.
    #target = open(hisfile,'w')

    global H
    L = {}
    for u in G: L[u] = {}

    for u in G:
        for v in G:
            if u >= v: continue
            L[u][v] = G.dmc_likelihood(u,v,qmod,qcon)
            print L[u][v]
        print '\n'
    level_counter = 0

    while (G.num_nodes >= 2):
    	# at least two nodes in the graph.
		# Get largest Luv.
        L_list = []
        L_prob = -10000000000


        norm_sum = 0
        num_pairs = 0
        for u in G:
        	for v in G:
        		if u >= v: continue
        		norm_sum = norm_sum + L[u][v]


        L_norm = {}
    	for u in G: L_norm[u] = {}

    	for u in G:
        	for v in G:
        		if u >= v: continue
        		L_norm[u][v] = L[u][v] / norm_sum

        # Finding cumulative probability
        L_cum = {}
        L_cum[num_pairs] = 0
        L_index = {}
        for u in G:
        	for v in G:
        		if u >= v: continue
        		num_pairs = num_pairs + 1
        		L_cum[num_pairs] = L_cum[num_pairs-1] + L_norm[u][v]
        		L_index[num_pairs] = (u,v)

        # Binary search

        random_number = random.random()
        start_index = 0
        end_index = num_pairs

        while start_index != end_index - 1:
        	pivot = int((start_index + end_index)/2)
        	if random_number >= L_cum[pivot] and random_number <= L_cum[pivot+1]:
        		break
        	elif random_number < L_cum[pivot]:
        		end_index = pivot
        	elif random_number > L_cum[pivot]:
        		start_index = pivot





        '''

        for u in G:
            for v in G:
                if u >= v: continue

                Luv = L[u][v]
                if Luv > L_prob:
                    L_list = [(u,v)]
                    L_prob = Luv
                elif Luv == L_prob:
                    L_list.append((u,v))


        '''
        u,v = L_index[end_index]

        # Choose random pair; assign random daddy.
        #pair = random.choice(L_list)
        #(u,v) = (pair[0],pair[1]) if random.random() > 0.5 else (pair[1],pair[0])

        # Nodes whose likelihood values need to be computed.
        altered = (G.neighbors(u) | G.neighbors(v) | set([u])) - set([v])

        # Prepare to delete v: add new edges in symmetric difference of v to u.
        for neighbor in G.neighbors(v):
            if u == neighbor: continue # Don't add self-edge.
            elif v == neighbor: continue # Don't add, will remove v anyways.
            elif G.has_edge(u,neighbor): continue # Edge already exists.
            else: G.add_edge(u,neighbor)
        G.remove_node(v)

        H.append((u,v))
        print "%s\t%s" %(u,v)

        # Fix up altered Luv values.
        for x in altered:
            for y in G:
                if x == y: continue
                L[min(x,y)][max(x,y)] = G.dmc_likelihood(x,y,qmod,qcon)

    last_node = G.nodes()[0]
    H.append((last_node,last_node))
    print "%s\t%s" %(last_node,last_node)





#==============================================================================
# The tree file H from the previous function is used here to create
# the nre file. Dendropy to create the phylogeny Tree T from reverse
# DMC algorithms
#==============================================================================
def generate_tree_dendro(H,outfile):
    length = len(H)
    length = length - 2
    '''

    tree = dendropy.Tree(stream=StringIO.StringIO(str(H[length])),schema="newick")
    for i in range(length-1,-1,-1):
        (u,v) = H[i]
        temp_tree = dendropy.Tree(stream=StringIO.StringIO(str(H[i])),schema="newick")
        current_parent = filter(lambda x: x.taxon.label == str(u), [y for y in tree.leaf_nodes()])
        temp_tree_nodes = [x for x in temp_tree.nodes()]
        current_parent[0].add_child(temp_tree_nodes[1])
        current_parent[0].add_child(temp_tree_nodes[2])


    print(tree.as_ascii_plot())
    '''
    (u,v) = H[length]
    counter = 1
    #[999]((Taxon1:0.3,Taxon2:0.14):0.5,(Taxon3:0.34, Taxon4:0.5):0.12);
    string_root = "(" + "Taxon" + str(u) + ":" + str(1) + "," + "Taxon" + str(v) + ":" + str(1) + ")" + ";"
    tree = dendropy.Tree.get(data=string_root, schema="newick")
    for i in range(length-1,-1,-1):
        (u,v) = H[i]
        string_temp = "(" + "Taxon" + str(u) + ":" + str(1) + "," + "Taxon" + str(v) + ":" + str(1) + ")" + ";"
        temp_tree = dendropy.Tree.get(data=string_temp, schema="newick")
        current_parent = filter(lambda x: x.taxon.label == "Taxon" + str(u), [y for y in tree.leaf_nodes()])
        temp_tree_nodes = [x for x in temp_tree.nodes()]
        current_parent[0].add_child(temp_tree_nodes[1])
        current_parent[0].add_child(temp_tree_nodes[2])
        current_parent[0].taxon.label = id(current_parent[0])  # used to be .oid here. Removing oid causes a maximum recursion depth error
        # This part runs when you do the second step, ie, use the -m dmc flag

    out = tree.as_string('newick')
    out =  out.replace('[&U]','[50]') ## 50 unit long sequences
    open(outfile,'w').write(out)    # The first T produced by the history
    print(tree.as_ascii_plot())



#==============================================================================
# Once the we have a Original graph G and measure the topological
# distance from it. The topolofical distance is their likelihood
# Make Tree from It
#==============================================================================

def measure_D_net(G,qmod,qcon):
    D_net_dic = {}
    D_net_ret = {}
    D_net = []
    for u in G: D_net_dic[u] = {}

    for u in sorted(G):
        key1 = "Taxon" + str(u)
        tmp_row = []
        for v in sorted(G):
            key2 = "Taxon" + str(v)
            if u < v: continue
            D_net_dic[u][v] = 1.0 - G.dmc_likelihood(u,v,qmod,qcon)
            tmp_row.append(D_net_dic[u][v])

            print D_net_dic[u][v],
        D_net.append(tmp_row)
        print '\n'


    names = []
    for u in G: names.append('Taxon'+str(u))
    print names
    print D_net
    D_net_final = _DistanceMatrix(names,D_net)
    #print D_net_final.names

    constructor = DistanceTreeConstructor()
    tree_dmc = constructor.upgma(D_net_final)
    #print tree_dmc
    Phylo.write(tree_dmc,'ph_dmc.nre','newick')

    return D_net_final

#==============================================================================
#                           Make D_Seq
#    The fasta file is made from
#'./indel-seq-gen --matrix HKY --outfile ' + dna_out + ' < ' + outfile
#==============================================================================



def D_seq_matrix(fasta_file):
    aln = AlignIO.read(fasta_file, 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree_seq = constructor.upgma(dm)
    #print tree_dmc
    Phylo.write(tree_seq,'ph_seq.nre','newick')
    print dm.names
    return dm


#==============================================================================
#                Make D_F =  \alpha D_Seq + (1 - \alpha) D_Seq
#==============================================================================

# If we wanted to change "how much" of each distance matrix is used in the final combo tree, I think it would be here.
def D_F_matrix(D_Seq,D_net,final_tree, alpha):

    names_Seq = D_Seq.names
    names_Net = D_net.names
    D_F = []
    D_F_names = []

    for key1 in names_Net:
        i = names_Net.index(key1)
        #print key1
        temp_row = []
        for j in range(0,i+1):


            key2 = names_Net[j]
            #print key2,
            if key1 in names_Net and key2 in names_Seq:
                if not key1 in D_F_names:
                    D_F_names.append(key1)
                i1 = names_Net.index(key1)
                j2 = names_Net.index(key2)                              # should be 1-alpha * D_net and alpha * D_seq
                new_val = ((1-alpha) * D_net[key1,key2]) + (alpha * D_Seq[key1,key2])  # alpha is set to 0.5
                #print new_val,                                          # we can change alpha to choose how much of D_Seq and D_net we want to use
                temp_row.append(new_val)                                 # although shouldn't the formula in new_val be the same as what's on line 540?
        #print temp_row
        D_F.append(temp_row)

    print D_F

    D_F_final = _DistanceMatrix(D_F_names,D_F)

    constructor = DistanceTreeConstructor()
    tree_D_F = constructor.upgma(D_F_final)
    #print tree_dmc
    Phylo.write(tree_D_F,final_tree,'newick')
    return D_F_final



#==============================================================================
#                                   MAIN
#==============================================================================
def main():

    #================================================================
    #                Handle arguments and options
    #================================================================
    start = time.time()
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(levelname)s: %(asctime)s -- %(message)s'
    )

    usage="usage: %prog [options] <network file> <model parameters>"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--model", action="store", dest="model", type="string",default="dmc",help="model to reverse with: dmc, ff, or pa")
    #parser.add_option("-e", "--edges", action="store", dest="edges", type="int", default=10,help="Number of edges")

    (options, args) = parser.parse_args()
    model = options.model
    # ===============================================================


    # =================== Run NetArch algorithm =====================
    if len(args) == 0:
        logging.critical("No graph specified. Exiting...")
        sys.exit(1)

    if model != "gen":
        G = Graph.read_graph(args[0])
        G1 = Graph.read_graph(args[0])

    # DMC
    if model == "dmc" or model == "DMC":
        if len(args) != 5:
            logging.critical("Illegal number of parameters (%i, expected 3): <network> <qmod> <qcon>. Exiting..." %(len(args)))
            sys.exit(1)

        qmod = float(args[1])
        qcon = float(args[2])
        outfile = str(args[3])
        alpha = float(args[4])
        #his = str(args[3])

        if qmod < 0 or qmod > 1:
            logging.critical("Parameter 'qmod' must be between 0 and 1. Exiting...")
            sys.exit(1)
        elif qcon < 0 or qcon > 1:
            logging.critical("Parameter 'qcon' must be between 0 and 1. Exiting...")
            sys.exit(1)

        logging.info("Running ReverseDMC...")
        print "#ANCHOR\tNODE_REMOVED"
        dmc_delorean(G,qmod,qcon) # add the history file to store the tree structure

        # Create the tree from the graph
        generate_tree_dendro(H,'rev_dmc.nre')
        #generate_tree(H)
        # From the graph generate the distance matrix using DMC measure
        D_net = measure_D_net(G1,qmod,qcon)
        # indel sequence generator will take
        # over from here and it will create
        # a fasta file that contain the dequences
        # with Taxons made from tree in the outfile,
        dna_out = outfile.replace('.nre','') + '_DNA'  # so this step is renaming 'treefile'.nre to 'treefile'_DNA.fasta
        # rather, this step is changing the object which just stores the file name.

        ind_seq = './indel-seq-gen --matrix HKY --outfile ' + dna_out + ' < ' + outfile  # This isn't a bash-style write arrow--the help page for indel-seq-gen says to use this when referring to the tree file.
        os.system(ind_seq)  # indel-seq-gen generates a collection of different files, a .seq file being one of them.


        new_seq = dna_out + ".seq"  # the dna_out object is just "ph1_DNA". there is no actual data in it. But rather there is a ph1_DNA file that is also generated

        #fasta_file = new_seq
        os.rename(new_seq,new_seq.replace('seq','fasta'))  # this line changes the .seq extension to .fasta
        fasta_file = new_seq.replace('seq','fasta')        # os.rename should be making ph1_DNA.fasta


        # Pass fasta file here to generate the
        # sequence distance
        D_seq = D_seq_matrix(fasta_file)

        # We have two distances we can make
        # convex combination of those and then
        # use them to generate the new distance matrix
        # Create a tree from it UPGMA model would work
        # here ..
        final_tree = new_seq.replace('.seq','') + '_F' + '.nre'
        D_F = D_F_matrix(D_seq,D_net,final_tree, alpha)
        os.system('rm *.anc_tree *.ma *.trace *.root *.scale_tree')

    elif model == "gen":
        nodes = int(args[0])
        edges = int(args[1])
        outfile = str(args[2])
        qmod = float(args[3])
        qcon = float(args[4])
        treefile = str(args[5])
        #G = Graph.rand_graph(nodes,edges,outfile)
        G = forward_dmc(nodes,edges,outfile,qmod,qcon,treefile)


    else:
        logging.critical("Invalid model: %s. Exiting..." %(model))
        sys.exit(1)


    # ========================= Finish ===========================
    logging.info("Time to run: %.2f (mins)" %((time.time()-start) / 60))


if __name__ == "__main__":
    main()

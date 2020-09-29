import networkx as nx
import sys
#import igraph as gr

inputfile=sys.argv[1]

#read the input .edg file and save it as Graph
G=nx.Graph()
fin=open(inputfile,'r')

count=1
nodehash={}
for line in fin:
    splitted=line.split()
    if splitted[0] not in nodehash:
         nodehash[splitted[0]]=count  
         count=count+1
    if splitted[1] not in nodehash:
         nodehash[splitted[1]]=count
         count=count+1
fin.close()


nodecount=len(nodehash.keys())
fin=open(inputfile,'r')
for line in fin:
    splitted=line.split()
    if (len(splitted)!=0):
       node1=nodehash[splitted[0]]
       node2=nodehash[splitted[1]]
       G.add_node(node1)
       G.add_node(node1)
       G.add_edge(node1,node2)
fin.close()  


nx.write_edgelist(G,"modified_out.txt")



#G2=gr.Read(klass, f, format=None, *args, **kwds)
#vertex_cluster=G2.community_fastgreedy(self, weights=None)

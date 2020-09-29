import networkx as nx
import sys
#import igraph as gr

inputfile="modularity.mod"
outfile="metric.mod"

#read the input .edg file and save it as Graph
G=nx.Graph()
fin=open(inputfile,'r')
fout=open(outfile,'w')
for line in fin:
    if line.startswith("solve;"):
        #Extra knowledge
    fout.write(line+"\n") 
  
fin.close()
fout.close()

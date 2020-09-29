import networkx as nx
import sys

inputfile=sys.argv[1]
outputfile=sys.argv[2]

#read the input .edg file and save it as Graph
G=nx.Graph()
fin=open(inputfile,'r')

readmode=0
cost={}
distance={}
labels=[]
for line in fin:
    if ( line.startswith("Labels") ):
        splitted=line.split()
        labels=splitted[1:]
        continue
    if ( line.startswith("Cost") ):
        readmode=1
        continue
    if ( line.startswith("Distance") ):
        readmode=2
        continue
    
    if(readmode==0):   
        splitted=line.split()
        if (len(splitted)!=0):
           G.add_node(splitted[0])
           G.add_node(splitted[1])
           G.add_edge(splitted[0],splitted[1])
    
     elif(readmode==1):
         splitted=line.split()
         cost[splitted[0]]=splitted[1:]
         
     elif(readmode==2):    
         splitted=line.split()
         distance[splitted[0]]=splitted[1:]
          
fin.close()

nodes=len(G)
edges=G.size()

if (len(cost.keys())!=nodes ):
   print "Number of nodes is not equal to size of cost matrix"
   sys.exit()

if (len(distance.keys())!=len(labels) ):
   print "Number of labels is not equal to size of distance matrix"
   sys.exit()


#edgestr will hold all edges in the graph as an adjacency matrix for data file
edgestr="";
for node in G.nodes():
    edgestr=edgestr+" "+str(node)
edgestr=edgestr+" := \n"

total=len(G.nodes())
count=0
for node1 in G.nodes():
    edgestr=edgestr+str(node1)+" "
    for node2 in G.nodes():
        if G.has_edge(node1,node2)==True:
           edgestr=edgestr+str(1)+" "
        else:
           edgestr=edgestr+str(0)+" "
    count=count+1
    if count==total:
       edgestr=edgestr+"; \n"
    else:
       edgestr=edgestr+" \n"

#coststr will hold all the costs associated with assigning each vertex a label. This matrix will be size of nodesXlabels.
coststr="";
for label in range(1,len(labels)+1): 
    coststr=coststr+" "+str(label)
coststr=coststr+" := \n"

total=len(G.nodes())
count=0
for node1 in G.nodes():
    coststr=coststr+str(node1)+" "
    for label in range(1,len(labels)+1):
           coststr=coststr+str(cost[node1][label-1])+" "        
    count=count+1
    if count==total:
       coststr=coststr+"; \n"
    else:
       coststr=coststr+" \n"

#distancestr will hold all the distances associated with labels to neighbour vertices. This matrix will be size of labelsXlabels.
distancestr="";
for label in range(1,len(labels)+1): 
    distancestr=distancestr+" "+str(label)
distancestr=distancestr+" := \n"

total=len(labels)
count=0
for label in range(1,len(labels)+1):
    distancestr=distancestr+str(label)+" "
    for label2 in range(1,len(labels)+1):
           distancestr=distancestr+str(cost[label][label2-1])+" "        
    count=count+1
    if count==total:
       distancestr=distancestr+"; \n"
    else:
       distancestr=distancestr+" \n"


#Create and write to the output data file
fout = open(outputfile, "w")

str1="param N := "+str(nodes)+" ;\n"
fout.write(str1)
str1="param M := "+str(edges)+" ;\n"
fout.write(str1)
str1="param L := "+str(len(labels))+" ;\n"
fout.write(str1)
str1="param EDGES : "+edgestr
fout.write(str1)
str1="param COST : "+coststr
fout.write(str1)
str1="param DISTANCE : "+distancestr
fout.write(str1)
str1="end; \n"
fout.write(str1)
fout.close()



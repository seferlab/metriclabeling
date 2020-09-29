# EMRE SEFER 06.04.2010
# Run as  "python ./associationparser.py gene_association.fb"
# To write output into a file use -> out.txt
import networkx as nx
import sys
import string
import random

#if ( len(sys.argv) != 2 ):
#     print "Illegal number of parameters: "
#     sys.exit(1)


#read the input ppi network file and save it as Graph

Ontology=nx.DiGraph()
Bio_Ontology=nx.DiGraph()
Molecular_Ontology=nx.DiGraph()
Cellular_Ontology=nx.DiGraph()

filename="gene_ontology_ext.obo.obo"
file=open(filename,'r')
flag=False
node=-10
part_count=0
edge=0
numnode=0
namespace=-1
name=""
Bio_name={}
Cellular_name={}
Molecular_name={}
for line in file:
    if line.strip()=="[Term]":
       flag=True
    elif (line.startswith("id:") and flag):
        splitted=line.split()
        node=int(splitted[1][3:])
        Ontology.add_node(node)
        numnode=numnode+1
    elif (line.startswith("name:") and flag):
        splitted=line.split()
        name=""
        for elem in splitted[1:]:
          name=name+elem+" "
    elif (line.startswith("namespace:") and flag):
        splitted=line.split()
        if splitted[1]=="biological_process":
            namespace=1
            Bio_Ontology.add_node(node)
            Bio_name[node]=name
        elif splitted[1]=="cellular_component":
            namespace=2
            Cellular_Ontology.add_node(node)
            Cellular_name[node]=name
        elif splitted[1]=="molecular_function":
            namespace=3
            Molecular_Ontology.add_node(node)
            Molecular_name[node]=name 
        else:
            print "There is an error here!!"        
    elif (line.startswith("is_a:") and flag):
        splitted=line.split()
        newnode=int(splitted[1][3:])
        Ontology.add_weighted_edges_from([(newnode,node,1.0)])
        if namespace==1:
           Bio_Ontology.add_weighted_edges_from([(newnode,node,1.0)])
        elif namespace==2:
           Cellular_Ontology.add_weighted_edges_from([(newnode,node,1.0)])
        elif namespace==3:
           Molecular_Ontology.add_weighted_edges_from([(newnode,node,1.0)])
        edge=edge+1
    elif (line.startswith("relationship:") and flag):
        splitted=line.split()
        if (splitted[1]=="part_of"):
           part_count=part_count+1
    elif not line.strip():
       flag=False 
       namespace=-1
    elif ( line.startswith("[Typedef]") and not flag): # For handling end of file
        break
file.close()

print "GO INFORMATION"
print "Ontology size is: "+str(len(Ontology))
#print "Ontology nodes are: "+str(len((Ontology.nodes(data=False))))
print "Part count is: "+str(part_count)
print "Edges are: "+str(edge)
#print "Number of nodes are: "+str(numnode)

print "Biosize is "+str(len(Bio_Ontology))
print "Cellularsize is "+str(len(Cellular_Ontology))
print "Molecularsize is "+str(len(Molecular_Ontology))

bioroot=-1
molecularroot=-1
cellularroot=-1

#Find roots for all types
for node in (Bio_Ontology.nodes(data=False)):
   if ( Bio_Ontology.out_degree(node,weighted=True)!=0 and Bio_Ontology.in_degree(node,weighted=True)==0 ):
      print "Root node for Biological sub"
      print node
      bioroot=node
      #print Bio_Ontology.out_degree(node,weighted=True)
      #print Bio_Ontology.in_degree(node,weighted=True) 

for node in (Molecular_Ontology.nodes(data=False)):
   if ( Molecular_Ontology.out_degree(node,weighted=True)!=0 and Molecular_Ontology.in_degree(node,weighted=True)==0 ):
      print "Root node for Molecular sub"
      print node
      molecularroot=node
      #print Molecular_Ontology.out_degree(node,weighted=True)
      #print Molecular_Ontology.in_degree(node,weighted=True) 

for node in (Cellular_Ontology.nodes(data=False)):
   if ( Cellular_Ontology.out_degree(node,weighted=True)!=0 and Cellular_Ontology.in_degree(node,weighted=True)==0 ):
      print "Root node for Cellular sub"
      print node
      cellularroot=node
      #print Cellular_Ontology.out_degree(node,weighted=True)
      #print Cellular_Ontology.in_degree(node,weighted=True) 

#for key in Bio_name.keys():
#    print str(key)+" --> "+Bio_name[key]
 
#for key in Cellular_name.keys():
#    print str(key)+" --> "+Cellular_name[key]

#for key in Molecular_name.keys():
#    print str(key)+" --> "+Molecular_name[key]

#Undirected graph version
Ontology2=Ontology.to_undirected()
Bio_Ontology2=Bio_Ontology.to_undirected()
Molecular_Ontology2=Molecular_Ontology.to_undirected()
Cellular_Ontology2=Cellular_Ontology.to_undirected()

#print Bio_Ontology.edges()
#print Bio_name.values()

max=0;
for elem in nx.connected_components(Bio_Ontology2):
    if len(elem)==19022:
        SubBio=nx.subgraph(Bio_Ontology2,elem)
print max
print len(SubBio)
        

#tot=0
#for key in SubBio.nodes():
#    if nx.shortest_path_length(SubBio,key,bioroot) >= 8:
#        tot=tot+1
#print tot

length=300
randindex=[]
FILE = open("out"+str(length)+".txt","w")
for i in range(0,length):
    val=random.randrange(0,len(Bio_Ontology2))
    while (val not in SubBio.nodes()) or val in randindex:
        val=random.randrange(0,len(Bio_Ontology2))
    randindex.append(val)
print "done"

tr=""
for i in range(0,len(randindex)):
    tr=tr+str(randindex[i])+" "
FILE.write(tr+"\n")    


distance = [float] * len(randindex)
for i in range(len(randindex)):
    distance[i] = [float] * (len(randindex))

for i in range(0,len(randindex)):
    distance[i][i]=0.0

for i in range(0,len(randindex)):
    for j in range(i+1,len(randindex)):
        path1=nx.shortest_path(Bio_Ontology2,randindex[i],bioroot)
        path2=nx.shortest_path(Bio_Ontology2,randindex[j],bioroot)
        inter = [val for val in path1 if val in path2]
        maxsize=-1
        maxelem=-1
        for elem in inter:
            if nx.shortest_path_length(Bio_Ontology2,elem,bioroot) > maxsize :
                maxsize=nx.shortest_path_length(Bio_Ontology2,elem,bioroot)
                maxelem=elem
        b=nx.shortest_path_length(Bio_Ontology2,maxelem,randindex[i])
        c=nx.shortest_path_length(Bio_Ontology2,maxelem,randindex[j])
        if maxsize==0:
            val=((b+c)/2)
        else:
            val=(b+c)/(2*maxsize+2)
        distance[i][j]=val
        distance[j][i]=val


for i in range(0,len(randindex)):
    tempstr=""
    for j in range(0,len(randindex)):
       tempstr=tempstr+str(distance[i][j])+" "
    FILE.write(tempstr+"\n")   
    
FILE.close()


exit(1)

print "\n"
print "RESULTS FOR A THALIANA FOR ATPID DATABASE"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.tair"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="march_2010_download.txt"
file=open(filename)
count=1
flag=0
#aliases={}
#experiments={}
for line in file:
    splitted=line.split("\t")
    first=splitted[0].lower()
    second=splitted[1].lower()    
    #splitted[4]="N/A"
    #splitted[5]="N/A"

    #if line.startswith("INTERACTOR_A") and flag==0 :
    #    flag=1
    #    continue
    #if flag==1: 
    #   splitted=line.split("\t")
    #   first=splitted[2].lower()
    #   second=splitted[3].lower()
    #   aliaslist1=splitted[4].split("|")
    #   aliaslist2=splitted[5].split("|")
    #   if splitted[6] in experiments:
    #      experiments[splitted[6]]=experiments[splitted[6]]+1
    #   else:
    #      experiments[splitted[6]]=1

    if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          #if splitted[4]!="N/A": 
          #    aliases[first]=map(string.lower, aliaslist1)
    if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          #if splitted[5]!="N/A": 
          #    aliases[second]=map(string.lower, aliaslist2)

file.close()

#print "Experiment Types are:"
#for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
#    print key+" "+str(val)



found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    splitted=line.split("\t")
    first=splitted[0].lower()
    second=splitted[1].lower()
    node1=proteins_to_count[first]
    node2=proteins_to_count[second]
    G.add_node(node1)
    G.add_node(node1)
    G.add_edge(node1,node2)
    #if line.startswith("INTERACTOR_A") and flag==0 :
    #    flag=1
    #    continue
    #if flag==1: 
    #   splitted=line.split("\t")
    #   if (len(splitted)!=0):
    #     first=splitted[2].lower()
    #     second=splitted[3].lower()
    #     node1=proteins_to_count[first]
    #     node2=proteins_to_count[second]
    #     G.add_node(node1)
    #     G.add_node(node1)
    #     G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    #else:
    #  flag=0
    #  if key in aliases:
    #    aliaslist=aliases[key]
    #    foundalias=""
    #    for alias in aliaslist:
    #      if alias in proteins_to_function.keys():
    #         foundalias=alias
    #         flag=1
    #         break
    #    if flag==1:   
    #      val=proteins_to_function[foundalias]
    #      flag1=0
    #      flag2=0
    #      flag3=0
    #      for elem in val:
    #        if elem in Bio_name.keys() and flag1==0:
    #          realbioannotated=realbioannotated+1
    #          flag1=1
    #        elif elem in Cellular_name.keys() and flag2==0:
    #          realcellannotated=realcellannotated+1
    #          flag2=1
    #       elif elem in Molecular_name.keys() and flag3==0:
    #          realmolannotated=realmolannotated+1
    #          flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    #else:
    #  flag=0
    #  if key in aliases:
    #    aliaslist=aliases[key]
    #    foundalias=""
    #    for alias in aliaslist:
    #      if alias in proteins_to_function.keys():
    #         foundalias=alias
    #         flag=1
    #         break
    #    if flag==1:   
    #      val=proteins_to_function[foundalias]
    #      flag1=0
    #      flag2=0
    #      flag3=0
    #      for elem in val:
    #        if elem in Bio_name.keys() and flag1==0:
    #          connbioannotated=connbioannotated+1
    #          flag1=1
    #        elif elem in Cellular_name.keys() and flag2==0:
    #          conncellannotated=conncellannotated+1
    #          flag2=1
    #        elif elem in Molecular_name.keys() and flag3==0:
    #          connmolannotated=connmolannotated+1
    #          flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR D MELAGONASTER"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.fb"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Drosophila_melanogaster-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)

found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))



biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)






print "\n"
print "RESULTS FOR S CEREVISIAE"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.sgd"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)

found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)




print "\n"
print "RESULTS FOR C ELEGANS"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.wb"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Caenorhabditis_elegans-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR A THALIANA"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.tair"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Arabidopsis_thaliana-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    # I commented these parts
    #splitted=line.split("\t")
    #first=splitted[0]
    #second=splitted[1]
    #splitted[4]="N/A"
    #splitted[5]="N/A"

    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR S POMBE"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.GeneDB_Spombe"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Schizosaccharomyces_pombe-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR MUS MUSCULUS"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.mgi"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Mus_musculus-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR HUMAN"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.goa_human"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Homo_sapiens-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)



print "\n"
print "RESULTS FOR RATTUS NORVEGICUS"
#S cerevisiae
G=nx.Graph()
proteins_to_function={}
#filename=sys.argv[1]
filename="gene_association.rgd"
file=open(filename)
for line in file:
   if line.startswith("!"): 
       continue
   parts=line.split("\t")
   #gene=parts[1][4:]
   gene=parts[2].lower()
   onto=""
   if ( parts[3][0:2]=="GO"): 
      onto=int(parts[3][3:])
   else:
      onto=int(parts[4][3:])  
   #Since each protein can have more than one function, for each gene we keep a list of functions it has
   # We must make sure that this function has not been added to list before.
   if ( proteins_to_function.has_key(gene) ):
       if ( not (onto in proteins_to_function[gene] ) ):
           proteins_to_function[gene].add(onto)
   else: 
       proteins_to_function[gene]=set()
       proteins_to_function[gene].add(onto)
file.close()

proteins_to_count={}
count_to_proteins={}
#filename=sys.argv[2]
filename="BIOGRID-ORGANISM-Rattus_norvegicus-3.0.65.tab.txt"
file=open(filename)
count=1
flag=0
aliases={}
experiments={}
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       first=splitted[2].lower()
       second=splitted[3].lower()
       aliaslist1=splitted[4].split("|")
       aliaslist2=splitted[5].split("|")
       if splitted[6] in experiments:
          experiments[splitted[6]]=experiments[splitted[6]]+1
       else:
          experiments[splitted[6]]=1

       if first not in proteins_to_count:
          proteins_to_count[first]=count  
          count_to_proteins[count]=first
          count=count+1
          if splitted[4]!="N/A": 
              aliases[first]=map(string.lower, aliaslist1)
       if second not in proteins_to_count:
          proteins_to_count[second]=count
          count_to_proteins[count]=second
          count=count+1
          if splitted[5]!="N/A": 
              aliases[second]=map(string.lower, aliaslist2)

file.close()

print "Experiment Types are:"
for (val,key) in sorted([(value,key) for (key,value) in experiments.items()], reverse=True):
    print key+" "+str(val)


found=0
notfound=0
for elem in proteins_to_count.keys():
    if elem in proteins_to_function.keys():     
       found=found+1
    else:
       flag=0
       if elem in aliases:
          aliaslist=aliases[elem]
          for alias in aliaslist:
            if alias in proteins_to_function.keys():
              flag=1
          if flag==1:   
             found=found+1
          else:
             notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)


found=0
notfound=0
for elem in proteins_to_function.keys():
    if elem in proteins_to_count.keys():     
       found=found+1
    else:  
       notfound=notfound+1

print "Found number is" + str(found)
print "Unfound number is" +  str(notfound)



#nodecount=len(nodehash.keys())
file=open(filename,'r')
flag=0
for line in file:
    if line.startswith("INTERACTOR_A") and flag==0 :
        flag=1
        continue
    if flag==1: 
       splitted=line.split("\t")
       if (len(splitted)!=0):
         first=splitted[2].lower()
         second=splitted[3].lower()
         node1=proteins_to_count[first]
         node2=proteins_to_count[second]
         G.add_node(node1)
         G.add_node(node1)
         G.add_edge(node1,node2)
file.close() 

list=nx.connected_components(G)
max=0
maxlist=[]
for elem in list:
    if len(elem) > max:
       max=len(elem)
       maxlist=elem

print "Number of proteins on largest connected component is: "+str(max)
print "Number of proteins in all network is: "+str(len(proteins_to_count.keys()))
print "Number of interactions on ppi: "+str(G.number_of_edges())
Sub=nx.subgraph(G, maxlist )
print "Number of interactions on connected ppi: "+str(Sub.number_of_edges())

print "Clustering coefficient is: "+str(nx.average_clustering(Sub))
print "Transitivity is: "+str(nx.transitivity(Sub))
print "Graph clique number: "+str(nx.graph_clique_number(Sub))
print "Diameter is: "+str(nx.diameter(Sub))
print "Eccentrcity is: "+str(len(nx.center(Sub)))

print "Clustering coefficient is: "+str(nx.average_clustering(G))
print "Transitivity is: "+str(nx.transitivity(G))
print "Graph clique number: "+str(nx.graph_clique_number(G))


biofunc=0
cellfunc=0
molfunc=0

allfunc=set()
#ATcount=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    allfunc.update(val)
    #if key.startswith("AT"):
    #   ATcount=ATcount+1

print "All function count is:"+str(len(allfunc))
#print "AT count is: "+str(ATcount)

for elem in allfunc:
       if elem in Bio_name.keys():
          biofunc=biofunc+1
       elif elem in Cellular_name.keys():
          cellfunc=cellfunc+1
       elif elem in Molecular_name.keys():
          molfunc=molfunc+1

print "Bio num: "+str(biofunc)
print "Cell num "+str(cellfunc)
print "Molecular num "+str(molfunc)

bioannotated=0
cellannotated=0
molannotated=0
for key in proteins_to_function.keys():
    val=proteins_to_function[key]
    flag1=0
    flag2=0
    flag3=0
    for elem in val:
        if elem in Bio_name.keys() and flag1==0:
            bioannotated=bioannotated+1
            flag1=1
        elif elem in Cellular_name.keys() and flag2==0:
            cellannotated=cellannotated+1
            flag2=1
        elif elem in Molecular_name.keys() and flag3==0:
            molannotated=molannotated+1
            flag3=1

realbioannotated=0
realcellannotated=0
realmolannotated=0

for key in proteins_to_count.keys():
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            realbioannotated=realbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            realcellannotated=realcellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            realmolannotated=realmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              realbioannotated=realbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              realcellannotated=realcellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              realmolannotated=realmolannotated+1
              flag3=1

connbioannotated=0
conncellannotated=0
connmolannotated=0

proteinlist=[]
for elem in maxlist:
    proteinlist.append(count_to_proteins[elem])

for key in proteinlist:
    if key in proteins_to_function.keys():
        val=proteins_to_function[key]
        flag1=0
        flag2=0
        flag3=0
        for elem in val:
          if elem in Bio_name.keys() and flag1==0:
            connbioannotated=connbioannotated+1
            flag1=1
          elif elem in Cellular_name.keys() and flag2==0:
            conncellannotated=conncellannotated+1
            flag2=1
          elif elem in Molecular_name.keys() and flag3==0:
            connmolannotated=connmolannotated+1
            flag3=1
    else:
      flag=0
      if key in aliases:
        aliaslist=aliases[key]
        foundalias=""
        for alias in aliaslist:
          if alias in proteins_to_function.keys():
             foundalias=alias
             flag=1
             break
        if flag==1:   
          val=proteins_to_function[foundalias]
          flag1=0
          flag2=0
          flag3=0
          for elem in val:
            if elem in Bio_name.keys() and flag1==0:
              connbioannotated=connbioannotated+1
              flag1=1
            elif elem in Cellular_name.keys() and flag2==0:
              conncellannotated=conncellannotated+1
              flag2=1
            elif elem in Molecular_name.keys() and flag3==0:
              connmolannotated=connmolannotated+1
              flag3=1       

print "Number of annotated proteins is: "+str(len(proteins_to_function.keys()))
print "Number all proteins seen in the ppi network is: "+str(len(proteins_to_count.keys()))
print "IN ONTOLOGY; Number of proteins that are Biological function annotated: "+str(bioannotated)
print "IN ONTOLOGY; Number of proteins that are Cell annotated:  "+str(cellannotated)
print "IN ONTOLOGY; Number of proteins that are Molecular annotated: "+str(molannotated)
print "IN PPI; Number of proteins that are Biological function annotated: "+str(realbioannotated)
print "IN PPI; Number of proteins that are Cell annotated:  "+str(realcellannotated)
print "IN PPI; Number of proteins that are Molecular annotated: "+str(realmolannotated)
print "IN Largest Connected PPI; Number of proteins that are Biological function annotated: "+str(connbioannotated)
print "IN Largest Connected PPI; Number of proteins that are Cell annotated:  "+str(conncellannotated)
print "IN Largest Connected PPI; Number of proteins that are Molecular annotated: "+str(connmolannotated)






exit(1)

#TEST PHASE
print len(proteins_to_function.keys())
print len(proteins_to_count.keys()) 

count=0
for key in proteins_to_functions.keys():
    if key not in proteins_to_count.keys():
        count=count+1 

print "Number of nonmatch is: "+str(count)




for key in  proteins_to_function.keys():
   myset= proteins_to_function[key]
   sub=""
   for elem in myset:
      sub=sub+" "+elem
   print key+" -> "+sub


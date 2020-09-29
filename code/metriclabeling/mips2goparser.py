# EMRE SEFER 06.04.2010
# Run as  "python ./associationparser.py gene_association.fb"
# To write output into a file use -> out.txt
import networkx as nx
import sys
import string

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




filename="mips2go.txt"
file=open(filename,'r')
mipsfonk={}
count=0
linecount=0
unknowncount=0
cellcount=0
biocount=0
molecularcount=0
othercount=0
for line in file:
    if ( not line.startswith("!") ):
          linecount=linecount+1
          splitted=line.split()
          mipscode=splitted[0].split(":")
          flag=0
          for part in splitted:
                if (part.startswith("GO:0")):
                    mipsfonk[mipscode[1]]=part
                    if (int(part[3:]) in Bio_name.keys() ):
                       biocount=biocount+1
                    elif (int(part[3:]) in Cellular_name.keys() ):
                       cellcount=cellcount+1
                    elif (int(part[3:]) in Molecular_name.keys() ): 
                       molecularcount=molecularcount+1
                    else:
                       othercount=othercount+1
                    flag=1
          if flag==0:
             for part in splitted:
                if (part.startswith("GO:.")):
                    unknowncount=unknowncount+1
             count=count+1
file.close()

print "Count number budur: "+str(count)
print "Linecount number budur: "+str(linecount)
print "Unknowncount number budur: "+str(unknowncount)
print biocount
print molecularcount
print cellcount 
print othercount
#print "Hasteki element sayisi budur: "+str(len(mipsfonk.keys()))"

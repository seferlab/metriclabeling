
import networkx as nx
import sys
import string
import random
import numpy

nodenum=20
edgenum=80
labelnum=90
edges = numpy.zeros((nodenum,nodenum))
cost = numpy.zeros((labelnum,nodenum))
distance = numpy.zeros((labelnum,labelnum))

for i in range(0,nodenum): 
    for j in range(0,i+1):
       edges[i,j]=random.randrange(0,2)       
       edges[j,i]=edges[i,j]

for i in range(0,nodenum):
       edges[i,i]=1

 
for i in range(0,labelnum): 
    for j in range(0,nodenum):
       cost[i,j]=2+random.random() 

for i in range(0,labelnum): 
    for j in range(0,labelnum):
       distance[i,j]=2+3*random.random() 


filename="metric2.data"
file=open(filename,"w")
file.write("param N := "+str(nodenum)+";\n")
file.write("param M := "+str(edgenum)+";\n")
file.write("param LabelNum := "+str(labelnum)+";\n")
     
file.write("param EDGES: ")
for i in range(0,nodenum):
     file.write(str(i)+" ")
file.write(" := \n")

for i in range(0,nodenum): 
    file.write(str(i)+" ")
    for j in range(0,nodenum):
        file.write(str(edges[i][j])+" ")  
    if i!=(nodenum-1):
        file.write(" \n")
    else:
        file.write(" ; \n")  

file.write("param COST: ")
for i in range(0,nodenum):
         file.write(str(i)+" ")
file.write(" := \n")

for i in range(0, labelnum): 
        file.write(str(i)+" ")
        for j in range(0,nodenum):
             file.write(str(cost[i][j])+" ")  
        if i!=(labelnum-1):
            file.write(" \n")
        else:
            file.write(" ; \n") 

file.write("param DISTANCE: ")
for i in range(0,labelnum):
         file.write(str(i)+" ")
file.write(" := \n")

for i in range(0, labelnum): 
        file.write(str(i)+" ")
        for j in range(0,labelnum):
             file.write(str(distance[i][j])+" ")  
        if i!=(labelnum-1):
           file.write(" \n")
        else:
            file.write(" ; \n") 



file.close()
  




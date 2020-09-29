import sys


inputfile=sys.argv[1]
outputfile="out.txt"
count=0
depohash={}
#read the input .edg file and save it as Graph
fin=open(inputfile,'r')
fout=open(outputfile,'w')
for line in fin:
    if line.startswith("!"):
       continue
    splitted=line.split()
    count=count+1
    depohash[splitted[2]]=splitted[3]
    #fout.write(splitted[2]+" --> "+splitted[3]+"\n") 
  
print count

list=sorted(depohash.keys())
for elem in list:
   fout.write(elem+" --> "+depohash[elem]+"\n") 

fin.close()
fout.close()


import sys

length=300
distance = [float] * length
for i in range(length):
    distance[i] = [float] * (length)

filename="out300.txt"
f = open(filename, "r")
x=0
for line in f:
   if x==0:
       randindex=line.rstrip().split(" ")    
   else:
       temp=line.rstrip().split(" ")
       for i in range(0,len(temp)):
           distance[x-1][i]=float(temp[i])
   x=x+1
   
count=0
newdistance = [float] * len(randindex)
for i in range(len(randindex)):
    newdistance[i] = [float] * (len(randindex))

for i in range(0,len(randindex)):
    newdistance[i][i]=0.0
    
for i in range(0,len(randindex)):
    for j in range(i+1,len(randindex)):
        for k in range(j+1,len(randindex)):
             num1=i
             num2=j
             num3=k
             temp=0
             if (distance[num1][num2]+distance[num2][num3]<distance[num1][num3]):
                 temp=1;
             elif (distance[num1][num3]+distance[num2][num3]<distance[num1][num2]):
                 temp=1;
             elif (distance[num1][num2]+distance[num1][num3]<distance[num2][num3]):
                 temp=1;    
             count=count+temp   

print "Unsatisfying triangles are: "+str(count)

for i in range(0,len(randindex)):
    for j in range(i+1,len(randindex)):
        minval=distance[i][j]
        for k in range(0,len(randindex)):
            if k!=i and k!=j:
             num1=i
             num2=j
             num3=k
             if (distance[num1][num3]+distance[num2][num3]<minval):
                 minval=distance[num1][num3]+distance[num2][num3]         
        newdistance[num1][num2]=minval
        newdistance[num2][num1]=minval

errortotal=0.0
for i in range(0,len(randindex)):
    for j in range(0,len(randindex)):
        num1=i
        num2=j
        print repr(newdistance[num1][num2])+" " 
        errortotal=errortotal+(newdistance[num1][num2]-distance[num1][num2])*(newdistance[num1][num2]-distance[num1][num2])
    print "\n"
    
print errortotal/2

tr=""
FILE = open("newdistanceout300.out","w")
for i in range(0,len(randindex)):
    for j in range(0,len(randindex)):
       tr=tr+str(newdistance[i][j])+" "
    FILE.write(tr+"\n")  
FILE.close()

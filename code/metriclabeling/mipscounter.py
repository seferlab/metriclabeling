import sys

filename="MappingFunCat2.0_2.1.txt"
file=open(filename)
count=0
level1=0
level2=0
level3=0
hash1=[]
hash2=[]
hash3=[]
namehash1={}
namehash2={}
namehash3={}
for line in file:
   if count<3 : 
       count=count+1
   else:
      parts=line.split()
      name=parts[0]
      remain=""
      for elem in parts[1:]:
         remain=remain+" "+elem    
      num=name.count(".")

      if num==0:
         if remain.count("->")==0 :
             hash1.append(name)
             namehash1[name]=remain
             level1=level1+1
      if num==1:
         if remain.count("->")==0 :
             hash2.append(name)
             namehash2[name]=remain
             level2=level2+1
      if num==2:
         if remain.count("->")==0 :
             hash3.append(name)
             namehash3[name]=remain
             level3=level3+1
file.close()

print level1      
print level2
print level3

#print len(namehash1.keys())
#print len(namehash2.keys())
#print len(namehash3.keys())

#LEVEL 3 NEEDED ADDITIONS
extra=[]
extra.append("02.04")
extra.append("02.05")
extra.append("02.08")
extra.append("02.09")
extra.append("02.10")
extra.append("02.17")
extra.append("02.19")
extra.append("02.25")
extra.append("04.01")
extra.append("04.02")
extra.append("12.07")
extra.append("12.10")
extra.append("12.11")
extra.append("14.01")
extra.append("14.04")
extra.append("14.10")
extra.append("16.02")
extra.append("16.05")
extra.append("16.06")
extra.append("16.07")
extra.append("16.09")
extra.append("16.10")
extra.append("16.11")
extra.append("16.12")
extra.append("16.14")
extra.append("16.15")
extra.append("16.16")
extra.append("16.18")
extra.append("16.22")
extra.append("16.23")
extra.append("16.25")
extra.append("30.07")
extra.append("38.01")
extra.append("38.02")
extra.append("38.03")
extra.append("38.05")
extra.append("38.06")
extra.append("38.07")
extra.append("40.05")
extra.append("40.20")
extra.append("42.01")
extra.append("42.03")
extra.append("42.05")
extra.append("42.07")
extra.append("42.08")
extra.append("42.09")
extra.append("42.19")
extra.append("42.22")
extra.append("42.25")
extra.append("42.28")
extra.append("42.29")
extra.append("42.30")
extra.append("42.32")
extra.append("42.33")
extra.append("42.36")
extra.append("42.37")
extra.append("70.01")
extra.append("70.03")
extra.append("70.05")
extra.append("70.07")
extra.append("70.08")
extra.append("70.09")
extra.append("70.19")
extra.append("70.22")
extra.append("70.25")
extra.append("70.28")
extra.append("70.29")
extra.append("70.30")
extra.append("70.32")
extra.append("70.33")
extra.append("70.37")

extrahash={}
for elem in extra:
   file=open(filename)
   count=0
   for line in file:
      if count<3 : 
        count=count+1
      else:
        parts=line.split()
        name=parts[0]
        if name==elem: 
           remain=""
           for elem in parts[1:]:
             remain=remain+" "+elem  
             extrahash[name]=remain 
           break 
   file.close()


for key in sorted(extrahash.keys()):
    namehash3[key]=extrahash[key]

#print len(namehash3)

outfile = open("mipslevel1.txt","w")
for key in sorted(namehash1.keys()):
   outfile.write(key+" "+namehash1[key]+"\n")
outfile.close()

outfile = open("mipslevel2.txt","w")
for key in sorted(namehash2.keys()):
   outfile.write(key+" "+namehash2[key]+"\n")
outfile.close()

outfile = open("mipslevel3.txt","w")
for key in sorted(namehash3.keys()):
   outfile.write(key+" "+namehash3[key]+"\n")
outfile.close()

# LEVEL3 ekstradan 71 tane daha ekleme yukarda yapildi

#count=1
#for elem in hash2:
#   print str(count)+" "+hash2[count-1]
#   count=count+1

#print hash1
#print hash2
#print hash3

filename="Saccharomyces_cerevisiae.funcat.funcat"
file=open(filename)
annotatenum=0
level2={}
level3={}
for line in file:
   if (line.startswith(">")):
      annotatenum=annotatenum+1
      parts=line.split()
      mips=parts[1]
      func=mips.split(";")
      for elem in func:
         num=elem.count(".")
         if num>=2:
            levelparts=elem.split(".")
            level2[levelparts[0]+"."+levelparts[1]]=1
            level3[levelparts[0]+"."+levelparts[1]+"."+levelparts[2]]=1
         elif num==1:
             level2[elem]=1;
             #level3
         #elif num==0:
              
file.close()

print "Mips Cerevisiae annotation number:"+str(annotatenum)
print "Cerevisiae Level2 annotation number:"+str(len(level2.keys()))
print "Cerevisiae Level3 annotation number:"+str(len(level3.keys()))  


filename="athaliana_funcat2008."
file=open(filename)
level2={}
level3={}
flag=0
namehash={}
extracount=0
num2=0
num1=0
num0=0
for line in file:
   if (line.startswith("A") and flag==0):
      flag=1
   if flag==1:
      parts=line.split("|")
      namehash[parts[0]]=1
      elem=parts[1]
      num=elem.count(".")
      if num>=2:
            num2=num2+1
            levelparts=elem.split(".")
            level2[levelparts[0]+"."+levelparts[1]]=1
            level3[levelparts[0]+"."+levelparts[1]+"."+levelparts[2]]=1
      elif num==1:
            num1=num1+1
            level2[elem]=1;
            if elem in extra:
                 extracount=extracount+1
      elif num==0:
            num0=num0+1
  
file.close()

print "nums "+str(num2)+" "+str(num1)+" "+str(num0)
print extracount
print "Mips Athaliana annotation number:"+str(len(namehash))
print "Athaliana Level2 annotation number:"+str(len(level2.keys()))
print "Athaliana Level3 annotation number:"+str(len(level3.keys())) 


filename="spombe_funcat20062005."
file=open(filename)
level2={}
level3={}
namehash={}
for line in file:
      parts=line.split("|")
      namehash[parts[0]]=1
      elem=parts[1]
      num=elem.count(".")
      if num>=2:
            levelparts=elem.split(".")
            level2[levelparts[0]+"."+levelparts[1]]=1
            level3[levelparts[0]+"."+levelparts[1]+"."+levelparts[2]]=1
      elif num==1:
             level2[elem]=1;
             #level3
      #elif num==0:
              
file.close()

print "Mips Pombe annotation number:"+str(len(namehash))
print "Pombe Level2 annotation number:"+str(len(level2.keys()))
print "Pombe Level3 annotation number:"+str(len(level3.keys())) 

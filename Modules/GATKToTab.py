__author__ = 'Pier-Luc Plante'
__license__ = 'GPL'

"""
This program is created to be part of SNPs detection pipeline.
Usage: python BamToTab.py Reference.fasta GATK.out 
"""

import sys

#Minimum value of coverage to pass from uncertain to Observed or Variant
cut_Off_Uncertainty=5
#Class created to be part of the dictionnary. Only an assembly of int, str and char.
class LineOfTab(object):
    position = 0
    numA = 0
    numT = 0
    numG = 0
    numC = 0
    total = 0
    ref = ''
    now = ''
    type = ""

#TODO: create functions


if len(sys.argv)==0:
	print(__doc__)
	sys.exit(1)

################################################################################
### first step: calculate lenght of ref to create a dict of the good lenght. ###
################################################################################


#print("START TO IMPORT REF SEQ")	#for debogue
longOfSeq=0
seqRef=""
for lines in open(sys.argv[1]):
	if lines[0]!='>':
		longOfSeq=len(lines)-1+longOfSeq
		seqRef=seqRef+lines[0:len(lines)-1]
i=0
tabDict={}
while i<longOfSeq:
	tabDict[i]=LineOfTab()
	tabDict[i].position=i+1
	tabDict[i].ref=seqRef[i]
	i+=1

seqRef=""		#"free" seqRef memory? does it really work?

#print("Ref sequence is now in memory!")	#for debogue
################################################################################
#### Second step: Read the input from GATK file and put it in the dictionnary ###
################################################################################

gatk=open(sys.argv[2],'r');
line=gatk.readline();
if line.split()[0]!="Locus":
	print("Check GATK file");
	sys.exit(1);
while (1):
	line=gatk.readline();
	if line=="":
		break;
	else:
                position=int(line.split()[0].split(':')[1])-1;
                total=line.split()[1];
                a=int(line.split()[4].split(':')[1]);
                c=int(line.split()[5].split(':')[1]);
                t=int(line.split()[7].split(':')[1]);
                g=int(line.split()[6].split(':')[1]);
		tabDict[position].numA=a;
		tabDict[position].numC=c;
		tabDict[position].numT=t;
		tabDict[position].numG=g;
		tabDict[position].total=total;
		

#################################################
# start of the dict completion (type and total) #
#################################################
i=0
#print("start of the dict completion")	#for debogue
while i < len(tabDict):
	tabDict[i].total=tabDict[i].numA+tabDict[i].numG+tabDict[i].numC+tabDict[i].numT
	#fill the total attribut
	#if tabDict[i].total!=0:
	#	tabDict[i].numA=float(tabDict[i].numA)/tabDict[i].total
	#	tabDict[i].numT=float(tabDict[i].numT)/tabDict[i].total
	#	tabDict[i].numG=float(tabDict[i].numG)/tabDict[i].total###
	#	tabDict[i].numC=float(tabDict[i].numC)/tabDict[i].total
	#fill the now attribut
	if tabDict[i].total==0:
		tabDict[i].now=tabDict[i].ref
	elif (tabDict[i].numA>tabDict[i].numT and tabDict[i].numA>tabDict[i].numG and 
	tabDict[i].numA>tabDict[i].numC and tabDict[i].total!=0):
		tabDict[i].now='A'
	elif (tabDict[i].numT>tabDict[i].numA and tabDict[i].numT>tabDict[i].numG and 
	tabDict[i].numT>tabDict[i].numC and tabDict[i].total!=0):
		tabDict[i].now='T'
	elif (tabDict[i].numC>tabDict[i].numT and tabDict[i].numC>tabDict[i].numG and 
	tabDict[i].numC>tabDict[i].numA and tabDict[i].total!=0):
		tabDict[i].now='C'
	elif (tabDict[i].numG>tabDict[i].numT and tabDict[i].numG>tabDict[i].numA and
	tabDict[i].numG>tabDict[i].numC and tabDict[i].total!=0):
		tabDict[i].now='G'
	else:
		tabDict[i].now='N'
	#fill the type attribut
	if tabDict[i].total==0:
		tabDict[i].type="Reference"
	elif tabDict[i].total < cut_Off_Uncertainty and tabDict[i].ref==tabDict[i].now:
		tabDict[i].type="Reference low coverage"
	elif tabDict[i].total < cut_Off_Uncertainty and tabDict[i].ref!=tabDict[i].now:
		tabDict[i].type="Variant low coverage"
	elif tabDict[i].total >= cut_Off_Uncertainty and tabDict[i].ref==tabDict[i].now:
		tabDict[i].type="Observed"
	elif tabDict[i].total >= cut_Off_Uncertainty and tabDict[i].ref!=tabDict[i].now:
		tabDict[i].type="Variant"
	else:
		tabDict[i].type="Erreur"
	i+=1
	#print("position"+str(i-1)+" is completed in the dict")	#for debogue

##############################
# end of the dict completion #
##############################

#################
# print results##
#################

#print("start to print the results")	#for debogue
#position=0		
#numA=0
#numT=0
#numG=0
#numC=0
#total=0
#ref=''
#now=''
#type=""
i=0
print("Pos\tA\tT\tC\tG\tTotal\tReference\tNow\tType")
while i<len(tabDict):
	print(str(tabDict[i].position)+"\t"+str(tabDict[i].numA)+"\t"+str(tabDict[i].numT)+
	"\t"+str(tabDict[i].numC)+"\t"+str(tabDict[i].numG)+"\t"+str(tabDict[i].total)+
	"\t"+str(tabDict[i].ref)+"\t"+str(tabDict[i].now)+"\t"+str(tabDict[i].type))
	i+=1



	
	
	
	


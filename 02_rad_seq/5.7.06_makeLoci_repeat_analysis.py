#!/usr/bin/python

#This script iterates through a RADseq vcf file (including non-variant sites) to define RAD-loci. A RAD-locus is a continuous stretch of called genotypes with gaps no large than a distance threshold (set distthresh). Whenever there is a larger gap, a new locus is created. Loci with an average depth larger than twice the overall average depth or smaller than 3 times the overall average depth and loci shorter than 250 bp are removed. The output is a vcf file which has the locus number is the ID column of each site. 

#Input order: path to vcf file, path to output vcf file

import sys
import numpy as np
import os
import gzip

distthresh=100 #this is the distance threshold between sites to make a new locus rather than extending the previous one

if sys.argv[1][-3::]==".gz":
	vcf=gzip.open(sys.argv[1],"r")
else:
	vcf=open(sys.argv[1],"r")

temp_out=open("loci_temp.vcf","w")

if sys.argv[1][-3::]==".gz":
	out=gzip.open(sys.argv[2],"w")
else:
	out=open(sys.argv[2],"w")


RADloc=1
vcf_dict={}


line=vcf.readline().split("\t")
while True:
	if line[0][0]=="#":
		if line[0][1]!="#":
			header=line
		out.write("\t".join(line))
		line=vcf.readline().split("\t")
		continue
	pos=int(line[1])
	if RADloc not in vcf_dict:
		vcf_dict[RADloc]=[pos,[]]
	line[2]=str(RADloc)
	dp=[x.split(":")[3] for x in line[9::]]
	while "." in dp:
		dp.remove(".")
	dp=[int(x) for x in dp]
	mean_dp=np.mean(dp)
	vcf_dict[RADloc][1].append(mean_dp)
	temp_out.write("\t".join(line))
	next_line=vcf.readline()
	if next_line=="":
		vcf_dict[RADloc].append(pos)
		break
	next_line=next_line.split("\t")
	if line[0]!= next_line[0] or int(next_line[1])-int(line[1])>distthresh:
		vcf_dict[RADloc].append(pos)
		#print vcf_dict[RADloc]
		RADloc+=1
	line=next_line

temp_out.close()
print "Number of loci counted:", RADloc

all_loc_dp=[]
loc_len=[]
all_max_dp=[]
for RADloc in vcf_dict:
	locus=vcf_dict[RADloc]
	loc_len.append(int(locus[-1]) - int(locus[0]))
	loc_dp=locus[1]
	mean_loc_dp=np.mean(loc_dp)
	max_loc_dp=max(loc_dp)
	all_loc_dp.append(mean_loc_dp)
	all_max_dp.append(max_loc_dp)
	vcf_dict[RADloc].append((loc_len[-1],mean_loc_dp,max_loc_dp))

all_mean_dp=np.mean(all_loc_dp)
all_mean_len=np.mean(loc_len)
all_mean_max=np.mean(all_max_dp)

print "Overall mean depth:", all_mean_dp
print "Overall mean max dp:", all_mean_max
print "Overall mean length:", all_mean_len

statout=open(sys.argv[2]+"_locistats.csv","w")
for length,depth in zip(loc_len,all_loc_dp):
	statout.write("{0},{1}\n".format(length,depth))

loc_len=[]
all_loc_dp=[]

for RADloc in vcf_dict.keys():
	loc_length,mean_loc_dp,max_loc_dp=vcf_dict[RADloc][-1]
	if mean_loc_dp>all_mean_dp*2 or mean_loc_dp<all_mean_dp/3 or loc_length<150 or max_loc_dp>all_mean_max*2: #here the locus filters can be adjusted
		del vcf_dict[RADloc]
		continue
	loc_len.append(loc_length)
	all_loc_dp.append(mean_loc_dp)

good_loci=vcf_dict.keys()
temp_vcf=open("loci_temp.vcf","r")

for line in temp_vcf:
	linesp=line.split("\t")
	RADloc=int(linesp[2])
	if RADloc in good_loci:
		out.write(line)
out.close()

os.system("rm loci_temp.vcf")




print "Number of loci remain after depth and length filtering:", len(vcf_dict)
print "Mean length of locus after depth and length filtering:", np.mean(loc_len)
print "Mean depth of locus after depth and length filtering:", np.mean(all_loc_dp)

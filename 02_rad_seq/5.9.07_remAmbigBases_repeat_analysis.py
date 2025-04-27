#!/usr/bin/python

#This script removes any site from a VCF file where the reference allele is ambiguous (i.e. not A,T,C or G).

#Input order: path to input vcf file, path to output vcf file

import sys
import gzip



if sys.argv[1][-3::]==".gz":
	zipped=True
else:
	zipped=False

print "Zipped?", zipped
if zipped:
	f=gzip.open(sys.argv[1],"r")
	out=gzip.open(sys.argv[2],"w")
else:
	f=open(sys.argv[1],"r")
	out=open(sys.argv[2],"w")

if zipped:
	out=gzip.open(sys.argv[2],"w")
else:
	out=open(sys.argv[2],"w")
nucs=("A","T","C","G")

for line in f:
	if line[0]=="#":
		out.write(line)
	else:
		cols=line.split("\t")
		ref=cols[3]
		if ref in nucs:
			out.write(line)

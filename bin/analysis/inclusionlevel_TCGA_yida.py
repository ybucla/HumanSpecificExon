#! /u/home/y/ybwang/python

import os, re, sys
import numpy
from collections import defaultdict

def main():
	TCGAdir = '/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/TCGA_data/data/'
	junction = {};
	with open(TCGAdir + sys.argv[1]+'/input.txt','r') as f:
		for line in f:
			ele = line.rstrip().split("\t")
			gid = ele[1].replace("\"",'')
			(chr, strand) = ele[3:5]
			# whether there exists junction counts
			(inclusionCounts, skippingCounts) = ele[11:13]
			inclusionCounts = inclusionCounts.replace('NA','0')
			skippingCounts = skippingCounts.replace('NA','0')
			arr_inclu = [int(x) for x in inclusionCounts.split(",")]
			arr_skipp = [int(x) for x in skippingCounts.split(",")]
			if numpy.sum(arr_inclu)/2 + numpy.sum(arr_skipp) == 0:
				continue
			inclusion_level = numpy.sum(arr_inclu)/2/float(numpy.sum(arr_inclu)/2 + numpy.sum(arr_skipp))
			print chr+"\t"+ele[5]+"\t"+ele[6]+"\t"+gid+"\t0\t"+strand+"\t"+str(numpy.sum(arr_inclu))+"\t"+str(numpy.sum(arr_skipp)) + "\t"+str(inclusion_level)
			

if __name__ == '__main__':
	main()

#! /u/home/y/ybwang/python

from collections import defaultdict
from itertools import combinations_with_replacement
import glob
import os
import re
import sys

def main():
	dir = sys.argv[1]
	filelist = glob.glob(dir+"/*.localFDR.0.05.list")
	for file in filelist:
		print '>'+file
		source = file.replace('.localFDR.0.05.list','')
		exondict = defaultdict(set)
		pepdict = defaultdict(set)
		head = ''
		with open(file,'r') as f:
			for line in f:
				line = line.rstrip()
				if re.search("======(\w+)",line) is not None:
					head = source+'.'+re.search("======(\w+)",line).group(1)
				elif '#' not in line:
					ele = line.split("\t")
					if ele[4] == 'Y':
						pep = ele[6][int(ele[1]):int(ele[2])+1]
						arr = ele[7].split(';')
						for i in range(len(arr)-1):
							exondict[arr[i]].add(head)
							pepdict[arr[i]].add(pep)
		for e in exondict:
			print e+"\t"+';'.join(exondict[e])+"\t"+';'.join(pepdict[e]) + "\t" + str(len(exondict[e]))
if __name__ == '__main__':
		main()

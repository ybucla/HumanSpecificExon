#!/usr/bin/python
from collections import defaultdict
import sys, glob, re
import os.path

def main():
	tissue = sys.argv[1]

	dir = '/u/nobackup/yxing/PROJECT/ybwang/HumanSpecificExon/asevent_Mathias_Wilhelm/'

	exonDict = defaultdict(dict)
	countDict = defaultdict(dict)
	tissuedir = glob.glob(dir+tissue+'/*')
	for d in tissuedir:
		basename = os.path.basename(d)
		file_fromGTF_SE = d + '/fromGTF.SE.txt'
		file_JC_raw_input_SE = d + '/JC.raw.input.SE.txt'
		with open(file_fromGTF_SE, 'r') as f:
			for line in f:
				if re.match('^ID', line) is not None:
					continue
				ele = line.rstrip().split("\t")
				exonDict[basename][ele[0]] = "\t".join(ele[1::])
		with open(file_JC_raw_input_SE, 'r') as f:
			for line in f:
				if re.match('^ID', line) is not None:
					continue
				ele = line.rstrip().split("\t")
				(icount,scount) = (int(ele[1]), int(ele[2]))
				ilevel = 'NA' if icount/2+scount ==0 else str(icount/2/float(icount/2+scount))
				countDict[basename][ele[0]] = ele[1]+"\t"+ele[2]+"\t"+ilevel
	result = defaultdict(dict)
	for basename in exonDict:
		for i in exonDict[basename]:
			ele = exonDict[basename][i].split("\t")
			ilevel = countDict[basename][i]
			eEvent = ele[2]+"\t"+ele[4]+"\t"+ele[5]+"\t"+ele[0].replace("\"",'')+"\t0\t"+ele[3]+"\t"+("\t".join(ele[6:10]))
			result[eEvent][basename]=ilevel
	for eEvent in result:
		(mean, num, inc, ski, ninc, nski) = (0,0,0,0,0,0)
		for basename in result[eEvent]:
			ele = result[eEvent][basename].split("\t")
			if ele[2] != 'NA':
				num = num+1
				mean = mean + float(ele[2])
			if ele[0] != 'NA':
				ninc = ninc + 1
				inc = inc + int(ele[0])
			if ele[1] != 'NA':
				nski = nski + 1
				ski = ski + int(ele[1])
			print eEvent+"\t"+result[eEvent][basename]+"\t"+basename
		mean = 'NA' if num == 0 else str(mean/float(num))
		inc = 'NA' if ninc == 0 else str(inc/float(ninc))
		ski = 'NA' if nski == 0 else str(ski/float(nski))
		print '#'+eEvent+"\t"+inc+"\t"+ski+"\t"+mean+"\taverage\n"
			
						
#--------------------------#
if __name__ == '__main__':
	main()

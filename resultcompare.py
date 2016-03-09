#!/usr/bin/python
import os

def main():
	percolatorResult()


# suf function
def percolatorResult():
	cometout = 'cometout/'
	dirlist = os.listdir(cometout)
	data = {}
	for i in dirlist:
		if i != 'GM18486.rna' and i != 'GM18909.rna':
			continue
		percolatorPath = cometout + i + '/0.05/*.0.05.pep.percolator'
		identityPath = 'massDb/identity/identify_' + i + '.fa'
		reg = 'chr'
		for t in os.popen('cut -f 2 ' + identityPath).readlines():
			reg = reg + '|' + t.strip()
		result = os.popen("awk -F '\\t' '$3 < 0.05' " + percolatorPath + " | awk -F '\\t' '$6 ~ /" + reg + "/' | cut -f 6 | sort | uniq").readlines()
		for t in result:
			if data.has_key(t):
				data[t].append(i)
			else:
				data[t] = [i]
	print data
		

if __name__=='__main__':
        main()

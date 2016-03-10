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
		#if i != 'GM18486.rna' and i != 'GM18909.rna':
			#continue
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
	num = {}
	for t in data:
		ele = data[t]
		for i in range(len(ele)):
			for j in range(i,len(ele)):
				key1 = ele[i] + "\t" + ele[j]
				key2 = ele[j] + "\t" + ele[i]
				if key1 in num:
					num[key1] = num[key1] + 1
					num[key2] = num[key1]
				else:
					num[key1] = 1
					num[key2] = num[key1]
		
	print "\t","\t".join(dirlist)
	for i in dirlist:
		print i,
		for j in dirlist:
			key = i + "\t" + j
			n = 0
			if key in num:
				n = num[key]
			print "\t",n,
		print ""


if __name__=='__main__':
        main()

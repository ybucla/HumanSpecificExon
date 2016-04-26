#!/usr/bin/python
from optparse import OptionParser
import os
from collections import defaultdict
import re
import subprocess

def main():
	usage = 'usage: %prog <options> star_output_folder'
	parser = OptionParser(usage)
	
	parser.add_option('-o', dest='output_file', default='putative_junc.fa', help='Output result filename [Default %default]')
	parser.add_option('-l', dest='flank', type='int', default=66, help='Extend flanking junction ends by this number of bp [Default %default]')
	parser.add_option('-g', dest='genome_file', default='/u/home/f/frankwoe/nobackup/hg19/hg19_by_chrom/', help='genomic fasta directory (by chromosomes) [Default %default]')
	parser.add_option('--min-junc-reads', dest='min_junc_reads', default=2, type='int', help='Minimum number of reads required spanning the junction [Default %default]')
	parser.add_option('--trim-RK', dest='trim_RK', default=False, action='store_true', help='Indicate whether trim both ends to first R or K [Default %default]')
	parser.add_option('--verbose', dest='verbose', default=False, action='store_true', help='Verbose mode -- print DNA to stdin [Default %default]')
	
	(options, args) = parser.parse_args()
	
	if len(args) < 1:
		parser.error('Missing required input.')
	
	args[0]=os.path.abspath(args[0])
	verbose = options.verbose
	
	tag = 0
	head = ''
	idsDict = defaultdict(list)
	with open(args[0], 'r') as f:		
		for line in f:
			line = line.rstrip()
			if line.find('======') != -1:
				tag = 1
				head = line.replace('======','')
			if tag == 1 and line.find('chr') == 0:
				if head == args[2]:
					idsDict[head].append(line)
	for key in idsDict:
		ids = defaultdict(list)
		for line in idsDict[key]:
			ele = line.split('_')
                        k = str(ele[0]) + '_' + str(ele[1]) + '_' + str(ele[2])
                        ids[k].append(line)
	star_path = args[1] # 'RNA/' + key + '/Aligned.out.sorted.bam'
	junctionReads(star_path, ids, key)
	#ids = massresult(args[1])
	#junction_reads = junctionReads(args[0] + '/Aligned.out.sorted.bam', ids)

def read_junctions(file, reads, min_reads):  # [chr, id, strand, left_junction, right_junction]
	n = 0;
	junctions=[]
	is_annotated = lambda x: False if x=='0' else True
	with open(file, 'r') as f:
		for line in f:
			ele = line.split('\t')
			key = '_'.join([ele[0], str(int(ele[1])-1), str(int(ele[2])+1) ])
			#if is_annotated(ele[5]) or len(reads[key]) < min_reads:
			tmpdict = {}
			for tmp in reads[key]:
				tmpdict[tmp[0]] = ''
			if len(tmpdict) < min_reads:
				continue
			n = n + 1
			read_ids, left, right = zip(*reads[key])
			left_max = max(left)
			right_max = max(right)
			if ele[3]=='1':
				junctions.append([ ele[0], key,  '+', [int(ele[1]) - left_max +1, int(ele[1]) -1 ], [int(ele[2]) +1, int(ele[2]) + right_max-1 ] ])
			elif ele[3]=='2':
				junctions.append([ ele[0], key,  '-', [int(ele[1]) - left_max +1, int(ele[1]) -1], [int(ele[2]) +1, int(ele[2]) + right_max-1 ] ])
		print "# junction more than 6 reads:\t" + str(n) + "\t" + str(len(junctions))
	return junctions


# Shayna's function for junctions from STAR	

def junctionReads(sam_fn, msids, title):  # bam file input
# store reads id that mapped to human_epcific_exon to a dict
	print '>',title
	exonid = defaultdict(list)
	d = subprocess.Popen('samtools view -b -q 255 ' + sam_fn + ' | bedtools intersect -a /u/home/y/ybwang/scratch/HumanSpecificExon/data/hg19_repeatMasker_Alu.sorted.bed -b stdin -wa -wb -split -sorted', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in d.stdout.readlines():
		line = line.rstrip()
		ele = line.split()
		rs = int(ele[7]) + 1
		readsid = str(rs)+'_'+ele[9].split("/")[0] # format xMyNzMa
		exonkey = str(ele[0]) + '_' + str(ele[1]) + '_' + str(ele[2]) + '_' + str(ele[3]) + '_' + str(ele[4]) + '_' + str(ele[5])
		if ele[3] == 'test':
			continue
		exonid[readsid].append(exonkey)
	retval = d.wait()
## get reads that span each junction
	junctionReads = defaultdict(list) # junctionReads[chr_start_end] =[list of read ids that span junction]
	p = subprocess.Popen('samtools view -q 255 ' + sam_fn, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		line = line.rstrip()
		if ((not line.startswith('@')) and ('N' in line.split()[5])): #not a header line and read spans a splice junction
			ele = line.split()
			bamid =  ele[3]+'_'+ele[0]
			if (bamid not in exonid):
				continue
			cigar = ele[5] # format xMyNzM
			exon_lengths = re.findall(r"(\d+M)",cigar) #[list of exon lengths]
			intron_lengths = re.findall(r"(\d+N)",cigar) #[list of intron lengths]
			pos = int(ele[3]) #read start location
# iterate through exons/introns to extract junction coordinates that read spans
			for index,item in enumerate(intron_lengths):
				jS = pos + int(exon_lengths[index][:-1]) - 1  # Zj: should add -1
				jE = jS + int(intron_lengths[index][:-1]) + 1 # Zj: should add +1 instead of -1
				key = str(ele[2]) + '_' + str(jS) + '_' + str(jE) #chr_start_end
				for t in exonid[bamid]:
					if key in msids:
						for k in msids[key]:
							print k,"\t",key,"\t",t,"\t",bamid
				junctionReads[key].append([ ele[0], int(exon_lengths[index][:-1]), int(exon_lengths[index+1][:-1]) ] )  # chr_start_end => [id, left_length, right_length]
				pos = jE
	retval = p.wait()
# end of looping through sam file
	return junctionReads	

def massresult(idfile):
	ids = defaultdict(list)
	with open(idfile, 'r') as f:
                for line in f:
			line = line.rstrip()
			ele = line.split('_')
			key = str(ele[0]) + '_' + str(ele[1]) + '_' + str(ele[2])
			ids[key].append(line)
	return ids	

if __name__=='__main__':
	main()

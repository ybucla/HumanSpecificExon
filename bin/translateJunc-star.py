#!/usr/bin/python
from optparse import OptionParser
import os
from collections import defaultdict
import re
import subprocess
import random

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
	junction_reads = junctionReads(args[0] + '/Aligned.out.sorted.bam')
	junctions_all = read_junctions(args[0] + '/SJ.out.tab', junction_reads, options.min_junc_reads)
	junctions = junctionFilter(junctions_all, args[0])
	junctions.append(['final',''])
	#print len(junctions)
	sj_dict = {}
	with open(args[1], 'r') as f:
		for line in f:
			line = line.rstrip()
			ele = line.split()			
			key = ele[0]+'_'+str(int(ele[1])-1)+'_'+str(int(ele[2])+1)+'_'+ele[3]
			sj_dict[key] = line
	last_chr=junctions[0][0]
	chr_juncs=[]
	with open(options.output_file, 'w') as out:
		np = 0;
		for junc in junctions:
			chr = junc[0]
			if chr == last_chr:
				chr_juncs.append(junc[1:])
				continue
			else:
				seqs = get_seqs(last_chr, chr_juncs, options.flank, options.genome_file,sj_dict)
				for i in range(len(chr_juncs)):
					id, strand, junc_left, junc_right = chr_juncs[i]
					seq = seqs[i]
					id_line = '>' + id + '_' + last_chr + '_' + strand + '_' + ','.join([str(x) for x in junc_left]) + '_' + ','.join([str(x) for x in junc_right])
					if verbose:
						print id_line
						print seq
					prot_seqs = translate_dna(seq)
					for j in range(len(prot_seqs)):
						prot_seq = prot_seqs[j]
						tr_prot_seq = trim_prot(prot_seq, options.trim_RK)
						if len(tr_prot_seq) < 5 or tr_prot_seq.upper() == tr_prot_seq:
							np = np + 1
							continue
						out.write(id_line + '_ORF:' + str(j) + '\n' + tr_prot_seq + '\n')
				chr_juncs = [junc[1:]]
				last_chr = chr
		print "# stop and 5 aa less pep:\t"+str(np)

def trim_prot(prot_seq, trim_RK = True):
	#print "protein len:\t", prot_seq,"\t",len(prot_seq),"\t",
	splice_site = ''.join(x for x in prot_seq if x.islower())  # must contain splice site
	
	if len(splice_site) < 1:
		return ''
	
	left, right = prot_seq.split(splice_site)
	left_start = len(left) - left[::-1].find('_')
	right_end = right.find('_')
	
	if left_start > len(left):  # if left has no stop codon
		left_start = 0
	
	if right_end < 0:  # if right has no stop codon
		right_end = len(right)
	
	prot_seq = prot_seq[left_start : len(left) + right_end + 1]
	#print prot_seq,"\t",len(prot_seq)
	if not trim_RK:
		return prot_seq
	
	start, end = 0, 0
	for ss in range(len(prot_seq)):
		if prot_seq[ss].upper()=='R' or prot_seq[ss].upper()=='K':
			if not '_' in prot_seq[ss:]:
				start = ss
				break
		if prot_seq[ss].upper() != prot_seq[ss]:
			start = ss
			break
	
	rev_prot = prot_seq[start:][::-1]
	for ee in range(len(rev_prot)):
		if rev_prot[ee].upper()=='R' or rev_prot[ee].upper()=='K':
			if not '_' in rev_prot[ee:]:
				end = ee
				break
		if rev_prot[ee].upper() != rev_prot[ee]:
			end = ee - 1
			break
	return rev_prot[ee:][::-1]
	
def translate_dna(seq):
	codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
	prot_seqs = []
	for start_site in range(3):
		cds = seq[start_site:]
		prot_seq = ''
		for n in range(0, len(cds), 3):
			current_triple = cds[n:n+3]
			if current_triple.upper() in codontable:
				letter = codontable[current_triple.upper()]
				if current_triple == current_triple.upper():
					prot_seq += letter
				else:
					prot_seq += letter.lower()
			else:
				#prot_seq += '-'
				continue
		prot_seqs.append(prot_seq)
	
	return prot_seqs

def get_seqs(chr, junctions, flank, file_dir,sj_dict):
	seqs = []
	genomic_seq = read_genome(file_dir, chr)
	for id, strand, junc_left, junc_right in junctions:
		tag_strand = '1'
		if strand == '-':
			tag_strand = '2'
		key = chr+'_'+str(junc_left[1])+'_'+str(junc_right[0])+'_'+tag_strand
		ele = sj_dict[key].split()
		left = flank + abs(junc_left[1] - junc_left[0]) + 1
		right = flank + abs(junc_right[1] - junc_right[0]) + 1
		if int(ele[11]) >0 and int(ele[11]) < left:			
			left_seq = [genomic_seq[x] for x in range(junc_left[1] - int(ele[11]), junc_left[1])]
			junc_left.insert(0, junc_left[1] - int(ele[11]) + 1)
		else:
			left_seq = [genomic_seq[x] for x in range(junc_left[0] - flank - 1, junc_left[1])]
			junc_left.insert(0, junc_left[1] - flank)
		if int(ele[12]) >0 and int(ele[12]) < right:	
			right_seq = [genomic_seq[x] for x in range(junc_right[0] - 1, junc_right[0] + int(ele[12]) - 1)]
			junc_right.append(junc_right[0] + int(ele[12]) - 1)
		else:
			right_seq = [genomic_seq[x] for x in range(junc_right[0] - 1, junc_right[1] + flank)]
			junc_right.append(junc_right[1] + flank)
		seq = [x.upper() for x in left_seq] + [right_seq[0].lower()] + [x.upper() for x in right_seq[1:]]
		if strand == '-':
			seq = rev_complement(seq)
		seq_str = ''.join(seq)
		#seq_str = seq_str.upper()
		seqs.append(seq_str)
	return seqs

def rev_complement(seq):
	rev_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', '|':'|'}
	rev_seq = seq[::-1]
	rev_comp_seq = map(lambda x: rev_dict[x], rev_seq)
	return rev_comp_seq

def read_genome(file_dir, chr):
	genome=[]
	file_path = os.path.abspath('%s/%s.fa' % (file_dir, chr))
	with open(file_path, 'r') as f:
		for line in f:
			line = line.rstrip()
			if line[0] == '>':
				chr=line[1:]
				continue
			else:
				genome.extend(line)
	return genome


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

def junctionReads(sam_fn):  # bam file input
# store reads id that mapped to human_epcific_exon to a dict
	exonid = {}
	d = subprocess.Popen('samtools view -b -q 255 ' + sam_fn + ' | bedtools intersect -a /u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5.unique.sorted.bed -b stdin -wa -wb -split -sorted', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in d.stdout.readlines():
		line = line.rstrip()
		ele = line.split()
		rs = int(ele[7]) + 1
		readsid = str(rs)+'_'+ele[9].split("/")[0] # format xMyNzM
		if ele[3] == 'test':
			continue
		exonid[readsid] = ''
	retval = d.wait()
## print count
	print '# reads mapping number:\t' + str(len(exonid))
## get reads that span each junction
	exon_juncid = {}
	junctionReads = defaultdict(list) # junctionReads[chr_start_end] =[list of read ids that span junction]
	p = subprocess.Popen('samtools view -q 255 ' + sam_fn, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	ntmp = 0
	for line in p.stdout.readlines():
		line = line.rstrip()
		if ((not line.startswith('@')) and ('N' in line.split()[5])): #not a header line and read spans a splice junction
			ele = line.split()
			bamid = ele[3]+'_'+ele[0]
			if (bamid not in exonid):
				continue
			ntmp += 1
			exon_juncid[bamid] = exonid[bamid]
			cigar = ele[5] # format xMyNzM
			exon_lengths = re.findall(r"(\d+M)",cigar) #[list of exon lengths]
			intron_lengths = re.findall(r"(\d+N)",cigar) #[list of intron lengths]
			pos = int(ele[3]) #read start location
# iterate through exons/introns to extract junction coordinates that read spans
			for index,item in enumerate(intron_lengths):
				jS = pos + int(exon_lengths[index][:-1]) - 1  # Zj: should add -1
				jE = jS + int(intron_lengths[index][:-1]) + 1 # Zj: should add +1 instead of -1
				key = str(ele[2]) + '_' + str(jS) + '_' + str(jE) #chr_start_end
				junctionReads[key].append([ ele[0], int(exon_lengths[index][:-1]), int(exon_lengths[index+1][:-1]) ] )  # chr_start_end => [id, left_length, right_length]
				pos = jE
	retval = p.wait()
# print count
	print '# junction reads mapping number:\t'+str(len(exon_juncid)) + "\t" + str(ntmp)
	print '# junction number:\t'+str(len(junctionReads))
# end of looping through sam file
	return junctionReads	

def junctionFilter(junctionReads, out):
	junctionReads_filter = []
	tmpbed = 'x' + str(random.randint(1000000, 9999999))+'-' + out.split('/')[-1]
	with open(tmpbed, 'w') as out:
		for i in junctionReads:
			id = i[0]+";"+i[1]+";"+i[2]+";"+str(i[3][0])+"_"+str(i[3][1])+"_"+str(i[4][0])+"_"+str(i[4][1])
			line_left = i[0]+"\t"+str( i[3][0] - 1 )+"\t"+str(i[3][1])+"\t"+id
			line_right = i[0]+"\t"+str( i[4][0] - 1 )+"\t"+str(i[4][1])+"\t"+id
			out.write(line_left+"\t0\t"+i[2]+"\n")
			out.write(line_right+"\t0\t"+i[2]+"\n")
	p = subprocess.Popen('bedtools coverage -a ' + tmpbed + ' -b /u/home/y/ybwang/nobackup-yxing-PROJECT/HumanSpecificExon/data/Ensembl_Alu_25bp_0.5.unique.sorted.bed -s', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	tmphash = defaultdict(list)
	for line in p.stdout.readlines():
                line = line.rstrip()
		ele = line.split()
		if int(ele[6]) > 0:
			if ele[3] not in tmphash:
				arr = ele[3].split(';')
				sites = arr[3].split('_')
				junctionReads_filter.append([arr[0],arr[1],arr[2],[int(sites[0]),int(sites[1])],[int(sites[2]),int(sites[3])]])
				tmphash[ele[3]] = ''
	retval = p.wait()
	os.remove(tmpbed)
	return junctionReads_filter
	
	
if __name__=='__main__':
	main()

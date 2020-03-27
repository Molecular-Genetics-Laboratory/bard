# --------------------------------------------------------------------------
# bard: Batch Analyzer of Ribo-seq Data
#
# Version 1.0
#
# Copyright (C) 2019, 2020  Somdeb Chattopadhyay (somdeb.ch@gmail.com)
#
# Molecular Genetics Laboratory
# National Institute of Immunology, New Delhi
# http://www.nii.res.in

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, see <http://www.gnu.org/licenses/>
#
#
#
# Usage:
# ~$ bard <JSON config file>
#
#---------------------------------------------------------------------------

from collections import defaultdict
from datetime import datetime as dt
import pandas as pd
from Bio import SeqIO, SeqUtils
import seaborn as sns
import numpy as np
import random
import pysam
import json
import sys
import os


import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
colors = plt.rcParams["axes.prop_cycle"]()

np.random.seed(12345)

# This line is automatically updated before each commit
versionstr = "bard v1.0 ID=14-29-15-27-03-2020"

codon_keys = {
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

codon_counts = {
	'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
	'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
	'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
	'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
	'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
	'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
	'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
	'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
	'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
	'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
	'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
	'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
	'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0
}

# List of pause score values per codon
codon_pauses = {
    'ATA':[], 'ATC':[], 'ATT':[], 'ATG':[],
    'ACA':[], 'ACC':[], 'ACG':[], 'ACT':[],
    'AAC':[], 'AAT':[], 'AAA':[], 'AAG':[],
    'AGC':[], 'AGT':[], 'AGA':[], 'AGG':[],
    'CTA':[], 'CTC':[], 'CTG':[], 'CTT':[],
    'CCA':[], 'CCC':[], 'CCG':[], 'CCT':[],
    'CAC':[], 'CAT':[], 'CAA':[], 'CAG':[],
    'CGA':[], 'CGC':[], 'CGG':[], 'CGT':[],
    'GTA':[], 'GTC':[], 'GTG':[], 'GTT':[],
    'GCA':[], 'GCC':[], 'GCG':[], 'GCT':[],
    'GAC':[], 'GAT':[], 'GAA':[], 'GAG':[],
    'GGA':[], 'GGC':[], 'GGG':[], 'GGT':[],
    'TCA':[], 'TCC':[], 'TCG':[], 'TCT':[],
    'TTC':[], 'TTT':[], 'TTA':[], 'TTG':[],
    'TAC':[], 'TAT':[], 'TAA':[], 'TAG':[],
    'TGC':[], 'TGT':[], 'TGA':[], 'TGG':[],
}

aa_to_codon = {
    'I' : ['ATA', 'ATC', 'ATT'],
    'M' : ['ATG'],
    'T' : ['ACA', 'ACC', 'ACG', 'ACT'],
    'N' : ['AAC', 'AAT'],
    'K' : ['AAA', 'AAG'],
    'S' : ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'R' : ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'L' : ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'P' : ['CCA', 'CCC', 'CCG', 'CCT'],
    'H' : ['CAC', 'CAT'],
    'Q' : ['CAA', 'CAG'],
    'V' : ['GTA', 'GTC', 'GTG', 'GTT'],
    'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
    'D' : ['GAC', 'GAT'],
    'E' : ['GAA', 'GAG'],
    'G' : ['GGA', 'GGC', 'GGG', 'GGT'],
    'F' : ['TTC', 'TTT'],
    'Y' : ['TAC', 'TAT'],
    'C' : ['TGC', 'TGT'],
    'W' : ['TGG'],
    '*' : ['TAA', 'TAG', 'TGA']
}


# Holds the global configuration
global_config       = {}
# GFF
annotation          = {}
# Keep a track of the warnings emitted
__warnings          = []
# A nested hash mapping the gene names
# to their endmap vectors, for each
# length of RPF. These are the raw index
# counts
per_gene_endmaps    = {}
# Metagene vectors per readlength
metagene_for_len    = {}
# 5' and 3' overlapping gene names
overlap_dict        = {}
# Coverage per gene
gene_coverages      = {}
# Genes with coverage above threshold
high_coverage_genes = []
# Offsets
offsets_e           = {}
offsets_p           = {}
offsets_a           = {}
# ribosome density over CDS
per_gene_densities  = {}
# Dictionary mapping the gene name to the corresponding
# cDNA sequence
transcripts_dict    = {} # name --> sequence
# Nested dictionary mapping a gene name to it's sequence and
# the corresponding positional pause vector
pps_vector_per_gene = {}

# Coverage profile
coverage_profile    = {}


# Template configuration file, printed out when
# no input is provided.

CONFIG_TMP = """{
    "coding_sequence_path": "",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "",
    "annotation_feature_tag": "",
    "bam_file_path": "",
    "readcov_from_terminal": "",
    "read_offset_terminal": "",
    "readcov_percentage_length": 100,
    "readcov_ignore_nts": 20,
    "endmap_upstream": 50,
    "endmap_downstream": 200,
    "offset_for_site": "P",
    "mapping_buffer": 25,
    "filtering_threshold": 0,
    "will_ignore_overlaps": true,
    "peak_scan_range": [
        -20,
        -10
    ],
    "use_readlengths": [
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33
    ],
    "genes_overlap_exception": "",
    "gene_list_file": "",
    "gene_list_action": ""
}
"""

CONFIG_HELP = """
This file details the configuration options.

1. coding_sequence_path:        absolute path to the CDS/cDNA file

2. coding_sequence_format:      the format of the above file. eg: fasta 

3. annotation_file_path:        absolute path to the annotation file (GTF/GFF)

4. annotation_feature_tag:      the tag in the GFF column 9 which uniquely identifies
                                the feature (eg: gene_id or ID or protein_id)

5. bam_file_path:               absolute path to the alignment file (BAM)

6. readcov_from_terminal:       "five_prime" default

7. read_offset_terminal:        "five_prime" or "three_prime". Which read terminal to
                                use to calculate offsets

8. readcov_percentage_length:   int: default 100.

9. readcov_ignore_nts:          #nts to ignore from the start or stop codon of the
                                gene whilst computing the read coverage per transcripts

10. endmap_upstream:            #nt upstream of CDS start to include in the metagene profile

11. endmap_downstream:          #nt downstream of CDS start to include in metagene profile

12. mapping_buffer:             default 25. #nt to compansate for the offset shift

13. offset_for_site:            "E", "P", "A"

14. will_ignore_overlaps:       Ignore overlapping genes? 

15. peak_scan_range:            [-20, -5] or [1, 30] or whatever depending on offsets selected

16. use_readlengths:            List of read lengths to include for analysis

17. gene_list_file:             File containing a list of genes

18. gene_list_action:           "include_only", "exclude_only", "exclude_balance"

19. genes_overlap_exception:    File containing list of overlapping genes which we will NOT ignore
"""

unique = lambda vector: list(set(vector))


def load_json(file="object.json"):

    with open(file, "r") as read_file:
        return(json.load(read_file))


def save_json(data, file="object.json"):

    with open(file, "w") as write_file:
        json.dump(data, write_file, indent=4)


def extended_sort(x, by, reverse=False):
	'''
	Takes two lists and sorts them based on the
	values in one of them (by)
	'''
	by, x = zip(*sorted(zip(by, x), reverse=reverse))
	return (by, x)


def randseq(length, gcbias=0.7):
	'''
	Generate a random DNA sequence
	'''
	nt = ['A', 'T', 'G', 'C']
	gc_prob = (gcbias)/2
	at_prob = (1-gcbias)/2
	sq = np.random.choice(
		nt, 
		size = length, replace = True, 
		p=[at_prob,at_prob,gc_prob,gc_prob]
	)
	return("".join(sq.tolist()))


def expand_to_codon(v):
    '''
    Expand amino acid letters in a list to codons.
    eg: ["C", "ATG"] becomes ["TGT", "TGC", "ATG"]
    '''
    valid_codon = list(codon_keys.keys())
    valid_amino = list(codon_keys.values())
    tmp         = []
    def expanded(j):
        if j in valid_codon:
            return([j])
        if j in valid_amino:
            return(aa_to_codon[j])
    def invalid(v):
        if ((v in valid_codon) or
           (v in valid_amino)):
                return False
        else:
            return True
    for i in v:
        if isinstance(i, list):
            tmp_1 = []
            for j in i:
                if invalid(j): continue
                for c in expanded(j):
                    tmp_1.append(c)
            tmp.append(tmp_1)
        else:
            if invalid(i): continue
            for c in expanded(i):
                tmp.append(c)
    return(tmp)


def extract_file_name(path):
	"""
	Given an absolute path, extracts the file name
	(without extension)
	NOTE: May fail on Windows
	"""
	return(os.path.split(path)[1].split(".")[0])


def reset():
	'''
	Purge global environment
	'''
	global_config.clear()
	gene_coverages.clear()
	annotation.clear()
	__warnings.clear()
	per_gene_endmaps.clear()
	metagene_per_len.clear()
	overlap_dict.clear()
	per_gene_densities.clear()
	transcripts_dict.clear()
	pps_vector_per_gene.clear()
	offsets_p.clear()
	offsets_e.clear()
	offsets_a.clear()


def groupN(x,n):
	'''
	Group list (x) contents into n members each
	eg: x=[1,2,3,4,5,6], n=3, becomes [[1,2,3],[4,5,6]]
	'''
	if n <= 0:
		return
	tmp = []
	for i in range(0, len(x)+n,n):
		j = i + n
		vec = x[i:j]
		if len(vec) != 0:
			tmp.append(vec)
		if j == len(x): break
	return(tmp)


def ungroup(x):
	'''
	Ungroup a 1-level deep nested list
	'''
	tmp = []
	for i in x:
		for j in i:
			tmp.append(j)
	return(tmp)


def position_of_codon(string, motifs, mode="codon"):
	'''
	Returns the zero-indexed position of
	a codon for a given ORF sequence.
	'''
	if len(string) % 3 != 0:
		return "err"
	codons = groupN(string,3)
	indexes = []
	i = 0
	for motif in motifs:
		for codon in codons:
			if codon == motif:
				if mode == "nt":
					indexes.append(i*3)
				if mode == "codon":
					indexes.append(i)
			i += 1
		i = 0
	return(indexes)


def chunk_from_id(vec, idx, length):
	if (idx+length+1) > len(vec):
		raise RuntimeError
	if (idx - length) < 0:
		raise RuntimeError
	return(vec[idx-length:idx+length+1])


def random_hexadecimal_string():
	'''
	Random number prepended to output filenames
	to avoid conflicts
	'''
	return("%05x" % random.randrange(16*300))


def value_from_dict(dict, key):
	'''
	Return 0 if a dictionary lacks a key, do not crash
	'''
	try:
		return(dict[key])
	except KeyError:
		return 0


def notify(event, level="notf", onetime=False, fatal=False):
	'''
	Runtime logs
	'''

	if onetime: # issue a notification only once (in a loop)
		if event not in __warnings:
			__warnings.append(event)
			pass
		else:
			return

	timestamp = str(dt.now().strftime(
		"%H.%M.%S"
	))

	# Open in append mode
	log_path = open(global_config["log_file"], "a")


	log_path.write("[{}] {}\n".format(timestamp, event))
	log_path.close()

	RED = "\033[31m"
	GRN = "\033[32m"
	YLW = "\033[33m"
	RST = "\033[m"

	if level=="notf":
		print("{}[NOTF: {}] {}{}".format(
		GRN,timestamp, RST, event
	))

	if level=="warn":
		print("{}[WARN: {}] {}{}".format(
		YLW, timestamp, RST, event
	))

	if level=="crit":
		print("{}[CRIT: {}] {}{}".format(
		RED, timestamp, RST, event
 	))
		if fatal:
			# Can call fatal only at critical log level
			# global_config["bamfile"].close()
			raise SystemExit()



def set_global_config(**kwargs):
	'''
	Used by init. Will be removed.
	'''
	for key, value in kwargs.items():
		global_config[key] = value


def set_config(config):
	'''
	Set the global configuration from a dictionary
	'''
	for key, value in config.items():
		global_config[key] = value


def get_offset(xvec, yvec, right_bound, left_bound):
    '''
    Returns the P-site offset in a metagene vector
    '''

    if isinstance(xvec, np.ndarray):
            xvec = xvec.tolist()
    if isinstance(yvec, np.ndarray):
            yvec = yvec.tolist()

    left   = xvec.index(left_bound)
    right  = xvec.index(right_bound)
    y_trim = yvec[left:right+1] # because python
    pos    = yvec.index(max(y_trim))
    offset = pos - xvec.index(0)
    return(xvec[pos], offset)


def offset_delta_for_metagene():
	'''
	Get psite offsets from the metagene profiles of each RPF length
	'''

	left_bound   = global_config["peak_scan_range"][0]
	right_bound  = global_config["peak_scan_range"][1]

	x_axis  = metagene_for_len["xaxis"]

	for rpf_length, metagene_vector in metagene_for_len.items():

		if rpf_length == "xaxis":
			continue

		_,delta = get_offset (
			xvec        = x_axis,
			yvec        = metagene_vector,
			left_bound  = left_bound,
			right_bound = right_bound
		)

		offsets_p[rpf_length] = delta

		save_file(offsets_p, "psite_offsets")

		# Counting on the fact that the 5' offsets
		# returned above are negative
		offsets_e[rpf_length] = delta + 3 
		offsets_a[rpf_length] = delta - 3
		save_file(offsets_e, "esite_offsets")
		save_file(offsets_a, "asite_offsets")


def read_genes_from_file(filename):
	'''
	Given a file with one gene per line
	return the genes as a list
	'''
	try:
		fh = open(filename, "r")
	except:
		raise FileNotFoundError

	genes = []
	for gene in fh:
		genes.append(gene.strip("\n").strip(" "))
		
	return(genes)


def parse_gff():
	'''
	Converts a GFFv3/GTF file into a dictionary.
	Not fully GFFv3 compliant. Has no concept of exons or such.
	'''
	annofile    = global_config["annotation_file_path"]
	# tag in GFF column 9 which contains gene name
	feature_tag = global_config["annotation_feature_tag"]

	annolist = []
	for record in open(annofile, 'r'):
		record = record.strip("\n")
		record = record.split("\t")
		if ((len(record) < 9) and ("#" in record[0])):
		# comment line. ignore
			continue
		if len(record) != 9:
		# GFF/GTF should only have 9 columns
			continue
		annolist.append(record)

	all_gene_names = []
        # We can sometimes find duplicate gene names in the annotation files.
        # In case a duplicate is detected, this function appends an integer
        # suffix to it's name, which is a count of #times it was encountered
	def check_duplicate_name(gene_name):

		if gene_name in all_gene_names:
			notify("Found duplicate annotation. Will tag and ignore.",
				level="warn",
				onetime=True)
			all_gene_names.append(gene_name)

			occurances = str(all_gene_names.count(gene_name)-1)
			gene_name  = gene_name + "_" + occurances
		else:
			all_gene_names.append(gene_name)

		return (gene_name)

	def extract_gene_name(index):
		found_tag = False
		for attribute in annolist[index][8].split(";"):
			# strip the feature tag and return gene name
			if feature_tag in attribute:
				found_tag = True
				return (
					check_duplicate_name (
						attribute
						.replace(feature_tag, "")
						.replace("\"","")
						.replace("=", "")
						.strip()
					)
				)
			else:
				pass

		if not found_tag:
			# No point continuing
			notify("{}th gene lacks feature tag".format(i), level="crit",
			 fatal=True)

	# Populate the annotation dictionary
	for i in range(len(annolist)):
		annotation[extract_gene_name(i)] = {
				"start"  : int(annolist[i][3]),
				"stop"   : int(annolist[i][4]),
				"strand" : annolist[i][6],
				"refname": annolist[i][0]
		}

	save_file(annotation, "Annotations")


def getindex(vector, value):
	'''
	Return the 0-based index of a numpy array which
	contains the specific value
	'''
	try:
		return(np.where(vector == value)[0][0])
	except:
		return


def script_init():
	'''
	Run some basic setup operations
	'''

	session_id = random_hexadecimal_string()

	prefix = extract_file_name (
		global_config["bam_file_path"]
	)

	suffix = "_{}_{}site_{}_offsets_at_{}".format(
		prefix,
		global_config["offset_for_site"],
		global_config["read_offset_terminal"],
		str(dt.now().strftime("%I-%M-%S%p"))
	)

	# Creating a directory, WITHOUT checking for pre-existing ones
	# with the same name. This is exactly what we want in case of
	# re runs, and we are including the timestamp in the name.
	basename = "{}_Results_{}_{}site_{}_on_{}".format(
				session_id, # prevents any screwy conflicts
				prefix,
				global_config["offset_for_site"],
				global_config["read_offset_terminal"],
				str(dt.now().strftime("%d-%b-%Y_at_%I-%M-%S%p"))
			  )

	# We should never see this
	if os.path.isdir(basename):
		notify("Are you running this script in parallel?", level="crit")
		notify("Please ensure the BAM files are named differently.",
			level="crit", fatal=True)

	img_dir  = "{}/{}".format(basename, "Plots")
	data_dir = "{}/{}".format(basename, "Data")

	os.makedirs(img_dir)
	os.makedirs(data_dir)

	# Set up the file paths
	set_global_config (
			# Keep track of how many times the metagene
			# function has been called. This is because
			# the first call will be without offsets, and
			# we need to keep track of it so that we can name
			# the output files accordingly.
			metagene_callcount = 0,
			suffix             = suffix,
			img_dir            = img_dir  +"/", # this breaks Windows compatibility
			data_dir           = data_dir +"/",

			img_id             = "_{}_{}".format(
				prefix, str(dt.now().strftime("%I-%M-%S%p"))
			),

			log_file  = "{}/{}_logfile_for_{}_{}.txt".format(
				basename,
				session_id,
				prefix,
				str(dt.now().strftime("%I-%M-%S%p"))
			),

			conf_file  = "{}/{}_runconfig_{}_{}.json".format(
				basename,
				session_id,
				prefix,
				str(dt.now().strftime("%I-%M-%S%p"))
			),

			session_id     = session_id,
			prefix         = prefix,
			script_version = versionstr
	)

	notify("ScriptID: {}".format(global_config["script_version"]))

	try:
		global_config["bamfile"] = pysam.AlignmentFile(
			global_config["bam_file_path"], "rb")

		total_mapped_reads = int(
			pysam.flagstat(global_config["bam_file_path"]
		).split("\n")[4].split("+")[0].strip(" "))

		
		if total_mapped_reads < 1e6:
			notify("{} mapped reads (< 1e6). Results will be unreliable".format(
				total_mapped_reads), level="crit")


		global_config["mapped_reads"] = total_mapped_reads

		notify("Found BAM file, seems OK ...", level="notf")
	except:
		notify("BAM file is unavailable.", level="crit", fatal=True)

	# Check if BAM index file exists. This can be fatal.
	if not os.path.isfile (
		os.path.split(global_config["bam_file_path"])[0]+"/"+extract_file_name(
		global_config["bam_file_path"])+".bam.bai"
	):
		notify("BAM file index does not exist. Please create one.",
			level="crit", fatal=True)
	else:
		notify("Found BAM index (did not check contents)")

	try:
		thresh = global_config["filtering_threshold"]
	except:
		notify("Filter threshold not specified, assuming 70%")
		global_config["filtering_threshold"] = 70

	# We want to automatically detect this
	if global_config["peak_scan_range"][0] > global_config["peak_scan_range"][1]:
		notify("peak_scan_range: should be set as [lowerbound, upperbound]",
			level="crit")
		notify("Example: [-20, -10] or [0, 30] etc, where 0 denotes start", 
			level="crit", fatal=True)

	# Check if we need to include any overlapping genes
	try:
		overlap_exception = read_genes_from_file(global_config["genes_overlap_exception"])
		set_global_config(overlap_exception = overlap_exception)
	except:
		set_global_config(overlap_exception = [])


	parse_gff()
	identify_overlaps()


def check_gene_list():
	"""
	Check if a gene list has been provided. If so,
	should we include or exclude the genes from the analysis?
	"""
	try:
		listfile   = global_config["gene_list_file"]
		listaction = global_config["gene_list_action"]
	except:
		notify("Gene list file or associated action not mentioned.\
			Skipping.", level="warn")
		return

	def remove_annotation(gene):
		try:
			del annotation[gene]
		except:
			notify("{} not found in annotation".format(
				gene), level="warn")
			notify("Ignoring {}".format(gene), level="warn")

	# Should we perform the analysis using only
	# these genes? Or completely excluding these
	# genes?
	total_annotation_size = len(annotation)
	allgenes = annotation.copy().keys() #will throw error otherwise
	actions  = ["include_only", "exclude_only", "exclude_balance"]
	genelist = []

	try:
		genelist = read_genes_from_file(listfile)
	except:
		notify("No gene list provided. Moving on.", level="warn")
		return

	if listaction not in actions:
		notify("What to do with the gene list?", level="warn")
		notify("Try: 'include_only', 'exclude_only' or 'exclude_balance'",
			level="warn")
		notify("Ignoring provided gene list", level="warn")
		return

	#ft = global_config["filtering_threshold"]
	#if ft != 0:
	#	notify("Filter threshold set to {}%, probably set to 0?".format(ft),
	#		level="warn")

	if listaction == "exclude_only":
		# Remove these genes from annotation
		for gene in genelist:
			remove_annotation(gene)

	if listaction == "exclude_balance":
		# Remove these genes from annotation,
		# and keep an equal number of randomly
		# selected genes from remaining set
		for gene in genelist:
			remove_annotation(gene)

		reduced_genes = annotation.copy().keys()

		balance_sample = random.sample(reduced_genes,
			len(genelist)
		)

		for gene in allgenes:
			if gene not in balance_sample:
				remove_annotation(gene)

	if listaction == "include_only":
		# Remove all genes except these genes
		# from the annotation
		for gene in allgenes:
			if gene not in genelist:
				remove_annotation(gene)

	notify("Internal annotation changed from {} to {} genes".format(
		total_annotation_size, len(annotation)
	))


def save_file(data, filename):
	try:
		save_json(
			data,
			"{}{}_{}{}.json".format(
				global_config["data_dir"],
				global_config["session_id"],
				filename,
				global_config["suffix"]
			)
		)
	except:
		notify("JSON marshal failure for object: {} (will skip)".format(filename),
			level="warn")

def generate_coverage_profiles(include_genes=[], disabled=False):
	if disabled:
		return

	for gene, details in annotation.items():
		
		if len(include_genes) != 0:
			if gene not in include_genes:
				continue

		coverage_profile[gene] = get_coverage (
			start=details["start"],
			stop=details["stop"],
			reference=details["refname"],
			strand=details["strand"],
			mode="profile"
		).copy()

	# dump the coverage information
	save_file(coverage_profile, "igv_profile")


def get_coverage(start, stop, reference, strand, mode="mean"):

    bamfile = global_config["bamfile"]

    def watson_only(read):
        if read.is_reverse:
            return False
        if not read.is_reverse:
            return True

    def crick_only(read):
        if read.is_reverse:
            return True
        if not read.is_reverse:
            return False

    if strand == "+":
        cv = bamfile.count_coverage(reference,
                                    start=start,
                                    stop=stop,
                                    read_callback=watson_only)
    if strand == "-":
        cv = bamfile.count_coverage(reference,
                                    start=start,
                                    stop=stop,
                                    read_callback=crick_only)

    if mode ==  "mean":
        try:
            mean_cov = np.sum(np.sum(cv, axis=0))/len(cv)
        except ZeroDivisionError:
            mean_cov = 0
        return(mean_cov)

    if mode == "profile":
        cv = np.sum(cv, axis=0)
        return(
            {
                "coverage": cv.tolist(),
                "position": np.arange(start, stop, 1).tolist()
            }
        )



def plot_metagene_per_readlength(plot_title="", figsize_x=10,
	figsize_y=16, labsize=10):

        up       = global_config["endmap_upstream"]
        dn       = global_config["endmap_downstream"]
        nt_range = global_config["peak_scan_range"]
        data     = metagene_for_len

        savedata = {}

        MAX_PLOTS = len(data)
        xaxis     = str(list(data.keys())[0]) # x-axis identifier
        x_vector  = data[xaxis][:(up+dn)] # trim it down
        savedata["xaxis"] = x_vector

        fig, ax   = plt.subplots(nrows=MAX_PLOTS-1, ncols=1, sharex=True,
        	sharey=True, figsize=(figsize_x,figsize_y))

        plotcount = 1 # for plot
        axiscount = 0 # for axis
        for rpf_len, metagene_vector in data.items():

                savedata[rpf_len] = []
                savedata[rpf_len].append([metagene_vector.tolist()])

                if rpf_len == xaxis: continue
                plotcount += 1
                if plotcount == MAX_PLOTS+1: break

                palette = next(colors)["color"]

                ax[axiscount].plot(
                	x_vector, metagene_vector,
                	label="{} mer".format(rpf_len),color=palette
                )

                ax[axiscount].xaxis.grid (
                	which="major", color="k",
                	linestyle="-", linewidth=0.2
                )

                ax[axiscount].legend(loc="upper right")

                # In case we are given a range of values within which to
                # idenify and annotate the p-site offsets
                if len(nt_range) == 2:
                        peakpos, offset = get_offset (
                                xvec        = x_vector,
                                yvec        = metagene_vector,
                                left_bound  = nt_range[0],
                                right_bound = nt_range[1]
                        )

                        if global_config["metagene_callcount"] == 0:

                            ax[axiscount].axvline(
                                x=peakpos, color="black",
                                linestyle="--"
                             )
                            savedata[rpf_len].append([peakpos])

                            ax[axiscount].text(
                                peakpos+0.4, 0, offset,
                                rotation=0, color="black",
                                weight="medium"
                            )

                ax[axiscount].axvline(x=0, color="red", linestyle="--")
                axiscount += 1

        fig.suptitle(plot_title)
        fig.text(0.5, 0.04, 'Nt. from CDS start', ha='center',
         fontsize=labsize)
        fig.text(0.04, 0.5, 'Mean normalized count', va='center',
         rotation='vertical', fontsize=labsize)

        fname_line = "" # Ugly.
        fname_bar  = ""

        if global_config["metagene_callcount"] == 0:
            fname_line = "metagene_kmer_no_offset"
            fname_bar  = "metagene_kmer_framing_no_offset"
        else:
            fname_line = "metagene_kmer_offset"
            fname_bar  = "metagene_kmer_framing_offset"


        save_file(savedata, fname_line)

        plt.tight_layout()
        plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))

        plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
        	global_config["session_id"],
            fname_line,
            global_config["img_id"],
            global_config["offset_for_site"]))

        plt.close()


        # Phasing plots
        fig, ax = plt.subplots(nrows=MAX_PLOTS-1, ncols=1, sharex=True,
        	sharey=True, figsize=(figsize_x,figsize_y))

        plotcount = 1 # for plot
        axiscount = 0 # for axis
        for rpf_len, metagene_vector in data.items():

                if rpf_len == xaxis: continue
                plotcount += 1
                if plotcount == MAX_PLOTS+1: break

                bars = ax[axiscount].bar(x_vector[50:100], metagene_vector[50:100],
                    label="{} mer".format(rpf_len),color="#11a0ff")

                for item in bars[::3]:
                    item.set_color("#ff5842")

                ax[axiscount].xaxis.grid(which="major", color="k",
                    linestyle="-", linewidth=0.2)
                ax[axiscount].legend(loc="upper right")

                axiscount += 1

        fig.suptitle(plot_title)
        fig.text(0.5, 0.04, 'Nt. from CDS start', ha='center',
            fontsize=labsize)
        fig.text(0.04, 0.5, 'Mean normalized count', va='center',
         rotation='vertical', fontsize=labsize)

        plt.tight_layout()
        plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))

        plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
            global_config["session_id"],
            fname_bar,
            global_config["img_id"],
            global_config["offset_for_site"]))

        plt.close()

        global_config["metagene_callcount"] += 1

def rpkm_per_gene(gene_annotation):

	"""
	RPKM for whole length of gene
	"""
	strandinfo           = gene_annotation["strand"]
	read_lengths         = global_config["use_readlengths"]
	genomic_start        = gene_annotation["start"]
	genomic_stop         = gene_annotation["stop"]
	gene_length          = genomic_stop - genomic_start
	reference_id         = gene_annotation["refname"]
	bamfile              = global_config["bamfile"]
	mapped_reads_library = global_config["total_mapped_reads"]

	mapped_reads_gene    = 0

	for read in bamfile.fetch(reference=reference_id, start=start, stop=stop):

		if read.reference_length not in read_lengths:
			continue
		if ((strandinfo == "+") and (read.is_reverse)):
			continue
		if ((strandinfo == "-") and (not read.is_reverse)):
			continue

		mapped_reads_gene += 1

	scaling_factor = mapped_reads_library/1e6 #per million mapped reads
	rpm = mapped_reads_gene/scaling_factor # normalize counts for read depth, reads per million
	rpkm = rpm/((gene_length)/1e3) # rpkm = rpm/gene length in kB
	return(rpkm)


# INPUT:  Parsed annotation dictionary for a specific gene
# OUTPUT: The mean gene coverage value
def gene_coverage(gene_annotation):
	# Return the read coverage of a particular gene so that we can sort it later
	strand 	           = gene_annotation["strand"]
	genomic_start      = gene_annotation["start"]
	genomic_stop       = gene_annotation["stop"]
	reference_id       = gene_annotation["refname"]
	gene_length        = genomic_stop - genomic_start
	ignore_terminal_nt = global_config["readcov_ignore_nts"]
	from_which_end     = global_config["readcov_from_terminal"]
	percentage_length  = global_config["readcov_percentage_length"]
	length_to_consider = (percentage_length/100) * gene_length

	# Is the gene long enough for coverage calculation?
	if ((gene_length >= (ignore_terminal_nt+2)) == False):
        	return(0)

	# Determine the genomic coordinate range within which to calculate the coverage.
	if strand == "+":
		if from_which_end == "five_prime":
			start = genomic_start + ignore_terminal_nt
			stop  = genomic_start + percentage_length
		if from_which_end == "three_prime":
			start = genomic_stop - percentage_length
			stop  = genomic_stop

	if strand == "-":
		if from_which_end == "five_prime":
			start = genomic_stop - percentage_length
			stop  = genomic_stop - ignore_terminal_nt
		if from_which_end == "three_prime":
			start = genomic_start
			stop  = genomic_start + percentage_length

	# Get the mean read coverage for a specific length of the genome
	mc = get_coverage(start=start, stop=stop, strand=strand,
						reference=reference_id, mode="mean")
	return (mc)


def calculate_gene_coverages(mode="mean_coverage"):

	tmp = {} # temporarily store coverages before sorting

	if mode == "mean_coverage":
		for gene_name, details in annotation.items():
			tmp[gene_name] = gene_coverage(details)

	if mode == "rpkm":
		for gene_name, details in annotation.items():
			tmp[gene_name] = rpkm_for_gene(details)

	tmp_tuple = sorted(tmp.items(),
		key=lambda item: item[1],
		reverse=True
	)

	for index in range(len(tmp_tuple)):
		gene_coverages[tmp_tuple[index][0]] = tmp_tuple[index][1]

	save_file(gene_coverages, "gene_coverages")


def plot_combined_metagene(peak=False, peak_range=[]):
    # Plots the metagene profile taking into account all lengths of
    # the RPF fragments.
    up = global_config["endmap_upstream"]
    dn = global_config["endmap_downstream"]
    xaxis = metagene_for_len["xaxis"][:(up+dn)]

    combined_vectors = []

    for length, vector in metagene_for_len.items():
        if length == "xaxis":
            continue
        combined_vectors.append(vector.tolist()[:(up+dn)])

    metagene_vector = np.sum(combined_vectors, axis=0)

    plt.plot(xaxis, metagene_vector, color="black")

    if peak and len(peak_range) != 0:

        peakpos,delta = get_offset(xvec=xaxis,
        	yvec=metagene_vector, right_bound=peak_range[1],
        	left_bound=peak_range[0])
        plt.axvline(x=delta, color="r")
        plt.text(peakpos+0.4, 0, delta,rotation=0, color="black",
        weight="bold") # 0.4 to keep the text away from the line

    start_data = { # include the plot configurations as well?
        'peakpos': peakpos,
        'delta': delta,
        'xaxis': xaxis.tolist(),
        "start_metagene": metagene_vector.tolist()
    }

    save_file(start_data, "start_position_peak")

    plt.xlabel("Nt. from CDS start")
    plt.ylabel("Total mapped reads")
    # plt.grid()

    # plt.tight_layout()
    plt.gcf().set_size_inches(18, 4)
    plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
    plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
        global_config["session_id"],
        "start_peak",
        global_config["img_id"],
        global_config["offset_for_site"]))

    plt.close()


    bars = plt.bar(xaxis[50:170], metagene_vector[50:170], color="#11a0ff")

    for item in bars[::3]:
        item.set_color("#ff5842")
    plt.xlabel("Nt. from CDS start")
    plt.ylabel("A.U.")

    # plt.tight_layout()
    plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
    plt.gcf().set_size_inches(18, 4)
    plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
         global_config["session_id"],
        "start_peak_framing",
        global_config["img_id"],
        global_config["offset_for_site"]))

    plt.close()



def get_max_coverage_genes(percentage=20):
	# Identify top N% highest coverage genes
	# We can't subset a dict,so we loop through the top N genes
	# in the sorted dict

	if len(gene_coverages) == 0:
		calculate_gene_coverages()

	MAX_COUNT = np.ceil((percentage/100)*len(gene_coverages))
	i = 0
	for gene_name, gene_coverage in gene_coverages.items():
		i += 1
		high_coverage_genes.append(gene_name)
		if i == MAX_COUNT: break


	save_file(high_coverage_genes, "high_coverage_genes")

	return(high_coverage_genes)



def skip_read(read,      # the read sequence
              skip_dict, # codons, if found in E/P/A site would be skipped
              offset,    # The p site offset for the read
              strand):   # Strand the read maps to
    # Should a read be skipped or not?
    # Skip only if read has codon in any of it's sites
    skip    = True
    E, P, A = "", "", "" # the E/P/A site codons of a read

    if strand == "+": # read maps to watson
        E = read[offset-3: offset]
        P = read[offset: offset+3]
        A = read[offset+3: offset+6]

    if strand == "-": # read maps to watson
        # reverse the codon
        offset = len(read)-offset
        A = read[offset-6: offset-3][::-1]
        P = read[offset-3: offset][::-1]
        E = read[offset: offset+3][::-1]

    if ((E in skip_dict["E"]) or
        (P in skip_dict["P"]) or
        (A in skip_dict["A"])):
            skip = False

    return(skip)


# INPUT:  1. Gene annotation dictionary for a specific gene
#         2. The PRF terminal we want to map (5' or 3')
# OUTPUT: 1. Defaultdict, mapping the length of an RPF to the
#            list of genomic coordinates of it's 5'/3' terminals
#	     eg: endmap[28] = [1000,1002,1050,1004, ...]
#         2. Small config dictionary containing the start and
#	         stop coordinates after padding (upstream+buffer)
#            along with the strand information
def isolate_readprofile_by_size(gene_annotation,  apply_offset = False,
							ignore_reads_by_nt = {5: [], 3: []},
							codons_to_skip = {'E':[], 'P':[], 'A': []}):
	'''

	Codons_to_skip is a list of codons in either the p,a or e site of a read to skip.
	Contains nested lists, for eg: ['AGG', 'E'] will skip all reads with an 'AGG'
	codon in the 'E' site. WARNING: This only works correctly when the offsets
	supplied are p-site offsets. Supplying an A-site offset and then passing a
	['AGG', 'A'] will only cause a further shift of 3nt and the wrong codon will
	be selected.

	TL;DR:  If you want to skip reads with ABC codon in E site, you have two options:
	        1. Use the p-site offset calculation, and pass in ['ABC', 'E']
	        2. Use the e-site offset calculatio, and pass in ['ABC', 'P']

	Ignore reads per nucleotide is a dictionary of the following structure:
		{
			5: ["A", "G" ...], # reads with these nucleotides at the 5' end will be ignored
			3: ["A", "G" ...]  # reads with these bases at the 3' end will be ignored
		}

	'''

	# print("entered function, stratified mappings")
	# Defaultdict to store a list of coordinates to which a
	# kmer of a particular RPF length is mapping to.
	stratified_mappings = defaultdict(list)
	# The start and stop values after the strand specific
	# buffers and upstream/downstream regions have been specified.
	# Also contains the strand information, for use when we want
	# to convert the mapping to a vector as in the previous plot.
	config              = {}
	genomic_start       = gene_annotation["start"]
	genomic_stop        = gene_annotation["stop"]
	gene_length         = abs(genomic_start - genomic_stop)
	strandinfo          = gene_annotation["strand"]
	reference_id        = gene_annotation["refname"]
	upto_upstream_nt    = global_config["endmap_upstream"]
	upto_downstrm_nt    = global_config["endmap_downstream"]
	bamfile             = global_config["bamfile"]
	_buffer_	        = global_config["mapping_buffer"]
	terminal_map        = global_config["read_offset_terminal"]
	read_lengths        = global_config["use_readlengths"]
	offset_site         = global_config["offset_for_site"]



	if offset_site == "P":
		offsets = offsets_p
	if offset_site == "A":
		offsets = offsets_a
	if offset_site == "E":
		offsets = offsets_e

	if apply_offset and (len(offsets) == 0):
		notify("Offsets won't be applied", type="warn", onetime=True)
		apply_offset = False

	if strandinfo == "+":
		# watson strand genes
		start = genomic_start - (upto_upstream_nt + _buffer_)
		stop  = genomic_start + (upto_downstrm_nt + _buffer_) + gene_length

	if strandinfo == "-":
		# crick strand genes
		start = genomic_stop  - (upto_downstrm_nt + _buffer_) - gene_length
		stop  = genomic_stop  + (upto_upstream_nt + _buffer_)

	# For return
	config["start"]  = start
	config["stop"]   = stop

	# print("converted coordinates")
	# Ignore genes too close to the actual genomic start coordinate
	# Including buffers and all.
	if (start <= 0):
		# we return an empty dictionary in that case
		# Will ignore genes very close to the chromosomal start
		return(stratified_mappings, config)

	# print("entering bam file")
	# Get reads mapping to a specific reference genomic position
	i, j = 0, 0

	for read in bamfile.fetch(reference=reference_id, start=start, stop=stop):

		# We need to know which strand a read is mapping to in order to figure out it's
		# 5' or 3' ends
		# watson_read = True

		if read.reference_length not in read_lengths:
			continue

		if ((strandinfo == "+") and (read.is_reverse)):
			continue

		if ((strandinfo == "-") and (not read.is_reverse)):
			continue

		# if read.is_reverse:
		# 	watson_read = False

		# Does not take into account the soft/hard clipping of the the read sequences
		# Probably should compare the cigars?
		# read_sequence = read.query_sequence

		# if watson_read:
		# 	base5  = read_sequence[0]
		# 	base3  = read_sequence[len(read_sequence)-1]
		# else:
		# 	base5  = read_sequence[len(read_sequence)-1]
		# 	base3  = read_sequence[0]

		# # Any reads with set nucleotides from the 5' or
		# # 3' ends will be skipped.
		# if (len(ignore_read_by_nt[5]) != 0):
		# 	if base5 in ignore_read_by_nt[5]:
		# 		# skip the read
		# 		continue
		# if (len(ignore_read_by_nt[3]) != 0):
		# 	if base3 in ignore_read_by_nt[3]:
		# 		# skip the read
		# 		continue

		if ((len(codons_to_skip["E"]) != 0) or
			(len(codons_to_skip["P"]) != 0) or
			(len(codons_to_skip["A"]) != 0)):

			seq   = read.query_sequence
			delta = abs(value_from_dict(offsets_p, # we only consider the p-site offset
						read.reference_length))

			if read.is_reverse:
				if skip_read(seq, codons_to_skip, delta, "-"):
					i += 1
					continue

			if not read.is_reverse:
				if skip_read(seq, codons_to_skip, delta, "+"):
					j += 1
					continue


		terminal_coordinate = 0

        # if gene is in positive strand, the five_prime terminal is different
        # for watson and crick reads. Include that information

        # 5' coordinate of watson strand read
		if terminal_map == "five_prime" and not read.is_reverse:
			terminal_coordinate = (read.reference_start+1)
		# 5' coordinate of crick strand read
		if terminal_map == "five_prime" and read.is_reverse:
			terminal_coordinate = read.reference_end
        # 3' coordinate of watson strand read
		if terminal_map == "three_prime" and not read.is_reverse:
			terminal_coordinate = read.reference_end
		# 5' coordinate of crick strand read
		if terminal_map == "three_prime" and read.is_reverse:
			terminal_coordinate = (read.reference_start+1)
		# why are we using np.array() here again?
		if apply_offset == False:
			stratified_mappings[read.reference_length].append (
				np.array(terminal_coordinate)
			)

		if apply_offset == True and terminal_map == "five_prime":

			if not read.is_reverse:
				stratified_mappings[read.reference_length].append (
					np.array(terminal_coordinate) +
					abs(value_from_dict(offsets,
						read.reference_length))
				)
			if read.is_reverse:
				stratified_mappings[read.reference_length].append (
					np.array(terminal_coordinate) -
					abs(value_from_dict(offsets,
						read.reference_length))
				)

		if apply_offset == True and terminal_map == "three_prime":

			if not read.is_reverse:
				stratified_mappings[read.reference_length].append (
					np.array(terminal_coordinate) -
					abs(value_from_dict(offsets,
						read.reference_length))
				)
			if read.is_reverse:
				stratified_mappings[read.reference_length].append (
					np.array(terminal_coordinate) +
					abs(value_from_dict(offsets,
						read.reference_length))
				)

	return (stratified_mappings, config)


# INPUT:  1.A defaultdict of RPF lengths, mapping to a list
#           of genomic coordinates of it's terminals
#	  2.The annotation data of the gene in question
#         3.The buffered start and stop coordinates over which
#           to calculate the vectors
# OUTPUT: A dictionary mappning, for one specific gene, the RPF
#         length to a vector of adjusted genomic positions containing
#         the counts of the number of RPF terminals mapping to
#         those positions
def convert_readmap_to_vectors(mapdict, config, gene_annotation):
	# Convert the defaultdict of stratified kmer mappings
	# into a dictionary of padded vectors

	# buffered genomic coordinates from the start codon,
	# considering the upstream and downstream codons from
	# start
	start            = config["start"]
	stop             = config["stop"]
	# Gene annotation details
	strand 	         = gene_annotation["strand"]
	genomic_start    = gene_annotation["start"]
	genomic_stop     = gene_annotation["stop"]
	gene_length      = abs(genomic_start - genomic_stop)
	# How many nt up and downstream of the start codon are
	# we supposed to consider
	upto_upstream_nt = global_config["endmap_upstream"]
	upto_downstrm_nt = global_config["endmap_downstream"]


	# Store the genonic endmap coordinates converted so
	# that the start codon gets index 0.
	tmp_maps         = {}
	# Convert the defaultdict list to numpy arrays for
	# easier hadling
	for length, mappings in mapdict.items():
		tmp_maps[length] = np.array(mappings)

	# Reset genomic mappings so that start codon gets index 0
	if strand == "+":
		for length, mappings in tmp_maps.items():
			tmp_maps[length] = mappings - genomic_start
	if strand == "-":
		for length, mappings in tmp_maps.items():
			tmp_maps[length] = genomic_stop - mappings

	# A numpy array holding the indexes around the genomic region
	indexes  = np.array([i for i in range(-upto_upstream_nt,
			upto_downstrm_nt+gene_length+1,1)
	])
	# Just set the index for the metagene hash for plotting,
	# we will populate it later down the line
	metagene_for_len["xaxis"] = indexes
	# Zeroed vector, holding counts for the number of RPF terminals
	# mapping to a particular position in index.
	template = np.zeros(len(indexes))

	# Convert endmaps to vectors
	for length, mappings in tmp_maps.items():

		nullvec = template.copy()

		for maps in mappings:
			# Ignore all RPFs whose terminals are mapping between
			# the nt_upstream/downstream and buffer coordinates
			if ((maps >= (-upto_upstream_nt)) and
				(maps <= upto_downstrm_nt + gene_length)):

				index = getindex(indexes, maps)
				if index == None:
					# Element doesn't exist
					continue
				# print(length, maps, index)
				nullvec[index] = nullvec[index] + 1

		tmp_maps[length] = nullvec.copy()

		nullvec = np.delete(nullvec,
			[i for i in range(len(nullvec))]
		)

	return(tmp_maps)


def identify_overlaps():

	# This is junk. We can do better.

	overlap_fiveprime    = [] # genes overlapping at their 5' ends
	overlap_threeprime   = [] # genes overlapping at their 3' ends

        # Convert the dictionary into a list of lists
	annotation_list      = []
	for gene_name, details in annotation.items():
		temp = [gene_name, details]
		annotation_list.append(temp)


	for i in range(len(annotation_list)):

		a = i
		b = a + 1

		if b == len(annotation_list):
			break

		gene_upstream   = annotation_list[a]
		gene_downstream = annotation_list[b] # <-- does this overlap?
                # A gene is overlapping from the five prime end
                # if it's start genomic coordinate is less than the
                # stop genomic coordinate of the downstream gene.
                # list of gene_name and dict of annotation
		stop_gene_upstream     = gene_upstream[1]["stop"]
		start_gene_downstream  = gene_downstream[1]["start"]

                # We need to consider the strandedness when we assign overlaps
                # they need to be same, otherwise we ignore
		upstream_gene_strand   = gene_upstream[1]["strand"]
		downstream_gene_strand = gene_downstream[1]["strand"]

		if ((start_gene_downstream < stop_gene_upstream) and
                   (upstream_gene_strand == downstream_gene_strand)):
                    # The gene is overlapping
                    # We record the name of the gene
			overlap_fiveprime.append(gene_downstream[0])
			overlap_threeprime.append(gene_upstream[0])

		a += 1
		b += 1

	if ((len(overlap_fiveprime) != 0) or
		(len(overlap_threeprime) != 0)):
		notify("There are overlapping genes", level="warn", onetime=True)

	# genes overlapping at their 5' end
	overlap_dict["five_prime"]  = overlap_fiveprime
	# genes overlapping at their 3' end
	overlap_dict["three_prime"] = overlap_threeprime



def map_gene_to_endmaps(apply_offset=False,
				only_reads_with_codons = {'E':[], 'P':[], 'A':[]}):
	#
	ignore_overlap = global_config["will_ignore_overlaps"]

	try:
		overlap_exception = global_config["overlap_exception"]
	except:
		overlap_exception = []

	# Diagnostic purposes only.
	track_overlaps  = []
	track_coverages = []
	overall_skipped = []
	genes_processed = []
	genes_used      = []

	notify("{} genes in high coverage list, out of {}".format(len(high_coverage_genes),
		len(annotation.keys())), level="notf")

	for gene_name, details in annotation.items():

		genes_processed.append(gene_name)

		if gene_name not in high_coverage_genes:
			track_coverages.append(gene_name)
			overall_skipped.append(gene_name)
			continue

		if ((ignore_overlap) and (gene_name not in overlap_exception)):
			# The "start" of a gene on the crick strand
			# is actually it's genomic stop coordinate
			if ((details["strand"] == "-") and
			(gene_name in overlap_dict["three_prime"])):
				track_overlaps.append(gene_name)
				overall_skipped.append(gene_name)
				continue

			if ((details["strand"] == "+") and
			(gene_name in overlap_dict["five_prime"])):
				track_overlaps.append(gene_name)
				overall_skipped.append(gene_name)
				continue

		genes_used.append(gene_name)

		maps, cnfs = isolate_readprofile_by_size(details, apply_offset,
									 codons_to_skip = only_reads_with_codons)
		endmaps    = convert_readmap_to_vectors(maps,cnfs,details)
		per_gene_endmaps[gene_name] = endmaps

	notify("Intersection: {}".format(len(set(high_coverage_genes) & set(genes_used))))
	notify("{} genes skipped due to overlaps".format(len(track_overlaps)), level="warn")
	notify("{} genes skippped due to low coverage".format(len(track_coverages)), level="warn")
	notify("{} genes processed in total".format(len(genes_processed)), level="notf")
	notify("{} genes skipped out of {} in annotation in total".format(len(overall_skipped),
		len(annotation.keys())), level="warn")
		
	# print out the gene names which were in the high coverage list but skipped
	save_file(unique(track_overlaps),  "SKIPPED_GENE_OVERLAPS_ENDMAPS.json")
	save_file(unique(track_coverages), "SKIPPED_GENE_COVERAGE_ENDMAPS.json")
	save_file(unique(genes_used),      "GENES_USED_ENDMAPS.json")


def normalize_counts(vector):
	# Normalize a numpy array by it's sum
	if sum(vector) == 0:
		return vector
	else:
		return(vector/sum(vector))


def frame_for_coordinate(genomic_start, genomic_stop, query_index, strand):

	if (query_index > genomic_stop) or (query_index < genomic_start):
		return 0

	if strand == "+":
		delta = query_index - genomic_start
		return((delta % 3)+1)
	if strand == "-":
		delta = genomic_stop - query_index
		return((delta % 3)+1)


def reading_frame(ignore_read_by_nt = {5: [], 3: []}):

	reads_in_frame   = {}
	kmer_include     = global_config["use_readlengths"]

	offset_direction = global_config["read_offset_terminal"]

	offset_site      = global_config["offset_for_site"]

	if offset_site == "P":
		offset_dict = offsets_p
	if offset_site == "A":
		offset_dict = offsets_a
	if offset_site == "E":
		offset_dict = offsets_e

	reads_in_frame["total"] = 0
	reads_in_frame[1]       = 0
	reads_in_frame[2]       = 0
	reads_in_frame[3]       = 0

	total, processed, skipped, skipped_other = 0, 0, 0, 0
	for gene_name, details in annotation.items():

		genomic_start = details["start"] # 1 indexed
		genomic_stop  = details["stop"]
		strand        = details["strand"]
		refname       = details["refname"]
		bamfile       = global_config["bamfile"]

		# Ignore genes which are too short or have low coverage
		if ((abs(genomic_start-genomic_stop) < 100) or
			(gene_name not in high_coverage_genes)):
			continue

		# 27 and 12 from DOI: 10.7554/eLife.42591.001
		if strand == "+":
			start = genomic_start + 27
			stop  = genomic_stop  - 12

		if strand == "-":
			start = genomic_start + 12
			stop  = genomic_stop  - 27

		for read in bamfile.fetch(reference=refname, start=start, stop=stop):
			total += 1

			# we need the strand of the read to figure out it's 5' or 3' ends
			watson_read = True

			if ((strand == "+") and (read.is_reverse)):
				skipped_other += 1
				continue
			if ((strand == "-") and (not read.is_reverse)):
				skipped_other += 1
				continue

			# Skip reads with a specific base at their 5' or 3' terminals
			read_sequence = read.query_sequence

			if read.is_reverse:
				watson_read = False

			if watson_read:
				base5  = read_sequence[0]
				base3  = read_sequence[len(read_sequence)-1]
			else:
				base5  = read_sequence[len(read_sequence)-1]
				base3  = read_sequence[0]

			# Any reads with set nucleotides from the 5' or
			# 3' ends will be skipped.
			if (len(ignore_read_by_nt[5]) != 0):
				if base5 in ignore_read_by_nt[5]:
					skipped += 1
					continue

			if (len(ignore_read_by_nt[3]) != 0):
				if base3 in ignore_read_by_nt[3]:
					skipped += 1
					continue

			# Only calculate frame info for a particular length
			if ((len(kmer_include) != 0) and
			(read.reference_length not in kmer_include)):
				skipped_other += 1
				continue

			read_map = 0
			offset   = 0

			if len(offset_dict) != 0:
				try:
					offset = offset_dict[read.reference_length]
				except:
					notify("Some offsets forced to 0, {} length RPF has no offset".
						format(read.reference_length), level="warn")
					# We ignore a read which has no corresponding offset
					skipped_other += 1
					continue

			if (strand == "+") and (offset_direction == "five_prime"):
					read_map   = (read.reference_start+1) + abs(offset)

			if (strand == "+") and (offset_direction == "three_prime"):
					read_map   = read.reference_end  - abs(offset)

			if (strand == "-") and (offset_direction == "five_prime"):
					read_map   = read.reference_end  - abs(offset)

			if (strand == "-") and (offset_direction == "three_prime"):
					read_map   = (read.reference_start+1) + abs(offset)

			read_frame = frame_for_coordinate(genomic_start, genomic_stop, read_map, strand)
			processed += 1

			try:
				reads_in_frame[read_frame] += 1
				reads_in_frame["total"]    += 1
			except:
				skipped_other += 1
				pass

	notify("{} reads processed out of {} reads total ({}%)".format(processed, total,
		round((processed/total)*100,2)))

	notify("{} ({}%)reads skipped due to nt mismatch, {} ({}%) skipped due to other issues".format(
		skipped,
		round((skipped/total)*100,2),skipped_other,round((skipped_other/total)*100,2)))

	save_file(reads_in_frame, "global_reading_frame")

	return(reads_in_frame)


def metagene_over_codon(codon_list, mode="combined"):
	# Driver code

	# Plot a combined metagene for a list of codons
	# Inp: ["ABC", "DEF", "GHI", "JKL"]
	# Out: Single plot combining all codons
	if mode == "combined":
		plot_metagene_over_codon(codon_list)

	# Plot separate metagenes for a list of codons
	# Inp: ["ABC", "DEF", "GHI", "JKL"]
	# Out: Four separate plots for the codons
	if mode == "separate":
		for codon in codon_list:
			plot_metagene_over_codon([codon])

	# Plot individual combined metagenes for a nested
	# list of codon lists
	# Inp: [["ABC", "DEF"], ["GHI", "JKL"]]
	# Out: Two separate plots for ABC&DEF combined
	#      and GHI&JKL combined.
	if mode == "composite":
		for nested_list in codon_list:
			plot_metagene_over_codon(nested_list)


def plot_metagene_over_codon(codon_list):

	codon_positions = {}
	cmetagene       = [] # metagene over a specific codon

	if len(transcripts_dict) == 0:
		retrieve_coding_sequence(
			global_config["coding_sequence_path"],
			global_config["coding_sequence_format"]
		)

	i = 0
	x = ['-2','E','P','A','2']

	for gene_name, details in annotation.items():

		if gene_name not in high_coverage_genes:
			continue

		try:
			pps = pps_vector_per_gene[gene_name]["pps"]
		except:
			notify("{} does not have a pps vector".format(gene_name),
				level="warn")
			continue

		try:
			seq = pps_vector_per_gene[gene_name]["seq"]
		except:
			notify("{} does not have a transcript sequence".format(gene_name),
				level="warn")
			continue

		if i == 401:
			break

		i += 1

		tmp = {}
		# tmp["seq"] = seq
		tmp["pps"] = pps
		tmp["pos"] = position_of_codon("".join(ungroup(seq)), codon_list,
			mode="codon")
		codon_positions[gene_name] = tmp


		for gene, values in codon_positions.items():

			for position in values["pos"]:
				val = []
				try:
					val = chunk_from_id(values["pps"], position, int((len(x)-1)/2))
				except:
					continue
				cmetagene.append(val)

	data = {'position': x, 'positional_pause': cmetagene}
	save_file(data, "data_metagene_over_codon_{}".format("+".join(codon_list)))

	cmetagene = np.apply_along_axis(np.mean, 0, np.array(cmetagene))
	# cmetagene = np.sum(cmetagene, axis=0)

	plt.plot(x, cmetagene, color="black")
	plt.xlabel("Position (nt)")
	plt.ylabel("Mean positional pause score")
	plt.suptitle("Pauses over {} ({}site, {})".format("+".join(codon_list),
		global_config["offset_for_site"],
		global_config["prefix"]
	))

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"metagene_over_codon_{}".format("+".join(codon_list)),
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()



def plot_reading_frame(framedict):

	x = ["Frame 1", "Frame 2", "Frame 3"]
	y = np.array([framedict[1], framedict[2], framedict[3]]) / framedict["total"]

	save_file({'frame': x, 'fraction_reads': y.tolist()}, "data_global_reading_frame")

	plt.bar(x, y)[0].set_color('#d53f04')
	plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.ylabel("Fraction")
	plt.ylim(0,1)

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"reading_frame",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()



def plot_stop_peak(leftpos=10,rightpos=0):

	tmp = []
	for vector in per_gene_densities.values():
		if len(vector) < 150:
			continue

		tmp.append(vector[len(vector)-150:])

	tmp = np.array(tmp)
	metagene = np.sum(tmp, axis=0)
	x = np.array([i for i in range(0, 150,1)][::-1])
	plt.plot(x, metagene, color="black")

	peakpos, delta = get_offset(xvec=x,
		yvec=metagene, right_bound=rightpos,
		left_bound=leftpos)

	plt.axvline(x=peakpos, color="r")
	plt.text(peakpos+0.5, 0, delta,rotation=0, color="blue",
	weight="bold") # 0.4 to keep the text away from the line
	plt.xlabel("Nt. from CDS stop")
	plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.ylabel("A.U.")

	plt.xlim(150,0)

	# plt.tight_layout()
	plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.gcf().set_size_inches(18, 4)

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"stop_peak",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()


	data = {
		"peakpos": peakpos,
		"delta": delta,
		"positions": x.tolist(),
		"stop_pos_metagene": metagene.tolist()
	}
	save_file(data, "data_stop_metagene")


def calculate_pause_score():

	stopcodon  = set(["TAA", "TAG", "TGA"])

	cds_path   = global_config["coding_sequence_path"]
	cds_format = global_config["coding_sequence_format"]

	retrieve_coding_sequence(cds_path, cds_format)
	skipped_transcripts = []
	skipped_density     = []
	skipped_mismatch    = []
	skipped_truncated   = []
	skippped_coverage   = []

	for gene_name, details in annotation.items():

		if gene_name not in high_coverage_genes:
			skippped_coverage.append(gene_name)
			continue

		try:
			seq = transcripts_dict[gene_name]
		except:
			skipped_transcripts.append(gene_name)
			continue

		try:
			rds = per_gene_densities[gene_name]
		except:
			skipped_density.append(gene_name)
			continue

		if len(seq) != len(rds):
			skipped_mismatch.append(gene_name)
			continue

		if (len(seq) %3 != 0) or (len(rds) %3 != 0):
			skipped_truncated.append(gene_name)
			continue

		seq = groupN(seq, 3) # split to codons
		rds = groupN(rds, 3) # density per CODON
		pps = [] # positional pause scores per codonn

		# Ignore gene sequences with premature stop codons
		if (bool(stopcodon & set(seq[:len(seq)-1]))):
			notify("{} has premature stop codon. Skipping.".format(gene_name),
			 level="warn")
			continue


		seq = seq[1:len(seq)-1] # remove first and last codon
		rds = rds[1:len(rds)-1]

		# print((len(seq)*3), (len(rds)*3), (len(seq)*3) % 3, (len(rds)*3) % 3)

		codon_vs_density = zip(seq, rds)

		# calculate pause score
		sum_codon_densities = np.sum(rds)

		if sum_codon_densities == 0:
			# This gene has no reads mapped to it. Skipping.
			continue

		for codon, density in codon_vs_density:

			positional_pause_score = (sum(density)) / (sum_codon_densities/len(rds))
			pps.append(positional_pause_score)

			codon_pauses[codon].append(positional_pause_score)

		# append to global dict
		pps_vector_per_gene[gene_name] = {
			"seq": seq,
			"pps": pps
		}

	# Average out the per codon pause scores
	# to get the final pause score
	for codon, pauselist in codon_pauses.items():
		if len(pauselist) == 0:
			codon_pauses[codon] = 0
		else:
			# Igore all zero values
			pauselist = np.array(pauselist)
			pauselist = pauselist[pauselist != 0]
			codon_pauses[codon] = np.median(pauselist)


	save_file(pps_vector_per_gene, "positional_pauses_universal")
	save_file(codon_pauses, "per_codon_pause_scores_universal")

	notify("{} genes not in transcripts".format(len(skipped_transcripts)), level="warn")
	notify("{} genes not in densities".format(len(skipped_density)), level="warn")
	notify("{} genes sequences-density mismatch".format(len(skipped_mismatch)), level="warn")
	notify("{} genes not modulo 3 in length".format(len(skipped_truncated)))
	notify("{} genes skipped due to coverage".format(len(skippped_coverage)))

	# Write out genes skipped in pause score
	save_file(unique(skipped_transcripts), "SKIPPED_GENE_NO_SEQUENCE_PPS.json")
	save_file(unique(skipped_density),     "SKIPPED_GENE_NO_DENSITY_PPS.json")
	save_file(unique(skipped_mismatch),    "SKIPPED_GENE_SEQ_RDS_LEN_MISMATCH_PPS.json")
	save_file(unique(skipped_truncated),   "SKIPPED_GENE_NOT_MULTIPLE3_PPS.json")
	save_file(unique(skippped_coverage),   "SKIPPED_GENE_LOW_COVERAGE_PPS.json")




def plot_aa_pause(aa):
    amino_acid = aa["AA"]
    aapause   = aa["Pause"]
    aatable = pd.DataFrame({
        "AA": amino_acid,
        "Type": ['Hydrophobic',
     'PolarNeutral',
     'PolarNegative',
     'PolarNegative',
     'Hydrophobic',
     'Hydrophobic',
     'PolarPositive',
     'Hydrophobic',
     'PolarPositive',
     'Hydrophobic',
     'Hydrophobic',
     'PolarNeutral',
     'Hydrophobic',
     'PolarNeutral',
     'PolarPositive',
     'PolarNeutral',
     'PolarNeutral',
     'Hydrophobic',
     'Hydrophobic',
     'Hydrophobic'],
        "pause": aapause
    })

    aatable = aatable.sort_values(by=['Type'])
    #aatable.to_csv("file.csv", index=False)
    fig, ax = plt.subplots()

    def change_width(ax, new_value) :
    	# 
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - new_value
            patch.set_width(new_value)
            patch.set_x(patch.get_x() + diff * 0.7) #recenter

    sns.set_style("white")
    palette = ["#3498db", "#95a5a6", "#e74c3c", "#2ecc71"]
    sns.set_palette(palette)
    avg = np.mean(list(aatable["pause"]))
    std = np.std(list(aatable["pause"]))
    plt.axhline(avg, color="black")
    plt.axhline(avg+std, color="grey")
    plt.axhline(avg-std, color="grey")
    sns.barplot(x="AA", y="pause", data=aatable, hue="Type")
    plt.xlabel("Amino acid")
    plt.ylabel("Mean pause score")
    change_width(ax, .45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
    plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
    global_config["session_id"],
    "pause_score_aa_annotated",
     global_config["img_id"],
     global_config["offset_for_site"]))

    plt.close()


def plot_pauses():

	# Plot positional pauses metagene first
	pps_metagene = []
	for vector in pps_vector_per_gene.values():
		if len(vector["pps"]) < 150:
			continue
		pps_metagene.append(vector["pps"][:150])
	xvec = [i for i in range(1,151,1)]
	pps_metagene = np.mean(pps_metagene, axis=0)
	plt.plot(xvec, pps_metagene,color="k")
	plt.xlabel("#codons from CDS start")
	plt.ylabel("PPS (Avg over all genes)")
	# plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	# plt.tight_layout()
	plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.gcf().set_size_inches(18, 4)

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"pps_over_genes",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()


	pps_metagene_data = {
		"XAxis": xvec,
		"PPScore": pps_metagene.tolist()
	}

	# PLot the pause score heat map
	codon = []
	pause = []

	for key, value in codon_pauses.items():

		if key in ["TAA", "TAG", "TGA"]:
			continue

		codon.append(key)
		pause.append(value)

	# fig, ax = plt.subplots()
	# seaborn.heatmap([vec1])

	# plt.bar(codon, pause, color="dodgerblue")
	# plt.axhline(y=np.mean(pause)+np.std(pause))
	# plt.axhline(y=np.mean(pause))
	# plt.axhline(y=np.mean(pause)-np.std(pause))
	# plt.xticks(rotation=90)
	# plt.show()

	# sns.heatmap([np.asarray(pause)], xticklabels=codon, cmap='RdBu_r', vmin=-0.5, vmax=2.5)
	sns.heatmap([np.asarray(pause)], xticklabels=codon, cmap='RdBu_r', vmax=2.5, center=1, vmin=-0.5) 

	plt.tight_layout()
	plt.gcf().set_size_inches(18, 5)
	plt.title("{} ({} site, {})".format(
		global_config["prefix"],global_config["offset_for_site"], global_config["read_offset_terminal"]
		)
	)

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"pause_score_codon",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()


	codon_pause_data = {
		"Codon": codon,
		"Pause": pause
	}

	# Show pause scores per amino acid
	# append pause score per codon to the appropriate AA, then avg

	pause_per_aa = {}

	for codon, avgpause in codon_pauses.items():

		AA  = codon_keys[codon]

		if AA in pause_per_aa.keys():
			pause_per_aa[AA].append(avgpause)
		else:
			pause_per_aa[AA] = [avgpause]

	del pause_per_aa["*"]

	aa = [] # amino acid
	ps = [] # avg pause per AA

	for key, value in pause_per_aa.items():
		aa.append(key)
		ps.append(np.mean(value))

	# Sort the lists for easier plotting. This is an extended sort.
	aa, ps = zip(*sorted(zip(aa, ps)))

	aa_pause_data = {
		"AA": aa,
		"Pause": ps
	}

	plot_aa_pause(aa_pause_data)

	# sns.heatmap([np.asarray(ps)], xticklabels=aa, cmap="RdBu_r", vmin=-0.5, vmax=2.5)
	sns.heatmap([np.asarray(ps)], xticklabels=aa, cmap="RdBu_r", center=1, vmax=2.5, vmin=-0.5)

	plt.tight_layout()
	
	plt.title("{} ({} site, {})".format(
		global_config["prefix"],global_config["offset_for_site"], global_config["read_offset_terminal"]
		)
	)

	plt.gcf().set_size_inches(18, 5)

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"pause_score_aa",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()



	save_file(pps_metagene_data, "positional_pauses")
	save_file(codon_pause_data, "per_codon_pause_scores")
	save_file(aa_pause_data, "per_aa_pause_scores")


def read_length_histogram():

	readlens = []
	bamfile  = global_config["bamfile"]

	for gene, details in annotation.items():

		start  = details["start"]
		stop   = details["stop"]
		strand = details["strand"]
		refer  = details["refname"]

		for read in bamfile.fetch(start=start, stop=stop, reference=refer):

			if strand == "+" and read.is_reverse:
				continue
			if strand == "-" and not read.is_reverse:
				continue
			# if read.reference_length > 50:
			# 	notify("Encountered read lengths > 50nt", level="warn", onetime=True)
			# 	continue

			readlens.append(read.reference_length)

	#bins = [i for i in range(20,41,1)]

	plt.hist(readlens, edgecolor="black",color="#808080")
	# plt.xticks(bins)
	plt.title("Read length distribution")
	plt.xlabel("Read lengths")
	plt.ylabel("Frequency")
	plt.title("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"read_length_histogram",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()

	save_file(readlens, "readlengths")


def read_terminal_stats():
	# Plot the fraction of total
	# reads with A/G/T/C in their terminals

	terminal_basecount = {
		"totalreads": 0,
		"5": {
			"A": 0,
			"G": 0,
			"T": 0,
			"C": 0,
			"N": 0
		},
		"3": {
			"A": 0,
			"G": 0,
			"T": 0,
			"C": 0,
			"N": 0
		}
	}

	bamfile    = global_config["bamfile"]

	for gene_name, details in annotation.items():

		start  = details["start"]
		stop   = details["stop"]
		strand = details["strand"]
		refer  = details["refname"]

		for read in bamfile.fetch(start=start, stop=stop, reference=refer):

			if strand == "+" and read.is_reverse:
				continue
			if strand == "-" and not read.is_reverse:
				continue
			if read.reference_length > 50:
				continue

			seq  = read.query_sequence

			if strand == "+":
				fp_base = seq[0]
				tp_base = seq[len(seq)-1]
				terminal_basecount["5"][fp_base] += 1
				terminal_basecount["3"][tp_base] += 1
			if strand == "-":
				fp_base = seq[len(seq)-1]
				tp_base = seq[0]
				terminal_basecount["5"][fp_base] += 1
				terminal_basecount["3"][tp_base] += 1

			terminal_basecount["totalreads"] += 1

	save_file(terminal_basecount, "terminal_base_fraction")

	xaxis = ["A", "T", "G", "C"]

	# Calculate the fraction of total
	A5 = terminal_basecount["5"]["A"]/terminal_basecount["totalreads"]
	G5 = terminal_basecount["5"]["G"]/terminal_basecount["totalreads"]
	T5 = terminal_basecount["5"]["T"]/terminal_basecount["totalreads"]
	C5 = terminal_basecount["5"]["C"]/terminal_basecount["totalreads"]

	A3 = terminal_basecount["3"]["A"]/terminal_basecount["totalreads"]
	G3 = terminal_basecount["3"]["G"]/terminal_basecount["totalreads"]
	T3 = terminal_basecount["3"]["T"]/terminal_basecount["totalreads"]
	C3 = terminal_basecount["3"]["C"]/terminal_basecount["totalreads"]

	y3 = [A3, T3, G3, C3]
	y5 = [A5, T5, G5, C5]


	# Plot for five prime end
	plt.bar(xaxis, y5)
	plt.xlabel("Base at read terminal")
	plt.ylabel("Read fraction")
	plt.title("Five prime end ({}) ({})".format(global_config["prefix"],global_config["offset_for_site"]))
	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"5p_basestat",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()

	# Plot for three prime end
	plt.bar(xaxis, y3)
	plt.xlabel("Base at read terminal")
	plt.ylabel("Read fraction")
	plt.title("Three prime end ({}) ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"3p_basestat",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()



def calculate_asymmetry_scores():
	'''
	Asymmetry score per gene =

		log2(sum(densities_2nd_half_of_gene)/sum(densities_1st_half))
	'''
	def gethalf(vector):
		
		pivot = int(np.floor(len(vector)/2))
		return(vector[:pivot], vector[pivot:])

	# Dictionary containing just the asymmetry score per gene
	asymmetry = {}
	tmp_asc   = [] # temporarily store scores in a list for box plot

	for gene_name, densities in per_gene_densities.items():

		if len(densities) < 100:
			notify("Asymmetry: {} too short. Skipped.".format(
				gene_name,
				level="warn"
			))
			continue

		# Ignore the first and last 50nt from the start and ends of the
		# gene to remove any biases.
		density_5p, density_3p = gethalf(densities[50:len(densities)-50])
		rho_sum_5p             = sum(density_5p)
		rho_sum_3p             = sum(density_3p)

		if ((rho_sum_3p == 0) or (rho_sum_5p == 0)):
			notify("Asymmetry: No density profile for {}. Skipped.".format(gene_name),
					level="warn")
			continue

		asymmetry_score      = np.log2(rho_sum_3p/rho_sum_5p)
		asymmetry[gene_name] = asymmetry_score
		tmp_asc.append(asymmetry_score)

	save_file(asymmetry, "asymmetry_scores")

	plt.boxplot(tmp_asc, vert=False, notch=True)
	plt.xlabel("log2(sum(#rpfs 3' half)/sum(#rpfs 5' half))")
	plt.title("Asymmetry {} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))
	plt.axvline(x=0, color="red", linestyle="--")
	plt.tight_layout()

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"asymmetry_score",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()



def consistency_score_per_transcript():
	'''
	Reads mapping to the second half of a gene vs the first
	NOTE: Taking the whole gene, not trimming anything
	'''
	bamfile = global_config["bamfile"]
	five_prime_counts  = []
	three_prime_counts = []

	consistency_dict = {} # temporarily hold the data before dumping

	for gene_name, details in annotation.items():

		count_5p, count_3p = 0, 0

		start  = details["start"]
		stop   = details["stop"]
		strand = details["strand"]
		refer  = details["refname"]

		pivot  = int(np.ceil((stop-start)/2)) + start

		for read in bamfile.fetch(reference=refer, start=start, stop=pivot):

			if strand == "+" and read.is_reverse:
				continue
			if strand == "-" and not read.is_reverse:
				continue

			count_5p += 1


		for read in bamfile.fetch(reference=refer, start=pivot, stop=stop):

			if strand == "+" and read.is_reverse:
				continue
			if strand == "-" and not read.is_reverse:
				continue

			count_3p += 1

		if ((count_5p == 0) or (count_3p == 0)):
			continue

		five_prime_counts.append(count_5p)
		three_prime_counts.append(count_3p)

		consistency_dict[gene_name] = {
			"5_prime_half": count_5p,
			"3_prime_half": count_3p
		}

	five_prime_counts  = np.log2(np.array(five_prime_counts)) # y
	three_prime_counts = np.log2(np.array(three_prime_counts)) # x
	pearson = np.corrcoef(three_prime_counts,five_prime_counts)[1, 0]

	m, c = np.polyfit (
		three_prime_counts,
		five_prime_counts,
		1 # 1st degree polynomial fit
	)

	consistency_dict["metadata"] =  {
		"pearson_coeff": pearson,
		"slope": m,
		"intercept": c
	}

	save_file(consistency_dict, "consistency data")

	regression_function = np.poly1d(m, c)

	plt.plot(
		three_prime_counts, # x
		five_prime_counts, # y
		"bo",
		three_prime_counts,   # x
		regression_function(three_prime_counts), # fx
		"-k",
		markersize=1,
		alpha=0.5
	)

	plt.title("y = ({})x + {}, Pearson's R = {}".format(
		round(m,2), round(c,2),
		round(pearson, 2)
	))

	plt.xlabel("log2(sum(read count 3' half))")
	plt.ylabel("log2(sum(read count 5' half))")
	plt.suptitle("{} ({})".format(global_config["prefix"],	global_config["offset_for_site"]))

	plt.savefig("{}{}_{}_{}_{}.tif".format(global_config["img_dir"],
		global_config["session_id"],
		"consistency_plot",
		global_config["img_id"],
		global_config["offset_for_site"]))

	plt.close()


def construct_metagene_vector():

	up = global_config["endmap_upstream"]
	dn = global_config["endmap_downstream"]
	lengths = global_config["use_readlengths"]

	tmp_list = []

	for length in lengths:
		for gene_name in annotation.keys():

			try:
				# chop out unwanted values
				vector = per_gene_endmaps[gene_name][length][:(up+dn)]
			except:
				# Gene doesn't have RPFs of a particular length
				continue

			tmp_list.append(normalize_counts(vector))

		metagene_for_len[length] = np.apply_along_axis (
			sum, 0, np.array(tmp_list)
		)

		del tmp_list[:]


def retrieve_coding_sequence(fname, ftype="fasta"):
	if len(transcripts_dict) == 0 :
		cds =  list(SeqIO.parse(fname, ftype))
		for i in range(len(cds)):
			transcripts_dict[cds[i].id] = str(cds[i].seq)


def calculate_densities_over_genes():
	"""
	Takes in the offset applied per_gene_endmap vector and does a column wise
	summation for all read lengths, basically giving us the metagene profile
	for each gene in question, returning a new dictionary with gene name -->
	mapping to it's metagene profile (for all read lengths), with offsets applied.

	WARNING: Do not run this function unless the offsets have been applied
	"""

	up = global_config["endmap_upstream"]
	# #nt downstream CDS start (yes) to display in combined metagene
	dn = global_config["endmap_downstream"]

	lengths = global_config["use_readlengths"]

	skipped = []
	notify("{} genes in per_gene_endmaps".format(len(per_gene_endmaps.keys())))
	for gene, endmaps_per_readlength in per_gene_endmaps.items():
		# WARNING: These values are not normalized by any means.
		# Normalize for length before any calculations can be done.
		tmp = [] # store the endmaps per readlength.
		for length, endmap in endmaps_per_readlength.items():
			if length not in lengths:
				continue # only calculate for the particular set of RPFs
			tmp.append(endmap)

		if len(tmp) == 0:
			notify("{}: no coverage data, can't determine densities".format(gene), level="warn")
			continue

		# print("{}, {}, {}".format(gene, len(tmp), tmp[1:5]))

		tmp = np.array(tmp)
		tmp = np.apply_along_axis(sum, 0, tmp)
		# trimming out the flanking regions we added 
		tmp = tmp[up:len(tmp)-dn]
		per_gene_densities[gene] = tmp.tolist().copy()
		tmp = tmp.tolist()
		del tmp[:]
	notify("{} genes in per_gene_densities".format(len(per_gene_densities)))
	save_file(unique(skipped), "SKIPPED_GENES_NO_COVERAGE_CALCDENSITY.json")


def main():

	print("\n----------------------------------------------------------------")
	print("               {}\n".format(versionstr))
	print(" This program is distributed in the hope that it will be useful,")
	print(" but WITHOUT ANY WARRANTY; without even the implied warranty of")
	print("  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the")
	print("         GNU General Public License for more details.")
	print("----------------------------------------------------------------\n")

	try:
		set_config(load_json(sys.argv[1]))
	except:
		# Check if a dummy config file already exists in the directory.
		# If not we create a new one.
		if not os.path.exists("config_template.json"):
			fh = open("config_template.json", "w")
			ft = open("config_template_help.txt", "w")
			fh.write(CONFIG_TMP)
			ft.write(CONFIG_HELP)
			fh.close()
			ft.close()

		print(" Malformed or non-existant configuration file\n")
		print(" Try using the templates if you need help\n")
		print(" Bye!\n")

		raise SystemExit()

	print("\nInitializing ...")
	script_init()

	# notify("Assuming BAM file has single end reads, with all unmapped reads removed",
	# 		level="warn")
	# notify("Assuming all reads map perfectly at their 5' or 3' ends", level="warn")
	# notify("Be sure to preprocess your BAM file with, eg: samtools, to ensure the above")
	check_gene_list()


	notify("Determining terminal base fractions")
	read_terminal_stats()


	notify("Found {} genes in annotation file".format(len(annotation)))


	notify("Plotting read length histogram")
	read_length_histogram()

	notify("Filtering top {}% genes by coverage".format(
		global_config["filtering_threshold"])
	)

	get_max_coverage_genes(percentage=global_config["filtering_threshold"])
	notify("Reduced number of genes to {}".format((len(high_coverage_genes))))


	notify("Mapping genes to endmaps (no offset)")
	map_gene_to_endmaps(apply_offset=False)



	notify("Constructing metagene vector")
	construct_metagene_vector() # yeast



	notify("Plotting metagene")
	plot_metagene_per_readlength()



	notify("Determining offsets")
	offset_delta_for_metagene()



	notify("Mapping genes to endmaps (with offsets)")
	map_gene_to_endmaps(apply_offset=True, only_reads_with_codons = {
			'E': [],
			'P': [],
			'A': []
		}
	)

	notify("Determining ribosome densities over genes")
	calculate_densities_over_genes()


	notify("Calculating pause scores")
	calculate_pause_score()


	notify("Plotting pause scores")
	plot_pauses()

	notify("Plotting stop peak")
	plot_stop_peak()

	notify("Calculating asymmetry scores")
	calculate_asymmetry_scores()

	notify("Consistency score per transcript is disabled")
	# consistency_score_per_transcript()


	notify("Metagene profile over codons is disabled")
	# metagene_over_codon(aa_to_codon["A"], mode="combined") # serine
	# metagene_over_codon([aa_to_codon["S"], aa_to_codon["A"]], mode="composite")
	try:
		cm_mode   = global_config["codon_metagene_mode"]
		codonlist = global_config["codon_metagene_list"]
		metagene_over_codon(expand_to_codon(codonlist), mode=cm_mode)
	except:
		notify("No metagene for codons specified, skipping", level="warn")


	notify("Determining reading frame")
	plot_reading_frame(reading_frame(ignore_read_by_nt={5:[], 3:[]}))


	notify("Constructing metagene vector")
	construct_metagene_vector()


	save_file(per_gene_densities, "per_gene_densities")

	notify("Plotting metagene")
	plot_metagene_per_readlength()


	notify("Plotting combined metagene")
	plot_combined_metagene(peak=True, peak_range=[5,25])

	notify("Coverage profile disabled")
	generate_coverage_profiles(disabled=True)

	notify("Cleaning up")
	global_config["bamfile"].close()
	del global_config["bamfile"] # can't be saved as plain text
	# Save the configuration for reference
	save_json(global_config, global_config["conf_file"])

	notify("Done!")


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		global_config["bamfile"].close()
		notify("Keyboard interrupt. Aborting.", level="crit", fatal=True)

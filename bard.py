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
# bard is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# bard is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses>
#
#
#
# Usage:
# ~$ python bard.py <JSON config file>
#
# ---------------------------------------------------------------------------
print("Starting, please wait ...", end="\r")
from collections import defaultdict
from datetime import datetime as dt
from Bio import SeqIO, SeqUtils
import seaborn as sns
import numpy as np
import inspect
import random
import pysam
import json
import sys
import os


import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
colors = plt.rcParams["axes.prop_cycle"]()

np.random.seed(12345)

# This line is automatically updated before each commit
# Do not edit
versionstr = "bard v1.0 ID=14-29-15-27-03-2020"

codon_to_aa = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}

codon_counts = {
    "TTT": 0,
    "TTC": 0,
    "TTA": 0,
    "TTG": 0,
    "CTT": 0,
    "CTC": 0,
    "CTA": 0,
    "CTG": 0,
    "ATT": 0,
    "ATC": 0,
    "ATA": 0,
    "ATG": 0,
    "GTT": 0,
    "GTC": 0,
    "GTA": 0,
    "GTG": 0,
    "TAT": 0,
    "TAC": 0,
    "TAA": 0,
    "TAG": 0,
    "CAT": 0,
    "CAC": 0,
    "CAA": 0,
    "CAG": 0,
    "AAT": 0,
    "AAC": 0,
    "AAA": 0,
    "AAG": 0,
    "GAT": 0,
    "GAC": 0,
    "GAA": 0,
    "GAG": 0,
    "TCT": 0,
    "TCC": 0,
    "TCA": 0,
    "TCG": 0,
    "CCT": 0,
    "CCC": 0,
    "CCA": 0,
    "CCG": 0,
    "ACT": 0,
    "ACC": 0,
    "ACA": 0,
    "ACG": 0,
    "GCT": 0,
    "GCC": 0,
    "GCA": 0,
    "GCG": 0,
    "TGT": 0,
    "TGC": 0,
    "TGA": 0,
    "TGG": 0,
    "CGT": 0,
    "CGC": 0,
    "CGA": 0,
    "CGG": 0,
    "AGT": 0,
    "AGC": 0,
    "AGA": 0,
    "AGG": 0,
    "GGT": 0,
    "GGC": 0,
    "GGA": 0,
    "GGG": 0,
}

# List of pause score values per codon
codon_pauselist = {
    "E": {
        "ATA": [],
        "ATC": [],
        "ATT": [],
        "ATG": [],
        "ACA": [],
        "ACC": [],
        "ACG": [],
        "ACT": [],
        "AAC": [],
        "AAT": [],
        "AAA": [],
        "AAG": [],
        "AGC": [],
        "AGT": [],
        "AGA": [],
        "AGG": [],
        "CTA": [],
        "CTC": [],
        "CTG": [],
        "CTT": [],
        "CCA": [],
        "CCC": [],
        "CCG": [],
        "CCT": [],
        "CAC": [],
        "CAT": [],
        "CAA": [],
        "CAG": [],
        "CGA": [],
        "CGC": [],
        "CGG": [],
        "CGT": [],
        "GTA": [],
        "GTC": [],
        "GTG": [],
        "GTT": [],
        "GCA": [],
        "GCC": [],
        "GCG": [],
        "GCT": [],
        "GAC": [],
        "GAT": [],
        "GAA": [],
        "GAG": [],
        "GGA": [],
        "GGC": [],
        "GGG": [],
        "GGT": [],
        "TCA": [],
        "TCC": [],
        "TCG": [],
        "TCT": [],
        "TTC": [],
        "TTT": [],
        "TTA": [],
        "TTG": [],
        "TAC": [],
        "TAT": [],
        "TAA": [],
        "TAG": [],
        "TGC": [],
        "TGT": [],
        "TGA": [],
        "TGG": [],
    },
    "P": {
        "ATA": [],
        "ATC": [],
        "ATT": [],
        "ATG": [],
        "ACA": [],
        "ACC": [],
        "ACG": [],
        "ACT": [],
        "AAC": [],
        "AAT": [],
        "AAA": [],
        "AAG": [],
        "AGC": [],
        "AGT": [],
        "AGA": [],
        "AGG": [],
        "CTA": [],
        "CTC": [],
        "CTG": [],
        "CTT": [],
        "CCA": [],
        "CCC": [],
        "CCG": [],
        "CCT": [],
        "CAC": [],
        "CAT": [],
        "CAA": [],
        "CAG": [],
        "CGA": [],
        "CGC": [],
        "CGG": [],
        "CGT": [],
        "GTA": [],
        "GTC": [],
        "GTG": [],
        "GTT": [],
        "GCA": [],
        "GCC": [],
        "GCG": [],
        "GCT": [],
        "GAC": [],
        "GAT": [],
        "GAA": [],
        "GAG": [],
        "GGA": [],
        "GGC": [],
        "GGG": [],
        "GGT": [],
        "TCA": [],
        "TCC": [],
        "TCG": [],
        "TCT": [],
        "TTC": [],
        "TTT": [],
        "TTA": [],
        "TTG": [],
        "TAC": [],
        "TAT": [],
        "TAA": [],
        "TAG": [],
        "TGC": [],
        "TGT": [],
        "TGA": [],
        "TGG": [],
    },
    "A": {
        "ATA": [],
        "ATC": [],
        "ATT": [],
        "ATG": [],
        "ACA": [],
        "ACC": [],
        "ACG": [],
        "ACT": [],
        "AAC": [],
        "AAT": [],
        "AAA": [],
        "AAG": [],
        "AGC": [],
        "AGT": [],
        "AGA": [],
        "AGG": [],
        "CTA": [],
        "CTC": [],
        "CTG": [],
        "CTT": [],
        "CCA": [],
        "CCC": [],
        "CCG": [],
        "CCT": [],
        "CAC": [],
        "CAT": [],
        "CAA": [],
        "CAG": [],
        "CGA": [],
        "CGC": [],
        "CGG": [],
        "CGT": [],
        "GTA": [],
        "GTC": [],
        "GTG": [],
        "GTT": [],
        "GCA": [],
        "GCC": [],
        "GCG": [],
        "GCT": [],
        "GAC": [],
        "GAT": [],
        "GAA": [],
        "GAG": [],
        "GGA": [],
        "GGC": [],
        "GGG": [],
        "GGT": [],
        "TCA": [],
        "TCC": [],
        "TCG": [],
        "TCT": [],
        "TTC": [],
        "TTT": [],
        "TTA": [],
        "TTG": [],
        "TAC": [],
        "TAT": [],
        "TAA": [],
        "TAG": [],
        "TGC": [],
        "TGT": [],
        "TGA": [],
        "TGG": [],
    },
}

# Amino acid to list of synonymous codons
aa_to_codon = {
    "I": ["ATA", "ATC", "ATT"],
    "M": ["ATG"],
    "T": ["ACA", "ACC", "ACG", "ACT"],
    "N": ["AAC", "AAT"],
    "K": ["AAA", "AAG"],
    "S": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
    "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
    "L": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
    "P": ["CCA", "CCC", "CCG", "CCT"],
    "H": ["CAC", "CAT"],
    "Q": ["CAA", "CAG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "D": ["GAC", "GAT"],
    "E": ["GAA", "GAG"],
    "G": ["GGA", "GGC", "GGG", "GGT"],
    "F": ["TTC", "TTT"],
    "Y": ["TAC", "TAT"],
    "C": ["TGC", "TGT"],
    "W": ["TGG"],
    "*": ["TAA", "TAG", "TGA"],
}

# Holds the global configuration
global_config = {}
# GFF
annotation = {}
# Keep a track of the warnings emitted
__warnings = []
# gene name --> read length -->
# list of read alignment coordinates
endmap_vectors = {}
endmap_vectors["E"] = {}
endmap_vectors["P"] = {}
endmap_vectors["A"] = {}
# Metagene vectors for each read length
metagene_per_readlength = {}
# Names of genes overlapping by either
# their 5' or 3' ends
overlap_dict = {}
# gene name --> coverage value
gene_coverages = {}
# Genes with coverage above threshold
high_coverage_genes = []
# readlength --> signed integer offset
offsets_e = {}  # P-site
offsets_p = {}  # E-site
offsets_a = {}  # A-site
# gene name --> ribosome density vector
ribosome_dentities = {}
ribosome_dentities["E"] = {}
ribosome_dentities["P"] = {}
ribosome_dentities["A"] = {}
# gene name --> cds sequence
transcripts_dict = {}
# gene name --> cds sequence in codons &
# the positional pause vector
pps_vector_per_gene = {}
pps_vector_per_gene["E"] = {}
pps_vector_per_gene["P"] = {}
pps_vector_per_gene["A"] = {}
# gene name --> read coverage profile vector
coverage_profile = {}


# Template configuration, to be printed out when
# no input is provided.
CONFIG_TMP = """{
    "coding_sequence_path": "/home/user/path/to/cds.fa",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "/home/user/path/to/annotation.gff",
    "annotation_feature_tag": "ID",
    "bam_file_path": "/home/user/path/to/alignment.bam",
    "read_offset_terminal": "five_prime",
    "coverage_cutoff": 40,
    "coverage_metric": "reads_per_nt",
    "will_ignore_overlaps": true,
    "peak_scan_range": [-25, -5],
    "use_readlengths": [26,27,28,29,30,31,32,33],
    "gene_list_file": "/home/user/path/to/gene_list.txt",
    "gene_list_action": "include_only",
    "genes_overlap_exception": "/home/user/path/to/overlapping_genes_list.txt"
}
"""
# Associated help text.
CONFIG_HELP = """
#########################################################################################
############################## bard v1.0 quick start guide ##############################
#########################################################################################

In order to run bard, you'll have to supply it with a JSON configuration file, like so:


        ~$ python bard.py config.json


A JSON file is just a structured text file you can open and edit in your text editor.
It looks something like this:



        {
            "option1": value1,
            "option2": value2,
            "option3": value3               (.. and so on)
        }



Notice that the last option-value pair of a JSON file does not end with a comma. Also,
while the "options" are strings, the "values" can be of any type, including integers
floats, arrays, boolean, etc.


An empty configuration file (bard_config_template.json) has been saved to your working
directory in case you need it.


Wherever necessary, please provide the absolute paths. For example, instead of writing
"file.txt", use "/home/user/full/path/to/file.txt" instead.


The available options and possible values are elaborated below. The mandatory ones are
marked with a star(*):


--------------------------------------+--------------------------------------------------
            OPTIONS                   |                     VALUES
--------------------------------------+--------------------------------------------------
*1. coding_sequence_path (-c):        |  Absolute path to the CDS / cDNA file, containing
                                      |  the sequences for all genes. Make sure that the
                                      |  fasta headers exactly match the names in the GFF
--------------------------------------+--------------------------------------------------
*2. coding_sequence_format (-cf):     |  The format of the CDS/cDNA file. eg: "fasta"
--------------------------------------+--------------------------------------------------
*3. annotation_file_path (-g):        |  Absolute path to the annotation file (GTF/GFF)
--------------------------------------+--------------------------------------------------
*4. annotation_feature_tag (-f):      |  The tag in column 9 of the GFF/GTF file which
                                      |  uniquely identifies the feature (gene).
                                      |  (eg: "gene_id" or "ID" or "protein_id")
--------------------------------------+--------------------------------------------------
*5. bam_file_path (-b):               |  Absolute path to the alignment file (BAM)
--------------------------------------+--------------------------------------------------
*6. read_offset_terminal (-t):        |  Calculate the P/E/A site offsets for either
                                      |  "five_prime" or "three_prime" ends of the reads
--------------------------------------+--------------------------------------------------
*7. will_ignore_overlaps (-x):        |  Ignore overlapping genes? true or false(boolean)
--------------------------------------+--------------------------------------------------
*8. peak_scan_range (-s):             |  Scan for the initiation peak within this range
                                      |  of nucleotides w.r.t the start codon.
                                      |  For example, [-25, -5] checks for the initiation
                                      |  peak -25 to -5 nt upstream of the start codon,
                                      |  in case you are using 5' offsets. If you are
                                      |  using 3' offsets, you may want to use something
                                      |  like: [5, 30]. You'll need to figure out the
                                      |  exact values empirically.
--------------------------------------+--------------------------------------------------
9. use_readlengths (-r):                   |  Ribosome protected fragment (RPF/read) lengths
                                      |  to use in the analysis. By default, it looks at
                                      |  all read lengths from 20 to 40nt, and then uses
                                      |  some heuristics to figure out which ones to use.
--------------------------------------+--------------------------------------------------
10. gene_list_file (-gf):             |  Absolute path to a file containing a list of
                                      |  gene names, one per line. The names must match
                                      |  with those in the GFF and the CDS/cDNA multiple
                                      |  fasta file.
--------------------------------------+--------------------------------------------------
11. gene_list_action (-ga):           |  Specifies what to do with the list of genes.
                                      | -------------------
                                      |
                                      |  1."include_only":    Run the analysis only for
                                      |                       the genes mentioned in the
                                      |                       file
                                      |
                                      |  2."exclude_only":    Use all genes except the
                                      |                       ones mentioned in the file
                                      |
                                      |  3."exclude_balance": Don't use the genes in the
                                      |                       file, instead, from the
                                      |                       remaining, use an equal
                                      |                       number of randomly selected
                                      |                       genes.
--------------------------------------+--------------------------------------------------
12. genes_overlap_exception (-oe):    |  Absolute path to a file containing a list of
                                      |  gene names (all of which overlap with some other)
                                      |  which will NOT be ignored during the analysis.
                                      |  One gene name per line.
--------------------------------------+--------------------------------------------------
13. coverage_cutoff (-cv)             |  Exclude any gene which has a coverage value
                                      |  below this numerical threshold. Defaults to 10.
--------------------------------------+--------------------------------------------------
14. coverage_metric (-cm):            |  Specifies how the coverage value is calculated.
                                      |  Can be either "reads_per_nt" (on average) or
                                      |  "rpkm". Defaults to RPKM.
--------------------------------------+--------------------------------------------------
"""

unique = lambda v: list(set(v))
stats = lambda v: {"avg": round(np.mean(v), 3), "std": round(np.std(v), 3)}


def line_number():
    """
    Calling this function from a particular line
    in our program will return the corresponding
    line number. Useful for debugging.
    """
    return inspect.currentframe().f_back.f_lineno


def load_json(file="object.json"):
    """
    Load a JSON file
    """
    with open(file, "r") as read_file:
        return json.load(read_file)


def save_json(data, file="object.json"):
    """
    Save object as JSON
    """
    with open(file, "w") as write_file:
        json.dump(data, write_file, indent=4)


def extended_sort(x, by, reverse=False):
    """
    Takes two lists and sorts them based on the
    values in one of them (by)
    """
    by, x = zip(*sorted(zip(by, x), reverse=reverse))
    return (by, x)


def sort_dict(d, by="key", reverse=False):
    """
    Sorts a dictionary by either the key
    or the value
    """
    if by == "key":
        mode = 0
    if by == "value":
        mode = 1

    tmp = sorted(d.items(), key=lambda item: item[mode], reverse=reverse)
    td = {}
    for item in tmp:
        td[item[0]] = item[1]

    return td


def subset_integer_progression(v, step=1, continuity="none"):
    """
    Takes in a dictionary and returns the subset whose
    values show integer progression.

    Used to select the best read lengths for E/P/A site
    offset calculation.

    continuity can be:
        "increasing", "decreaing" or "none"
    """
    ip, ret = [], {}
    dsize = len(v)
    keys = list(v.keys())
    vals = list(v.values())

    for i in range(dsize):
        if i + 1 >= dsize:
            break
        # using absolute values here, because
        # the five prime offsets return negative
        if continuity == "none":
            delta = abs(vals[i] - vals[i + 1])
        if continuity == "increasing":
            delta = abs(vals[i + 1]) - abs(vals[i])
        if continuity == "decreasing":
            delta = abs(vals[i]) - abs(vals[i + 1])
        # 3' end can have same length cuts
        # Lets's see how this goes
        if delta == step or delta == 0:
            ip.append(keys[i])
            ip.append(keys[i + 1])

    ip = list(set(ip))

    # subset
    for k, x in v.items():
        if k not in ip:
            continue
        ret[k] = x

    return ret


def expand_to_codon(v):
    """
    Expand amino acid letters in a list to codons.
    eg: ["C", "ATG"] becomes ["TGT", "TGC", "ATG"]
    """
    valid_codon = list(codon_to_aa.keys())
    valid_amino = list(codon_to_aa.values())
    tmp = []

    def expanded(j):
        if j in valid_codon:
            return [j]
        if j in valid_amino:
            return aa_to_codon[j]

    def invalid(v):
        if (v in valid_codon) or (v in valid_amino):
            return False
        else:
            return True

    for i in v:
        if isinstance(i, list):
            tmp_1 = []
            for j in i:
                if invalid(j):
                    continue
                for c in expanded(j):
                    tmp_1.append(c)
            tmp.append(tmp_1)
        else:
            if invalid(i):
                continue
            for c in expanded(i):
                tmp.append(c)
    return tmp


def extract_file_name(path):
    """
    Given an absolute path, extracts the file name
    (without extension)
    NOTE: May fail on Windows
    """
    return os.path.split(path)[1].split(".")[0]


def reset():
    """
    Purge global environment
    """
    global_config.clear()
    gene_coverages.clear()
    annotation.clear()
    __warnings.clear()
    endmap_vectors.clear()
    metagene_per_len.clear()
    overlap_dict.clear()
    ribosome_dentities.clear()
    transcripts_dict.clear()
    pps_vector_per_gene.clear()
    offsets_p.clear()
    offsets_e.clear()
    offsets_a.clear()


def groupN(x, n):
    """
    Group list (x) contents into n members each
    eg: x=[1,2,3,4,5,6], n=3, becomes [[1,2,3],[4,5,6]]
        x=[1,2,3,4]    , n=3, becomes [[1,2,3],[4]]
    """
    if n <= 0:
        return
    tmp = []
    for i in range(0, len(x) + n, n):
        j = i + n
        vec = x[i:j]
        if len(vec) != 0:
            tmp.append(vec)
        if j == len(x):
            break
    return tmp


def close_enough(a, b, thresh=0.1):
    """
    Returns whether two numbers are equivalent
    within a certain tolerance.
    """
    return abs(a - b) <= thresh


def ungroup(x):
    """
    Ungroup a 1-level deep nested list
    """
    tmp = []
    for i in x:
        for j in i:
            tmp.append(j)
    return tmp


def position_of_codon(string, motifs, mode="codon"):
    """
    Returns the zero-indexed position of
    a codon in a given ORF.
    """
    if len(string) % 3 != 0:
        return "err"
    codons = groupN(string, 3)
    indexes = []
    i = 0
    for motif in motifs:
        for codon in codons:
            if codon == motif:
                if mode == "nt":
                    indexes.append(i * 3)
                if mode == "codon":
                    indexes.append(i)
            i += 1
        i = 0
    return indexes


def random_hexadecimal_string():
    """
    Random number prepended to filenames
    to avoid conflicts
    """
    return "%06x" % random.randrange(16 * 1e6)


def value_from_dict(dict, key):
    """
    Return 0 if a dictionary lacks a key, do not crash
    """
    try:
        return dict[key]
    except KeyError:
        return 0


def notify(event, level="notf", onetime=False, fatal=False):
    """
    Runtime logs
    """

    if onetime:  # issue a notification only once (in a loop)
        if event not in __warnings:
            __warnings.append(event)
            pass
        else:
            return

    timestamp = str(dt.now().strftime("%H.%M.%S"))

    # Open in append mode
    log_path = open(global_config["log_file"], "a")

    log_path.write("[{}] {}\n".format(timestamp, event))
    log_path.close()

    RED = "\033[31m"
    GRN = "\033[32m"
    YLW = "\033[33m"
    RST = "\033[m"

    if level == "notf":
        print("{}[NOTF: {}] {}{}".format(GRN, timestamp, RST, event))

    if level == "warn":
        print("{}[WARN: {}] {}{}".format(YLW, timestamp, RST, event))

    if level == "crit":
        print("{}[CRIT: {}] {}{}".format(RED, timestamp, RST, event))
        if fatal:
            # global_config["bamfile"].close()
            raise SystemExit()


def set_global_config(**kwargs):
    """
    Used by init. Allows specifying
    the key-value pairs within the
    function call itself
    """
    for key, value in kwargs.items():
        global_config[key] = value


def set_config(config):
    """
    Set the global configuration from a dictionary.
    This is used for the JSON config file.
    """
    for key, value in config.items():
        global_config[key] = value


def get_offset(xvec, yvec, right_bound, left_bound):
    """
    Returns the P-site offset in a metagene vector
    xvec: A vector of indices (coordinates)
    yvec: Metagene vector
    right_bound, left_bound: lower and upperbound
    coordinate ranges in which to search for the
    metagene maxima.
    """

    if isinstance(xvec, np.ndarray):
        xvec = xvec.tolist()
    if isinstance(yvec, np.ndarray):
        yvec = yvec.tolist()

    left = xvec.index(left_bound)
    right = xvec.index(right_bound)
    y_trim = yvec[left : right + 1]  # because python
    pos = yvec.index(max(y_trim))
    offset = pos - xvec.index(0)
    return (xvec[pos], offset, max(y_trim))


def offset_delta_for_metagene():
    """
    Get psite offsets from the metagene profiles of each RPF
    readlength
    """

    global offsets_p
    global offsets_e
    global offsets_a

    left_bound = global_config["peak_scan_range"][0]
    right_bound = global_config["peak_scan_range"][1]
    progression = global_config["offset_progression"]

    x_axis = metagene_per_readlength["xaxis"]
    peak_heights, n = {}, 1  # readlength --> peak_height

    for rpf_length, metagene_vector in metagene_per_readlength.items():

        if len(metagene_vector) == 1 and metagene_vector[0] == -1:
            notify(
                "Not calculating offset for RPF length {}".format(rpf_length),
                level="warn",
            )
            continue

        if rpf_length == "xaxis":
            continue

        _, delta, maxima = get_offset(
            xvec=x_axis,
            yvec=metagene_vector,
            left_bound=left_bound,
            right_bound=right_bound,
        )

        peak_heights[rpf_length] = maxima

        offsets_p[rpf_length] = delta
        # The 5' offsets returned above are negative.
        # The expression will be the same irrespective
        # of whether we're taking 5' or 3' offsets.
        offsets_e[rpf_length] = delta + 3
        offsets_a[rpf_length] = delta - 3

    # First we select readlengths by integer progression of their offsets.
    # Then filter for profile quality by selecting fragment lengths within
    # mean+/-stdev of peak heights

    if len(global_config["use_readlengths"]) == 0:
        # We select the best readlengths
        def subset_dict(off_map, ph_map, vals_lst):
            """
            Returns the subset of a dict whose
            values() are roughly equal to values
            provided in list
            """
            tmp_map = {}
            for i in vals_lst:
                for key, value in ph_map.items():
                    if close_enough(i, value, thresh=0.1):
                        try:
                            tmp_map[key] = off_map[key]
                        except:
                            # we have subset the offset dict.
                            # some keys will be unavailable
                            continue

            return tmp_map

        ph = np.array(list(peak_heights.values()))

        # Arbitrary empirical cutoffs. Just to reduce noise
        peaks_select = ph[(ph >= 0) & (ph >= (max(ph) - ((60 / 100) * max(ph))))]

        offsets_p = sort_dict(
            subset_dict(offsets_p, peak_heights, peaks_select), by="key", reverse=False
        )

        offsets_p = subset_integer_progression(offsets_p, continuity=progression)

        # Select readlengths with metagene peak maxima within
        # +/- N stdevs from mean

        # May be string or int. may cause bugs
        global_config["readlengths"] = list(offsets_p.keys())

        offsets_e = sort_dict(
            subset_dict(offsets_e, peak_heights, peaks_select), by="key", reverse=False
        )

        offsets_e = subset_integer_progression(offsets_e, continuity=progression)

        offsets_a = sort_dict(
            subset_dict(offsets_a, peak_heights, peaks_select), by="key", reverse=False
        )

        offsets_a = subset_integer_progression(offsets_a, continuity=progression)

    save_file(offsets_p, "psite_offsets")
    save_file(offsets_e, "esite_offsets")
    save_file(offsets_a, "asite_offsets")


def read_genes_from_file(filename):
    """
    Given a file with one gene per line,
    return the genes as a list
    """
    try:
        fh = open(filename, "r")
    except:
        raise FileNotFoundError

    genes = []
    for gene in fh:
        genes.append(gene.strip("\n").strip(" "))

    return genes


def parse_gff():
    """
    Converts a GFFv3/GTF file into a dictionary.
    """
    annofile = global_config["annotation_file_path"]
    # Feature tag in GFF column 9 which contains gene name
    feature_tag = global_config["annotation_feature_tag"]

    annolist = []
    for record in open(annofile, "r"):
        record = record.strip("\n")
        record = record.split("\t")
        if (len(record) < 9) and ("#" in record[0]):
            # comment line. ignore
            continue
        if len(record) != 9:
            # GFF/GTF should only have 9 columns
            continue
        annolist.append(record)

    all_gene_names = []
    # We can sometimes find duplicate gene names in the annotation files.
    # In case a duplicate is detected, this function appends an integer
    # suffix to it's name, which is a count of #times it was encountered.
    def check_duplicate_name(gene_name):

        if gene_name in all_gene_names:
            notify(
                "Found duplicate annotation. Will tag and ignore.",
                level="warn",
                onetime=True,
            )
            all_gene_names.append(gene_name)

            occurances = str(all_gene_names.count(gene_name) - 1)
            gene_name = gene_name + "_" + occurances
        else:
            all_gene_names.append(gene_name)

        return gene_name

    def extract_gene_name(index):
        found_tag = False
        for attribute in annolist[index][8].split(";"):
            # strip the feature tag and return gene name
            if feature_tag in attribute:
                found_tag = True
                return check_duplicate_name(
                    attribute.replace(feature_tag, "")
                    .replace('"', "")
                    .replace("=", "")
                    .strip()
                )
            else:
                pass

        if not found_tag:
            # No point continuing
            notify("{}th gene lacks feature tag".format(i), level="crit", fatal=True)

    # Populate the annotation dictionary
    for i in range(len(annolist)):
        annotation[extract_gene_name(i)] = {
            "start": int(annolist[i][3]),
            "stop": int(annolist[i][4]),
            "strand": annolist[i][6],
            "refname": annolist[i][0],
        }

    save_file(annotation, "Annotations")


def index_of_value(vector, value):
    """
    Return the 0-based index of a numpy array which
    contains the specific value
    """
    try:
        return np.where(vector == value)[0][0]
    except:
        return


def script_init():
    """
    Start initializations
    """
    session_id = random_hexadecimal_string()

    prefix = extract_file_name(global_config["bam_file_path"])

    suffix = "_{}_{}_offsets_{}".format(
        prefix,
        global_config["read_offset_terminal"],
        str(dt.now().strftime("%I-%M-%S%p")),
    )

    # Creating a directory, WITHOUT checking for pre-existing ones
    # with the same name. This is exactly what we want in case of
    # re runs, and we are including the timestamp in the name.
    basename = "{}_Results_{}_{}_offsets_on_{}".format(
        session_id,  # prevents any screwy conflicts
        prefix,
        global_config["read_offset_terminal"],
        str(dt.now().strftime("%d-%b-%Y_at_%I-%M-%S%p")),
    )

    # Filename conflicts. We should never see this.
    if os.path.isdir(basename):
        notify("Are you running this script in parallel?", level="crit")
        notify(
            "Please ensure the BAM files are named differently.",
            level="crit",
            fatal=True,
        )

    img_dir = "{}/{}".format(basename, "Plots")
    data_dir = "{}/{}".format(basename, "Data")

    os.makedirs(img_dir)
    os.makedirs(data_dir)

    # Set the progression type for p-site offsets
    if global_config["read_offset_terminal"] == "five_prime":
        progression = "increasing"

    if global_config["read_offset_terminal"] == "three_prime":
        progression = "decreasing"

    # Set up the file paths
    set_global_config(
        offset_progression=progression,
        # Keep track of how many times the metagene
        # function has been called. This is because
        # the first call will be without offsets, and
        # we need to keep track of it so that we can name
        # the output files accordingly.
        metagene_callcount=0,
        suffix=suffix,
        img_dir=img_dir + "/",  # this breaks Windows compatibility
        data_dir=data_dir + "/",
        img_id="_{}_{}".format(prefix, str(dt.now().strftime("%I-%M-%S%p"))),
        log_file="{}/{}_logfile_for_{}_{}.txt".format(
            basename, session_id, prefix, str(dt.now().strftime("%I-%M-%S%p"))
        ),
        conf_file="{}/{}_runconfig_{}_{}.json".format(
            basename, session_id, prefix, str(dt.now().strftime("%I-%M-%S%p"))
        ),
        session_id=session_id,
        prefix=prefix,
        script_version=versionstr,
    )

    notify("Started {}".format(global_config["script_version"]))
    notify("Session ID is: {}".format(session_id))

    try:
        global_config["bamfile"] = pysam.AlignmentFile(
            global_config["bam_file_path"], "rb"
        )

        total_mapped_reads = int(
            pysam.flagstat(global_config["bam_file_path"])
            .split("\n")[4]
            .split("+")[0]
            .strip(" ")
        )

        if total_mapped_reads < 1e6:
            notify(
                "{} mapped reads (< 1e6). Results may be unreliable".format(
                    total_mapped_reads
                ),
                level="crit",
            )

        global_config["mapped_reads"] = total_mapped_reads

        notify("Found BAM file ({}), seems OK".format(prefix), level="notf")
    except:
        notify("BAM file is unavailable.", level="crit", fatal=True)

    # Check if BAM index file exists. This can be fatal.
    # We should try to index it ourselves.
    if not os.path.isfile(
        os.path.split(global_config["bam_file_path"])[0]
        + "/"
        + extract_file_name(global_config["bam_file_path"])
        + ".bam.bai"
    ):
        notify(
            "BAM file index does not exist. Please create one.",
            level="crit",
            fatal=True,
        )
    else:
        notify("Found BAM index (did not check contents)")

    if ("coverage_cutoff" not in global_config) or (
        "coverage_metric" not in global_config
    ):

        notify("Will filter genes with coverage < 10 RPKM")
        global_config["coverage_cutoff"] = 10
        global_config["coverage_metric"] = "rpkm"

    if global_config["peak_scan_range"][0] > global_config["peak_scan_range"][1]:
        # swap
        global_config["peak_scan_range"][0],
        global_config["peak_scan_range"][1] = (global_config["peak_scan_range"][1],)
        global_config["peak_scan_range"][0]

        notify(
            "peak_scan_range: should be set as [lowerbound, upperbound]", level="warn"
        )
        notify(
            "Example: [-20, -10] or [0, 30] etc, where 0 denotes start", level="warn"
        )

    # For a gene, get read coverage from 5' or 3' end
    if "readcov_from_terminal" not in global_config:
        global_config["readcov_from_terminal"] = "five_prime"

    # Calculate coverage for N% of gene length
    if "readcov_percentage_length" not in global_config:
        global_config["readcov_percentage_length"] = 100

    # #nucleotides to ignore from the 5' or 3' end
    # when calculating the coverages
    if "readcov_ignore_nts" not in global_config:
        global_config["readcov_ignore_nts"] = 20

    # #nucleotides upstream of start codon to include in metagene
    if "endmap_upstream" not in global_config:
        global_config["endmap_upstream"] = 50

    # #nucleotides downstream of start codon to include in metagene
    if "endmap_downstream" not in global_config:
        global_config["endmap_downstream"] = 200

    # We calculate the read coverage for a few additional nucleotides
    # upstream of endmap_upstream and downstream of endmap_downstream,
    # so that we dont loose any information when rotating the metagene
    # vectors by their offsets
    if "mapping_buffer" not in global_config:
        global_config["mapping_buffer"] = 25

    # By default we use heuristics to determine the best RPF
    # readlengths to use for the analysis.
    # The user may supply an override if necessary.

    if len(global_config["use_readlengths"]) == 0:
        global_config["readlengths"] = [i for i in range(20, 41, 1)]
    else:
        global_config["readlengths"] = global_config["use_readlengths"]

    # Check if we need to include any overlapping genes
    try:
        overlap_exception = read_genes_from_file(
            global_config["genes_overlap_exception"]
        )
        set_global_config(overlap_exception=overlap_exception)
    except:
        set_global_config(overlap_exception=[])

    parse_gff()
    detect_overlaps()


def check_gene_list():
    """
    Check if a gene list has been provided. If so,
    should we include or exclude the genes from the analysis?
    """
    try:
        listfile = global_config["gene_list_file"]
        listaction = global_config["gene_list_action"]
    except:
        notify(
            "Gene list file or associated action not mentioned.\
                     Skipping.",
            level="warn",
        )
        return

    def remove_annotation(gene):
        try:
            del annotation[gene]
        except:
            notify("{} not found in annotation".format(gene), level="warn")
            notify("Ignoring {}".format(gene), level="warn")

    # Should we perform the analysis using only
    # these genes? Or completely excluding these
    # genes?
    total_annotation_size = len(annotation)
    allgenes = annotation.copy().keys()  # will throw error otherwise
    actions = ["include_only", "exclude_only", "exclude_balance"]
    genelist = []

    try:
        genelist = read_genes_from_file(listfile)
    except:
        notify("No gene list provided. Moving on.", level="warn")
        return

    if listaction not in actions:
        notify("What to do with the gene list?", level="warn")
        notify("Try: 'include_only', 'exclude_only' or 'exclude_balance'", level="warn")
        notify("Ignoring provided gene list", level="warn")
        return

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
        balance_sample = random.sample(reduced_genes, len(genelist))

        for gene in allgenes:
            if gene not in balance_sample:
                remove_annotation(gene)

    if listaction == "include_only":
        # Remove all genes except these genes
        # from the annotation
        for gene in allgenes:
            if gene not in genelist:
                remove_annotation(gene)

    notify(
        "Internal annotation changed from {} to {} genes".format(
            total_annotation_size, len(annotation)
        )
    )


def save_file(data, filename):
    """
    Dump an object to disk as JSON
    """
    try:
        save_json(
            data,
            "{}{}_{}{}.json".format(
                global_config["data_dir"],
                global_config["session_id"],
                filename,
                global_config["suffix"],
            ),
        )
    except:
        notify(
            "JSON marshal failed for object: {} (will ignore)".format(filename),
            level="warn",
        )


def generate_coverage_profiles(include_genes=[], disabled=False):
    """
    Extract the read coverage profile for all genes
    """
    if disabled:
        return

    for gene, details in annotation.items():

        if len(include_genes) != 0:
            if gene not in include_genes:
                continue

        coverage_profile[gene] = get_coverage(
            start=details["start"],
            stop=details["stop"],
            reference=details["refname"],
            strand=details["strand"],
            mode="profile",
        ).copy()

    # dump the coverage information
    save_file(coverage_profile, "read_coverage_profile")


def get_coverage(start, stop, reference, strand, mode="mean"):
    """
    Get the read coverage for a gene
    """

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
        cv = bamfile.count_coverage(
            reference, start=start, stop=stop, read_callback=watson_only
        )
    if strand == "-":
        cv = bamfile.count_coverage(
            reference, start=start, stop=stop, read_callback=crick_only
        )

    if mode == "mean":
        try:
            mean_cov = np.sum(np.sum(cv, axis=0)) / len(cv)
        except ZeroDivisionError:
            mean_cov = 0
        return mean_cov

    if mode == "profile":
        cv = np.sum(cv, axis=0)
        return {"coverage": cv.tolist(), "position": np.arange(start, stop, 1).tolist()}


def plot_metagene_per_readlength(plot_title="", figsize_x=10, figsize_y=8, labsize=14):
    """
    Plot the metagene vectors for each RPF read length
    (with and without offsets)
    """

    up = global_config["endmap_upstream"]
    dn = global_config["endmap_downstream"]
    nt_range = global_config["peak_scan_range"]
    data = metagene_per_readlength

    savedata = {}

    # MAX_PLOTS = len(data)
    MAX_PLOTS = len(offsets_p) + 1
    xaxis = str(list(data.keys())[0])  # x-axis identifier
    x_vector = data[xaxis][: (up + dn)]  # trim it down
    savedata["xaxis"] = x_vector.tolist()

    fig, ax = plt.subplots(
        nrows=MAX_PLOTS - 1,
        ncols=1,
        sharex=True,
        sharey=True,
        figsize=(figsize_x, figsize_y),
    )

    plotcount = 1  # for plot
    axiscount = 0  # for axis

    # TODO: why do the same thing twice?
    for rpf_len, metagene_vector in data.items():
        if rpf_len not in list(offsets_p.keys()):
            continue
        # This is a brain-dead up attempt to fix a bug. Don't read
        # too much into this.
        if len(metagene_vector) == 1 and metagene_vector[0] == -1:
            notify(
                "Can't compute metagene for RPF length {} (no data)".format(rpf_len),
                level="warn",
            )
            continue

        savedata[rpf_len] = []
        savedata[rpf_len].append([metagene_vector.tolist()])

        if rpf_len == xaxis:
            continue
        plotcount += 1
        if plotcount == MAX_PLOTS + 1:
            break

        palette = next(colors)["color"]
        ax[axiscount].plot(
            x_vector, metagene_vector, label="{} mer".format(rpf_len), color=palette
        )
        ax[axiscount].xaxis.grid(which="major", color="k", linestyle="-", linewidth=0.2)
        ax[axiscount].legend(loc="upper right")

        # In case we are given a range of values within which to
        # idenify and annotate the p-site offsets
        if len(nt_range) == 2:
            peakpos, offset, _ = get_offset(
                xvec=x_vector,
                yvec=metagene_vector,
                left_bound=nt_range[0],
                right_bound=nt_range[1],
            )

            if global_config["metagene_callcount"] == 0:

                ax[axiscount].axvline(x=peakpos, color="black", linestyle="--")
                savedata[rpf_len].append([peakpos])

                ax[axiscount].text(
                    peakpos - 12, 1, offset, rotation=0, color="black", weight="bold"
                )
        ax[axiscount].axvline(x=0, color="red", linestyle="--")
        axiscount += 1

    fig.suptitle(plot_title)
    fig.text(0.5, 0.04, "Nt from CDS start", ha="center", fontsize=labsize)
    fig.text(0.04, 0.5, "A.U.", va="center", rotation="vertical", fontsize=labsize)

    fname_line = ""  # Ugly. I know.
    fname_bar = ""

    if global_config["metagene_callcount"] == 0:
        fname_line = "metagene_kmer_no_offset"
        fname_bar = "metagene_kmer_framing_no_offset"
    else:
        fname_line = "metagene_kmer_offset"
        fname_bar = "metagene_kmer_framing_offset"

    save_file(savedata, fname_line)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            fname_line,
            global_config["img_id"],
        )
    )

    plt.close()

    # Phasing plots
    fig, ax = plt.subplots(
        nrows=MAX_PLOTS - 1,
        ncols=1,
        sharex=True,
        sharey=True,
        figsize=(figsize_x, figsize_y),
    )

    plotcount = 1  # for plot
    axiscount = 0  # for axis

    for rpf_len, metagene_vector in data.items():

        if rpf_len not in list(offsets_p.keys()):
            continue

        if len(metagene_vector) == 1 and metagene_vector[0] == -1:
            notify(
                "Can't compute metagene for RPF length {} (no data)".format(rpf_len),
                level="warn",
            )
            continue

        if rpf_len == xaxis:
            continue
        plotcount += 1
        if plotcount == MAX_PLOTS + 1:
            break

        bars = ax[axiscount].bar(
            x_vector[50:100],
            metagene_vector[50:100],
            label="{} mer".format(rpf_len),
            color="#11a0ff",
        )

        for item in bars[::3]:
            item.set_color("#ff5842")

        ax[axiscount].xaxis.grid(which="major", color="k", linestyle="-", linewidth=0.2)
        ax[axiscount].legend(loc="upper right")

        axiscount += 1

    fig.suptitle(plot_title)
    fig.text(0.5, 0.04, "Nt from CDS start", ha="center", fontsize=labsize)
    fig.text(0.04, 0.5, "A.U.", va="center", rotation="vertical", fontsize=labsize)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            fname_bar,
            global_config["img_id"],
        )
    )

    plt.close()

    global_config["metagene_callcount"] += 1


def get_gene_rpkm(gene_annotation):
    """
    Returns the read coverage for a gene in RPKM units
    """

    strand = gene_annotation["strand"]
    read_lengths = global_config["readlengths"]
    genomic_start = gene_annotation["start"]
    genomic_stop = gene_annotation["stop"]
    gene_length = genomic_stop - genomic_start
    reference_id = gene_annotation["refname"]
    bamfile = global_config["bamfile"]
    mapped_reads_library = global_config["mapped_reads"]

    mapped_reads_gene = 0

    for read in bamfile.fetch(
        reference=reference_id, start=genomic_start, stop=genomic_stop
    ):

        if read.reference_length not in read_lengths:
            continue
        if (strand == "+") and (read.is_reverse):
            continue
        if (strand == "-") and (not read.is_reverse):
            continue

        mapped_reads_gene += 1

    scaling_factor = mapped_reads_library / 1e6  # per million mapped reads
    rpm = (
        mapped_reads_gene / scaling_factor
    )  # normalize counts for read depth, reads per million
    rpkm = rpm / ((gene_length) / 1e3)  # rpkm = rpm/gene length in kB
    return rpkm


def gene_coverage(gene_annotation):
    """
    Wrapper function around get_coverage(). Modifies the
    gene cordinates to account for padding, %length,
    5' or 3' coverages, etc
    """

    strand = gene_annotation["strand"]
    genomic_start = gene_annotation["start"]
    genomic_stop = gene_annotation["stop"]
    reference_id = gene_annotation["refname"]
    gene_length = genomic_stop - genomic_start
    ignore_terminal_nt = global_config["readcov_ignore_nts"]
    from_which_end = global_config["readcov_from_terminal"]
    percentage_length = global_config["readcov_percentage_length"]
    length_to_consider = (percentage_length / 100) * gene_length

    # Is the gene long enough for coverage calculation?
    if (gene_length >= (ignore_terminal_nt + 2)) == False:
        return 0

    # Determine the genomic coordinate range within which to
    # calculate the coverage.
    if strand == "+":
        if from_which_end == "five_prime":
            start = genomic_start + ignore_terminal_nt
            stop = (
                genomic_start + percentage_length
            )  # w00t! using percentage directly??
        if from_which_end == "three_prime":
            start = genomic_stop - percentage_length
            stop = genomic_stop

    if strand == "-":
        if from_which_end == "five_prime":
            start = genomic_stop - percentage_length
            stop = genomic_stop - ignore_terminal_nt
        if from_which_end == "three_prime":
            start = genomic_start
            stop = genomic_start + percentage_length

    # Get the mean read coverage for a specific length of the genome
    mc = get_coverage(
        start=start, stop=stop, strand=strand, reference=reference_id, mode="mean"
    )
    return mc


def calculate_gene_coverages(mode="reads_per_nt"):
    """
    Wrapper function to go through the annotation and
    calculate the read coverage for each gene. It then
    sorts the genes based on the coverage value.
    """

    tmp = {}  # temporarily store coverages before sorting
    i, gs = 0, len(annotation.keys())
    if mode == "reads_per_nt":  # average reads per nucleotide
        for gene_name, details in annotation.items():
            tmp[gene_name] = gene_coverage(details)
            i += 1
            print(" >> {:.2f}%".format(round(i / gs * 100, 2)), end="\r")

    if mode == "rpkm":
        for gene_name, details in annotation.items():
            tmp[gene_name] = get_gene_rpkm(details)
            i += 1
            print(" >> {:.2f}%".format(round(i / gs * 100, 2)), end="\r")

    tmp_tuple = sorted(tmp.items(), key=lambda item: item[1], reverse=True)

    for index in range(len(tmp_tuple)):
        gene_coverages[tmp_tuple[index][0]] = tmp_tuple[index][1]

    save_file(gene_coverages, "gene_coverages")


def plot_initiation_peak(peak=False, peak_range=[]):
    """
    Plot the initiation peak
    """
    up = global_config["endmap_upstream"]
    dn = global_config["endmap_downstream"]
    xaxis = metagene_per_readlength["xaxis"][: (up + dn)]

    combined_vectors = []

    for length, vector in metagene_per_readlength.items():

        if length not in list(offsets_p.keys()):
            continue

        if len(vector) == 1 and vector[0] == -1:
            notify(
                "Not plotting metagene for RPF length {}".format(length), level="warn"
            )
            continue

        if length == "xaxis":
            continue
        combined_vectors.append(vector.tolist()[: (up + dn)])

    metagene_vector = np.sum(combined_vectors, axis=0)

    plt.plot(xaxis, metagene_vector, color="black")

    if peak and len(peak_range) != 0:

        peakpos, delta, _ = get_offset(
            xvec=xaxis, yvec=metagene_vector, right_bound=180, left_bound=2
        )
        plt.axvline(x=delta, color="r")
        plt.text(
            peakpos + 0.4, -2, delta, rotation=0, color="black", weight="bold"
        )  # 0.4 to keep the text away from the line

    start_data = {  # include the plot configurations as well?
        "peakpos": peakpos,
        "delta": delta,
        "xaxis": xaxis.tolist(),
        "start_metagene": metagene_vector.tolist(),
    }

    save_file(start_data, "start_position_peak")

    plt.xlabel("Nt from CDS start", fontsize=18)
    plt.ylabel("A.U.", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.tight_layout()
    plt.gcf().set_size_inches(19, 5)
    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "start_peak",
            global_config["img_id"],
        )
    )

    plt.close()

    bars = plt.bar(xaxis[50:170], metagene_vector[50:170], color="#11a0ff")

    for item in bars[::3]:
        item.set_color("#ff5842")

    plt.xlabel("Nt from CDS start", fontsize=18)
    plt.ylabel("A.U.", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.tight_layout()
    plt.gcf().set_size_inches(19, 5)
    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "start_peak_framing",
            global_config["img_id"],
        )
    )

    plt.close()


def get_max_coverage_genes():
    """
    Identify genes with read coverage above a certain
    cutoff.
    """
    coverage_metric = global_config["coverage_metric"]
    coverage_cutoff = global_config["coverage_cutoff"]

    if len(gene_coverages) == 0:
        calculate_gene_coverages(mode=coverage_metric)

    for gene_name, coverage in gene_coverages.items():
        if coverage >= coverage_cutoff:
            high_coverage_genes.append(gene_name)

    save_file(high_coverage_genes, "high_coverage_genes")


def skip_read(
    read,  # the read sequence
    skip_dict,  # codons, if found in E/P/A site would be skipped
    offset,  # The p site offset for the read
    strand,
):
    """
    Check whether a read has a particular codon in it's
    E/P/A site, and if so, whether it should be skipped
    during the analysis.
    """
    skip = True
    E, P, A = "", "", ""  # the E/P/A site codons of a read

    if strand == "+":  # read maps to watson
        E = read[offset - 3 : offset]
        P = read[offset : offset + 3]
        A = read[offset + 3 : offset + 6]

    if strand == "-":  # read maps to crick
        # reverse the codon
        offset = len(read) - offset
        A = read[offset - 6 : offset - 3][::-1]
        P = read[offset - 3 : offset][::-1]
        E = read[offset : offset + 3][::-1]

    if (E in skip_dict["E"]) or (P in skip_dict["P"]) or (A in skip_dict["A"]):
        skip = False

    return skip


# INPUT:  1. Gene annotation dictionary for a specific gene
#         2. The PRF terminal we want to map (5' or 3')
# OUTPUT: 1. Defaultdict, mapping the length of an RPF to a
#            list of it's 5' or 3' end alignment coordinates
#            eg: master_dict[28] = [1000, 1000, 1002,1050,1050,1050, ...]
#         2. Small config dictionary containing the start and
#            stop coordinates after padding (upstream+buffer)
#            along with the strand information
def terminal_alignment_positions_per_readlength(
    gene_annotation,
    apply_offset=False,
    offset_site="P",
    ignore_reads_by_nt={5: [], 3: []},  # This option is disabled
    codons_to_skip={"E": [], "P": [], "A": []},
):
    """
    For each read length, get a list of alignment positions of
    their 5' or 3' ends (for a specific gene)
    """
    master_dict = defaultdict(list)
    config = {}
    genomic_start = gene_annotation["start"]
    genomic_stop = gene_annotation["stop"]
    gene_length = abs(genomic_start - genomic_stop)
    strandinfo = gene_annotation["strand"]
    reference_id = gene_annotation["refname"]
    upto_upstream_nt = global_config["endmap_upstream"]
    upto_downstrm_nt = global_config["endmap_downstream"]
    bamfile = global_config["bamfile"]
    _buffer_ = global_config["mapping_buffer"]
    terminal_map = global_config["read_offset_terminal"]
    read_lengths = global_config["readlengths"]

    if offset_site == "P":
        offsets = offsets_p
    if offset_site == "A":
        offsets = offsets_a
    if offset_site == "E":
        offsets = offsets_e

    if apply_offset and (len(offsets) == 0):
        notify("Tried to apply offsets but none available", level="warn", onetime=True)
        apply_offset = False

    if strandinfo == "+":
        # watson strand genes
        start = genomic_start - (upto_upstream_nt + _buffer_)
        stop = genomic_start + (upto_downstrm_nt + _buffer_) + gene_length

    if strandinfo == "-":
        # crick strand genes
        start = genomic_stop - (upto_downstrm_nt + _buffer_) - gene_length
        stop = genomic_stop + (upto_upstream_nt + _buffer_)

    # Will be returned with the maps
    config["start"] = start
    config["stop"] = stop

    # Ignore genes too close to chromosomal start (including buffering)
    if start <= 0:
        return (master_dict, config)

    i, j = 0, 0

    for read in bamfile.fetch(reference=reference_id, start=start, stop=stop):

        if read.reference_length not in read_lengths:
            continue

        if (strandinfo == "+") and (read.is_reverse):
            continue

        if (strandinfo == "-") and (not read.is_reverse):
            continue

        # if read.is_reverse:
        #     watson_read = False

        # Does not take into account the soft/hard clipping of the the read sequences
        # Probably should compare the cigars?
        # read_sequence = read.query_sequence

        # if watson_read:
        #     base5  = read_sequence[0]
        #     base3  = read_sequence[len(read_sequence)-1]
        # else:
        #     base5  = read_sequence[len(read_sequence)-1]
        #     base3  = read_sequence[0]

        # # Any reads with set nucleotides from the 5' or
        # # 3' ends will be skipped.
        # if (len(ignore_read_by_nt[5]) != 0):
        #     if base5 in ignore_read_by_nt[5]:
        #        # skip the read
        #        continue
        # if (len(ignore_read_by_nt[3]) != 0):
        #     if base3 in ignore_read_by_nt[3]:
        #        # skip the read
        #        continue

        if (
            (len(codons_to_skip["E"]) != 0)
            or (len(codons_to_skip["P"]) != 0)
            or (len(codons_to_skip["A"]) != 0)
        ):

            seq = read.query_sequence
            delta = abs(
                value_from_dict(
                    offsets_p,  # we only consider the p-site offset
                    read.reference_length,
                )
            )

            if read.is_reverse:
                if skip_read(seq, codons_to_skip, delta, "-"):
                    i += 1
                    continue

            if not read.is_reverse:
                if skip_read(seq, codons_to_skip, delta, "+"):
                    j += 1
                    continue

        terminal_coordinate = 0

        # 5' coordinate of watson strand read
        if terminal_map == "five_prime" and not read.is_reverse:
            terminal_coordinate = read.reference_start + 1
        # 5' coordinate of crick strand read
        if terminal_map == "five_prime" and read.is_reverse:
            terminal_coordinate = read.reference_end
        # 3' coordinate of watson strand read
        if terminal_map == "three_prime" and not read.is_reverse:
            terminal_coordinate = read.reference_end
        # 5' coordinate of crick strand read
        if terminal_map == "three_prime" and read.is_reverse:
            terminal_coordinate = read.reference_start + 1
        # TODO: why are we using np.array() here again?
        if apply_offset == False:
            master_dict[read.reference_length].append(np.array(terminal_coordinate))

        if apply_offset == True and terminal_map == "five_prime":

            if not read.is_reverse:
                master_dict[read.reference_length].append(
                    np.array(terminal_coordinate)
                    + abs(value_from_dict(offsets, read.reference_length))
                )
            if read.is_reverse:
                master_dict[read.reference_length].append(
                    np.array(terminal_coordinate)
                    - abs(value_from_dict(offsets, read.reference_length))
                )

        if apply_offset == True and terminal_map == "three_prime":

            if not read.is_reverse:
                master_dict[read.reference_length].append(
                    np.array(terminal_coordinate)
                    - abs(value_from_dict(offsets, read.reference_length))
                )
            if read.is_reverse:
                master_dict[read.reference_length].append(
                    np.array(terminal_coordinate)
                    + abs(value_from_dict(offsets, read.reference_length))
                )

    return (master_dict, config)


# INPUT:  1. Output from terminal_alignment_positions_per_readlength()
#         2. The annotation data of the gene in question
#         3. The buffered start and stop coordinates over which
#            to calculate the vectors
#
# OUTPUT: Dictionary mapping gene_name --> readlength --> endmap_vector
def position_list_to_endmap_vector(mapdict, config, gene_annotation):
    """
    Returns a vector, each element denoting the number of 5' or 3'
    read terminals aligning per nucleotide position (represented
    by the element indices). Separate vectors for each read length.
    For a single gene.
    """
    start = config["start"]
    stop = config["stop"]
    # Gene annotation details
    strand = gene_annotation["strand"]
    genomic_start = gene_annotation["start"]
    genomic_stop = gene_annotation["stop"]
    gene_length = abs(genomic_start - genomic_stop)
    # How many nt up and downstream of the start codon are
    # we supposed to consider
    upto_upstream_nt = global_config["endmap_upstream"]
    upto_downstrm_nt = global_config["endmap_downstream"]
    # master
    tmp_maps = {}
    # Convert the defaultdict list to numpy arrays for
    # easier handling
    for length, mappings in mapdict.items():
        if length not in global_config["readlengths"]:
            continue
        tmp_maps[length] = np.array(mappings)

    # Reset genomic mappings so that start codon gets index 0
    if strand == "+":
        for length, mappings in tmp_maps.items():
            tmp_maps[length] = mappings - genomic_start
    if strand == "-":
        for length, mappings in tmp_maps.items():
            tmp_maps[length] = genomic_stop - mappings

    # A numpy array holding the indexes around the genomic region
    indexes = np.array(
        [i for i in range(-upto_upstream_nt, upto_downstrm_nt + gene_length + 1, 1)]
    )

    # Indexes for plotting the metagene
    metagene_per_readlength["xaxis"] = indexes
    # Zeroed vector, each element will be incremented by the read counts
    # Element indices equivalent to nucleotide positons
    template = np.zeros(len(indexes))

    # Convert endmaps to vectors
    for length, mappings in tmp_maps.items():
        nullvec = template.copy()

        for maps in mappings:
            if (maps >= (-upto_upstream_nt)) and (
                maps <= upto_downstrm_nt + gene_length
            ):
                index = index_of_value(indexes, maps)
                if index == None:
                    # Element doesn't exist
                    continue

                nullvec[index] = nullvec[index] + 1

        tmp_maps[length] = nullvec.copy()  # not ideal
        nullvec = np.delete(nullvec, [i for i in range(len(nullvec))])

    return tmp_maps


def detect_overlaps():
    """
    Identify overlapping genes. This is quite rudimentary.
    """
    overlap_fiveprime = []  # genes overlapping at their 5' ends
    overlap_threeprime = []  # genes overlapping at their 3' ends

    # Convert the dictionary into a list of lists
    annotation_list = []
    for gene_name, details in annotation.items():
        temp = [gene_name, details]
        annotation_list.append(temp)

    for i in range(len(annotation_list)):

        a = i
        b = a + 1

        if b == len(annotation_list):
            break

        gene_upstream = annotation_list[a]
        gene_downstream = annotation_list[b]  # <-- does this overlap?
        # A gene is overlapping from the five prime end
        # if it's start coordinate is less than the
        # stop coordinate of the downstream gene.
        stop_gene_upstream = gene_upstream[1]["stop"]
        start_gene_downstream = gene_downstream[1]["start"]

        # Genes need to be on the same strand to be considered
        upstream_gene_strand = gene_upstream[1]["strand"]
        downstream_gene_strand = gene_downstream[1]["strand"]

        if (start_gene_downstream < stop_gene_upstream) and (
            upstream_gene_strand == downstream_gene_strand
        ):
            # The gene is overlapping
            # We record the name of the gene
            overlap_fiveprime.append(gene_downstream[0])
            overlap_threeprime.append(gene_upstream[0])

        a += 1
        b += 1

    if (len(overlap_fiveprime) != 0) or (len(overlap_threeprime) != 0):
        notify("There are overlapping genes", level="warn", onetime=True)

    # genes overlapping at their 5' end
    overlap_dict["five_prime"] = overlap_fiveprime
    # genes overlapping at their 3' end
    overlap_dict["three_prime"] = overlap_threeprime


def map_gene_to_endmaps(
    apply_offset=False, only_reads_with_codons={"E": [], "P": [], "A": []}
):
    """
    Loops over each gene and calculates their endmap vectors.
    This is the driver code.
    """
    ignore_overlap = global_config["will_ignore_overlaps"]

    try:
        overlap_exception = global_config["overlap_exception"]
    except:
        overlap_exception = []

    # Diagnostic purposes only.
    track_overlaps = []
    track_coverages = []
    overall_skipped = []
    genes_processed = []
    genes_used = []

    notify(
        "{} genes in high coverage list, out of {}".format(
            len(high_coverage_genes), len(annotation.keys())
        ),
        level="notf",
    )

    notify("Starting endmap scan (this can take some time)")

    i, ant_sz = 0, len(high_coverage_genes)

    for gene_name, details in annotation.items():

        print(
            " >> {:.2f}% | {}".format(round(i / ant_sz * 100, 2), gene_name), end="\r"
        )

        genes_processed.append(gene_name)

        if gene_name not in high_coverage_genes:
            track_coverages.append(gene_name)
            overall_skipped.append(gene_name)
            continue

        if (ignore_overlap) and (gene_name not in overlap_exception):
            # The "start" of a gene on the crick strand
            # is actually it's genomic stop
            if (details["strand"] == "-") and (
                gene_name in overlap_dict["three_prime"]
            ):
                track_overlaps.append(gene_name)
                overall_skipped.append(gene_name)
                continue

            if (details["strand"] == "+") and (gene_name in overlap_dict["five_prime"]):
                track_overlaps.append(gene_name)
                overall_skipped.append(gene_name)
                continue

        genes_used.append(gene_name)

        # Get read terminal alignment coordinates
        maps_P, cnfs_P = terminal_alignment_positions_per_readlength(
            details,
            apply_offset,
            offset_site="P",
            codons_to_skip=only_reads_with_codons,
        )

        maps_E, cnfs_E = terminal_alignment_positions_per_readlength(
            details,
            apply_offset,
            offset_site="E",
            codons_to_skip=only_reads_with_codons,
        )

        maps_A, cnfs_A = terminal_alignment_positions_per_readlength(
            details,
            apply_offset,
            offset_site="A",
            codons_to_skip=only_reads_with_codons,
        )

        # Convert coordinate list to endmap vector
        endmaps_P = position_list_to_endmap_vector(maps_P, cnfs_P, details)
        endmaps_E = position_list_to_endmap_vector(maps_E, cnfs_E, details)
        endmaps_A = position_list_to_endmap_vector(maps_A, cnfs_A, details)

        endmap_vectors["E"][gene_name] = endmaps_E
        endmap_vectors["P"][gene_name] = endmaps_P
        endmap_vectors["A"][gene_name] = endmaps_A

        i += 1

        print("                                      ", end="\r")

    notify(
        "{} gene(s) skipped due to overlaps".format(len(track_overlaps)), level="warn"
    )
    notify(
        "{} gene(s) skippped due to low coverage".format(len(track_coverages)),
        level="warn",
    )
    notify("{} gene(s) processed in total".format(len(genes_processed)), level="notf")
    notify(
        "{} gene(s) skipped out of {} in annotation".format(
            len(overall_skipped), len(annotation.keys())
        ),
        level="warn",
    )

    # print out the gene names which were in the high coverage list but skipped
    save_file(unique(track_overlaps), "skipped_genes_overlaps_endmaps.json")
    save_file(unique(track_coverages), "skipped_genes_low_coverage_endmaps.json")
    save_file(unique(genes_used), "genes_used_endmaps.json")


def normalize_counts(vector):
    """
    Normalize a numpy array by it's sum
    """
    if sum(vector) == 0:
        return vector
    else:
        return vector / sum(vector)


def frame_for_coordinate(genomic_start, genomic_stop, query_index, strand):
    """
    Returns the reading frame for an arbitrary genomic coordinate
    """
    if (query_index > genomic_stop) or (query_index < genomic_start):
        return 0

    if strand == "+":
        delta = query_index - genomic_start
        return (delta % 3) + 1
    if strand == "-":
        delta = genomic_stop - query_index
        return (delta % 3) + 1


def reading_frame(ignore_read_by_nt={5: [], 3: []}):
    """
    Fraction of reads (p-sites) mapping to each reading frame.
    """
    reads_in_frame = {}
    kmer_include = global_config["readlengths"]

    offset_direction = global_config["read_offset_terminal"]

    offset_dict = offsets_p

    reads_in_frame["total"] = 0
    reads_in_frame[1] = 0
    reads_in_frame[2] = 0
    reads_in_frame[3] = 0

    total, processed, skipped, skipped_other = 0, 0, 0, 0

    for gene_name, details in annotation.items():

        genomic_start = details["start"]  # 1 indexed
        genomic_stop = details["stop"]
        strand = details["strand"]
        refname = details["refname"]
        bamfile = global_config["bamfile"]

        # Ignore genes which are too short or have low coverage
        if (abs(genomic_start - genomic_stop) < 100) or (
            gene_name not in high_coverage_genes
        ):
            continue

        # 27 and 12 from DOI: 10.7554/eLife.42591.001
        if strand == "+":
            start = genomic_start + 27
            stop = genomic_stop - 12

        if strand == "-":
            start = genomic_start + 12
            stop = genomic_stop - 27

        for read in bamfile.fetch(reference=refname, start=start, stop=stop):

            total += 1

            if (strand == "+") and (read.is_reverse):
                skipped_other += 1
                continue
            if (strand == "-") and (not read.is_reverse):
                skipped_other += 1
                continue

            # Skip reads with a specific base at their 5' or 3' terminals
            read_sequence = read.query_sequence

            if not read.is_reverse:
                base5 = read_sequence[0]
                base3 = read_sequence[len(read_sequence) - 1]
            else:
                base5 = read_sequence[len(read_sequence) - 1]
                base3 = read_sequence[0]

            # Any reads with set nucleotides from the 5' or
            # 3' ends will be skipped.
            if len(ignore_read_by_nt[5]) != 0:
                if base5 in ignore_read_by_nt[5]:
                    skipped += 1
                    continue

            if len(ignore_read_by_nt[3]) != 0:
                if base3 in ignore_read_by_nt[3]:
                    skipped += 1
                    continue

            # Only calculate frame info for a particular length
            if (len(kmer_include) != 0) and (read.reference_length not in kmer_include):
                skipped_other += 1
                continue

            read_map = 0
            offset = 0

            if len(offset_dict) != 0:
                try:
                    offset = offset_dict[read.reference_length]
                except:
                    notify(
                        "Fragment length {} has no offset".format(
                            read.reference_length
                        ),
                        level="warn",
                        onetime=True,
                    )
                    # Ignore
                    skipped_other += 1
                    continue

            if (strand == "+") and (offset_direction == "five_prime"):
                read_map = (read.reference_start + 1) + abs(offset)

            if (strand == "+") and (offset_direction == "three_prime"):
                read_map = read.reference_end - abs(offset)

            if (strand == "-") and (offset_direction == "five_prime"):
                read_map = read.reference_end - abs(offset)

            if (strand == "-") and (offset_direction == "three_prime"):
                read_map = (read.reference_start + 1) + abs(offset)

            read_frame = frame_for_coordinate(
                genomic_start, genomic_stop, read_map, strand
            )
            processed += 1

            try:
                reads_in_frame[read_frame] += 1
                reads_in_frame["total"] += 1
            except:
                skipped_other += 1
                pass

    save_file(reads_in_frame, "global_reading_frame")

    return reads_in_frame


def metagene_over_codon(codon_list, mode="combined"):
    """
    Metagene over a particular codon. For all genes.
    Drive code.
    """
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


def chunk_from_id(vec, idx, length):
    """
    Return a portion of a vector +/-N
    elements from an index
    """
    if (idx + length + 1) > len(vec):
        raise RuntimeError
    if (idx - length) < 0:
        raise RuntimeError
    return vec[idx - length : idx + length + 1]


def plot_metagene_over_codon(codon_list, site="P"):
    """
    Plots out the metagene over a codon. Defaults
    to P-site offsets.
    """
    codon_positions = {}
    cmetagene = []  # metagene over a specific codon

    if len(transcripts_dict) == 0:
        retrieve_coding_sequence(
            global_config["coding_sequence_path"],
            global_config["coding_sequence_format"],
        )

    i = 0
    x = ["-2", "E", "P", "A", "2"]

    for gene_name, details in annotation.items():

        if gene_name not in high_coverage_genes:
            continue

        try:
            pps = pps_vector_per_gene[site][gene_name]["pps"]
        except:
            notify(
                "{} does not have a positional pause vector".format(gene_name),
                level="warn",
            )
            continue

        try:
            seq = pps_vector_per_gene[site][gene_name]["seq"]
        except:
            notify(
                "{} does not have a transcript sequence".format(gene_name),
                level="warn",
            )
            continue

        if i == 401:
            break

        i += 1

        tmp = {}
        # tmp["seq"] = seq
        tmp["pps"] = pps
        tmp["pos"] = position_of_codon("".join(ungroup(seq)), codon_list, mode="codon")
        codon_positions[gene_name] = tmp

        for gene, values in codon_positions.items():

            for position in values["pos"]:
                val = []
                try:
                    val = chunk_from_id(values["pps"], position, int((len(x) - 1) / 2))
                except:
                    continue
                cmetagene.append(val)

    data = {"position": x, "positional_pause": cmetagene}
    save_file(data, "data_metagene_over_codon_{}".format("+".join(codon_list)))

    cmetagene = np.apply_along_axis(np.mean, 0, np.array(cmetagene))
    # cmetagene = np.sum(cmetagene, axis=0)

    plt.plot(x, cmetagene, color="black")
    plt.xlabel("Position (nt)")
    plt.ylabel("Mean positional pause score")
    plt.suptitle(
        "Pauses over {} ({} offset)".format(
            "+".join(codon_list), global_config["prefix"]
        )
    )

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "metagene_over_codon_{}".format("+".join(codon_list)),
            global_config["img_id"],
        )
    )

    plt.close()


def plot_reading_frame(framedict):
    """
    Bar plot of reading frame fractions.
    """

    x = ["Frame 1", "Frame 2", "Frame 3"]
    y = np.array([framedict[1], framedict[2], framedict[3]]) / framedict["total"]

    save_file({"frame": x, "fraction_reads": y.tolist()}, "data_global_reading_frame")

    plt.bar(x, y, edgecolor="black", color="dodgerblue")[0].set_color("#d53f04")
    plt.title("{}".format(global_config["prefix"]))
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=12)
    plt.ylabel("Read fraction", fontsize=15)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "reading_frame",
            global_config["img_id"],
        )
    )

    plt.close()


def plot_stop_peak(leftpos=10, rightpos=0):
    """
    Plot metagene for termination peak (if any)
    """
    # Perhaps use the endmaps instead of density values.
    # This would bring both the start and stop metagenes
    # on the same scale
    tmp = []
    for vector in ribosome_dentities["P"].values():
        if len(vector) < 150:
            continue

        tmp.append(vector[len(vector) - 150 :])

    tmp = np.array(tmp)
    metagene = np.sum(tmp, axis=0)
    x = np.array([i for i in range(0, 150, 1)][::-1])
    plt.plot(x, metagene, color="black")

    peakpos, delta, _ = get_offset(
        xvec=x, yvec=metagene, right_bound=rightpos, left_bound=leftpos
    )

    plt.xlabel("Nt from CDS stop", fontsize=18)
    plt.ylabel("A.U.", fontsize=18)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.xlim(150, 0)

    plt.gcf().set_size_inches(19, 5)
    plt.tight_layout()

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "stop_peak",
            global_config["img_id"],
        )
    )

    plt.close()

    data = {
        "peakpos": peakpos,
        "delta": delta,
        "positions": x.tolist(),
        "stop_pos_metagene": metagene.tolist(),
    }
    save_file(data, "data_stop_metagene")


def calculate_pause_score(site="P"):
    """
    Get the pause scores for all codons. Defaults to
    P-site offsets.
    """

    stopcodon = set(["TAA", "TAG", "TGA"])

    cds_path = global_config["coding_sequence_path"]
    cds_format = global_config["coding_sequence_format"]

    retrieve_coding_sequence(cds_path, cds_format)
    skipped_transcripts = []
    skipped_density = []
    skipped_mismatch = []
    skipped_truncated = []
    skippped_coverage = []

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
            rds = ribosome_dentities[site][gene_name]
        except:
            skipped_density.append(gene_name)
            continue

        if len(seq) != len(rds):
            skipped_mismatch.append(gene_name)
            continue

        if (len(seq) % 3 != 0) or (len(rds) % 3 != 0):
            skipped_truncated.append(gene_name)
            continue

        seq = groupN(seq, 3)  # split to codons
        rds = groupN(rds, 3)  # density per codon
        pps = []  # positional pause scores per codonn

        # Ignore gene sequences with premature stop codons
        if bool(stopcodon & set(seq[: len(seq) - 1])):
            notify(
                "{} has premature stop codon. Skipping.".format(gene_name), level="warn"
            )
            continue

        seq = seq[1 : len(seq) - 1]  # remove first and last codon
        rds = rds[1 : len(rds) - 1]

        codon_vs_density = zip(seq, rds)
        sum_codon_densities = np.sum(rds)

        if sum_codon_densities == 0:
            # This gene has no reads mapped to it. Skipping.
            continue

        for codon, density in codon_vs_density:

            positional_pause_score = (sum(density)) / (sum_codon_densities / len(rds))
            pps.append(positional_pause_score)

            codon_pauselist[site][codon].append(positional_pause_score)

        # append to global dict
        pps_vector_per_gene[site][gene_name] = {"seq": seq, "pps": pps}

    # Average out the per codon pause scores to get the final pause score
    for codon, pauselist in codon_pauselist[site].items():
        if len(pauselist) == 0:
            codon_pauselist[site][codon] = 0
        else:
            pauselist = np.array(pauselist)
            # We created our original plots without removing
            # the zero values and using np.mean()
            #
            # This was done because it recapitulates the elife Paper
            # plots. Keeping the default for now. Need to discuss further
            #
            # Then we decided to use median instead of mean.
            # But this meant we got zero values in pause scores.
            # So we decided to remove the zero values.
            # pauselist = pauselist[pauselist != 0]
            codon_pauselist[site][codon] = np.mean(pauselist)

    save_file(pps_vector_per_gene, "positional_pauses_universal")
    save_file(codon_pauselist, "per_codon_pause_scores_universal")

    notify(
        "{} gene(s) skipped because sequence is unavailable".format(
            len(skipped_transcripts)
        ),
        level="warn",
    )

    notify(
        "{} gene(s) skipped because density vector is unavailable".format(
            len(skipped_density)
        ),
        level="warn",
    )

    notify(
        "{} gene(s) skipped due to sequences-density length mismatch".format(
            len(skipped_mismatch)
        ),
        level="warn",
    )

    notify(
        "{} gene(s) skiped because length not multiple of 3".format(
            len(skipped_truncated)
        )
    )

    notify("{} gene(s) skipped due to low coverage".format(len(skippped_coverage)))

    # Write out genes skipped in pause score
    save_file(
        unique(skipped_transcripts), "skipped_gene_no_sequence_pps_{}.json".format(site)
    )
    save_file(
        unique(skipped_density), "skipped_gene_no_density_pps_{}.json".format(site)
    )
    save_file(
        unique(skipped_mismatch),
        "skipped_gene_sequence_density_length_mismatch_{}.json".format(site),
    )
    save_file(
        unique(skipped_truncated),
        "skipped_gene_not_multiple_3_pps_{}.json".format(site),
    )
    save_file(
        unique(skippped_coverage), "skipped_gene_low_coverage_pps_{}.json".format(site)
    )


def plot_pauses_combined():
    """
    Plots the per-codon and per-amino acid
    pause score heatmaps for E/P/A sites
    """

    pps_metagene = []
    for vector in pps_vector_per_gene["P"].values():
        if len(vector["pps"]) < 150:
            continue
        pps_metagene.append(vector["pps"][:150])

    xvec = [i for i in range(1, 151, 1)]
    pps_metagene = np.mean(pps_metagene, axis=0)
    plt.plot(xvec, pps_metagene, color="k")
    plt.xlabel("Codons from CDS start", fontsize=18)
    plt.ylabel("Positional pause score", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.gcf().set_size_inches(19, 5)
    plt.tight_layout()

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "pps_over_genes",
            global_config["img_id"],
        )
    )

    plt.close()

    pps_metagene_data = {"XAxis": xvec, "PPScore": pps_metagene.tolist()}

    # For now we're keeping it like this, but codons can just be a list.
    # Need to ensure that the pause scores correspond exactly
    # with the particular codon.
    codon = {"E": [], "P": [], "A": []}  # 3 are redundant, 1 should suffice
    pause = {"E": [], "P": [], "A": []}

    for site in list(pause.keys()):
        for key, value in codon_pauselist[site].items():

            if key in ["TAA", "TAG", "TGA"]:
                continue

            codon[site].append(key)  # yes
            pause[site].append(value)

    ######### HEATMAP (per codon)##########
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.1)

    sns.heatmap(
        [pause["E"]],
        xticklabels=codon["P"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="E",
        ax=ax1,
    )

    sns.heatmap(
        [pause["P"]],
        xticklabels=codon["P"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="P",
        ax=ax2,
    )

    sns.heatmap(
        [pause["A"]],
        xticklabels=codon["P"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="A",
        ax=ax3,
    )

    plt.tight_layout()

    plt.gcf().set_size_inches(19, 5)
    plt.xticks(fontsize=15)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "pause_score_codon",
            global_config["img_id"],
        )
    )

    plt.close()

    ##### BARPLOTS (per codon) #####

    fig, ax = plt.subplots(3, 1, sharex=True)
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0.1)
    fig.text(0.09, 0.4, "Pause score", ha="center", rotation="vertical", fontsize=18)

    # Plot each graph, and manually set the y tick values
    ax[0].bar(codon["E"], pause["E"], color="dodgerblue", label="E site")
    ax[0].axhline(y=np.mean(pause["E"]), color="red")
    ax[0].axhline(
        y=np.mean(pause["E"]) + np.std(pause["E"]), color="black", linestyle="dashed"
    )
    ax[0].axhline(
        y=np.mean(pause["E"]) - np.std(pause["E"]), color="black", linestyle="dashed"
    )
    ax[0].legend(prop={"size": 15})

    ax[1].bar(codon["P"], pause["P"], color="dodgerblue", label="P site")
    ax[1].axhline(y=np.mean(pause["P"]), color="red")
    ax[1].axhline(
        y=np.mean(pause["P"]) + np.std(pause["P"]), color="black", linestyle="dashed"
    )
    ax[1].axhline(
        y=np.mean(pause["P"]) - np.std(pause["P"]), color="black", linestyle="dashed"
    )
    ax[1].legend(prop={"size": 15})

    ax[2].bar(codon["A"], pause["A"], color="dodgerblue", label="A site")
    ax[2].axhline(y=np.mean(pause["A"]), color="red")
    ax[2].axhline(
        y=np.mean(pause["A"]) + np.std(pause["A"]), color="black", linestyle="dashed"
    )
    ax[2].axhline(
        y=np.mean(pause["A"]) - np.std(pause["A"]), color="black", linestyle="dashed"
    )
    ax[2].legend(prop={"size": 15})

    plt.xticks(fontsize=15, rotation=90)
    plt.gcf().set_size_inches(19, 8)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "pause_score_codon_bar",
            global_config["img_id"],
        )
    )

    plt.close()

    # Pause score for an amino acid is the average pause score of all
    # it's synonymous codons
    pause_per_aa = {"E": {}, "P": {}, "A": {}}

    for site in ["E", "P", "A"]:
        for codon, avgpause in codon_pauselist[site].items():

            AA = codon_to_aa[codon]

            if AA in pause_per_aa[site].keys():
                pause_per_aa[site][AA].append(avgpause)
            else:
                pause_per_aa[site][AA] = [avgpause]

    del pause_per_aa["E"]["*"]  # use continue in if instead?
    del pause_per_aa["P"]["*"]
    del pause_per_aa["A"]["*"]

    aa = {"E": [], "P": [], "A": []}  # ugly ugly ugly ugly ugly ...
    ps = {"E": [], "P": [], "A": []}

    for site in ["E", "P", "A"]:
        for aminoacid, pausevec in pause_per_aa[site].items():
            aa[site].append(aminoacid)
            ps[site].append(np.mean(pausevec))

    # Sort the lists for easier plotting. This is an extended sort.
    aa["E"], ps["E"] = zip(*sorted(zip(aa["E"], ps["E"])))
    aa["P"], ps["P"] = zip(*sorted(zip(aa["P"], ps["P"])))
    aa["A"], ps["A"] = zip(*sorted(zip(aa["A"], ps["A"])))

    aa_pause_data = {"AA": aa, "Pause": ps}

    if (aa["E"] != aa["P"]) or (aa["E"] != aa["A"]):
        notify("Amino acid pause: AA order mismatch", level="warn")

    #### HEATMAP (amino acids)########
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.1)

    sns.heatmap(
        [ps["E"]],
        xticklabels=aa["E"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="E",
        ax=ax1,
    )

    sns.heatmap(
        [ps["P"]],
        xticklabels=aa["E"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="P",
        ax=ax2,
    )

    sns.heatmap(
        [ps["A"]],
        xticklabels=aa["E"],
        cmap="RdBu_r",
        center=1,
        vmax=2.5,
        vmin=-0.5,
        yticklabels="A",
        ax=ax3,
    )

    plt.tight_layout()

    plt.gcf().set_size_inches(19, 5)
    plt.xticks(fontsize=20)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "pause_score_aa",
            global_config["img_id"],
        )
    )

    plt.close()

    ######## BARPLOTS (per amino acid) ########################
    fig, ax = plt.subplots(3, 1, sharex=True)
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0.1)
    fig.text(0.09, 0.4, "Pause score", ha="center", rotation="vertical", fontsize=18)
    # Plot each graph, and manually set the y tick values

    ax[0].bar(aa["E"], ps["E"], color="dodgerblue", label="E site")
    ax[0].axhline(y=np.mean(ps["E"]), color="red")
    ax[0].axhline(
        y=np.mean(ps["E"]) + np.std(ps["E"]), color="black", linestyle="dashed"
    )
    ax[0].axhline(
        y=np.mean(ps["E"]) - np.std(ps["E"]), color="black", linestyle="dashed"
    )
    ax[0].legend(prop={"size": 12})

    ax[1].bar(aa["P"], ps["P"], color="dodgerblue", label="P site")
    ax[1].axhline(y=np.mean(ps["P"]), color="red")
    ax[1].axhline(
        y=np.mean(ps["P"]) + np.std(ps["P"]), color="black", linestyle="dashed"
    )
    ax[1].axhline(
        y=np.mean(ps["P"]) - np.std(ps["P"]), color="black", linestyle="dashed"
    )
    ax[1].legend(prop={"size": 12})

    ax[2].bar(aa["A"], ps["A"], color="dodgerblue", label="A site")
    ax[2].axhline(y=np.mean(ps["A"]), color="red")
    ax[2].axhline(
        y=np.mean(ps["A"]) + np.std(ps["A"]), color="black", linestyle="dashed"
    )
    ax[2].axhline(
        y=np.mean(ps["A"]) - np.std(ps["A"]), color="black", linestyle="dashed"
    )
    ax[2].legend(prop={"size": 12})

    plt.xticks(fontsize=18, rotation=0)
    plt.gcf().set_size_inches(19, 8)

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "pause_score_aa_bar",
            global_config["img_id"],
        )
    )

    plt.close()

    save_file(pps_metagene_data, "pps_metagene")
    save_file(aa_pause_data, "per_aa_pause_scores")


def read_length_histogram():
    """
    Plot a histogram of the read lengths
    """

    readlens = []
    bamfile = global_config["bamfile"]

    for gene, details in annotation.items():

        start = details["start"]
        stop = details["stop"]
        strand = details["strand"]
        refer = details["refname"]

        for read in bamfile.fetch(start=start, stop=stop, reference=refer):

            if strand == "+" and read.is_reverse:
                continue
            if strand == "-" and not read.is_reverse:
                continue
            if read.reference_length > 50:
                notify(
                    "Histogram: found read of length > 50nt",
                    level="warn",
                    onetime=True,
                )
                continue

            readlens.append(read.reference_length)

    # bins = [i for i in range(20,41,1)
    plt.hist(readlens, edgecolor="black", color="#1A9BE8")
    # plt.xticks(bins)
    plt.title("Read length distribution")
    plt.xlabel("Read length", fontsize=15)
    plt.ylabel("Frequency", fontsize=15)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.title("{}".format(global_config["prefix"]))

    plt.tight_layout()

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "read_length_histogram",
            global_config["img_id"],
        )
    )

    plt.close()

    save_file(readlens, "read_lengths")


def read_terminal_stats():
    """
    Plot the fraction of bases at the 5' and 3'
    ends of the reads
    """

    terminal_basecount = {
        "totalreads": 0,
        "5": {"A": 0, "G": 0, "T": 0, "C": 0, "N": 0},
        "3": {"A": 0, "G": 0, "T": 0, "C": 0, "N": 0},
    }

    bamfile = global_config["bamfile"]

    for gene_name, details in annotation.items():

        start = details["start"]
        stop = details["stop"]
        strand = details["strand"]
        refer = details["refname"]

        for read in bamfile.fetch(start=start, stop=stop, reference=refer):

            if strand == "+" and read.is_reverse:
                continue
            if strand == "-" and not read.is_reverse:
                continue
            if read.reference_length > 50:
                continue

            seq = read.query_sequence

            if strand == "+":
                fp_base = seq[0]
                tp_base = seq[len(seq) - 1]
                terminal_basecount["5"][fp_base] += 1
                terminal_basecount["3"][tp_base] += 1
            if strand == "-":
                fp_base = seq[len(seq) - 1]
                tp_base = seq[0]
                terminal_basecount["5"][fp_base] += 1
                terminal_basecount["3"][tp_base] += 1

            terminal_basecount["totalreads"] += 1

    save_file(terminal_basecount, "terminal_base_fraction")

    xaxis = ["A", "T", "G", "C"]

    # Calculate the fraction of total
    A5 = terminal_basecount["5"]["A"] / terminal_basecount["totalreads"]
    G5 = terminal_basecount["5"]["G"] / terminal_basecount["totalreads"]
    T5 = terminal_basecount["5"]["T"] / terminal_basecount["totalreads"]
    C5 = terminal_basecount["5"]["C"] / terminal_basecount["totalreads"]

    A3 = terminal_basecount["3"]["A"] / terminal_basecount["totalreads"]
    G3 = terminal_basecount["3"]["G"] / terminal_basecount["totalreads"]
    T3 = terminal_basecount["3"]["T"] / terminal_basecount["totalreads"]
    C3 = terminal_basecount["3"]["C"] / terminal_basecount["totalreads"]

    y3 = [A3, T3, G3, C3]
    y5 = [A5, T5, G5, C5]

    # set width of bar
    barWidth = 0.25

    # Set position of bar on X axis
    r1 = np.arange(len(y3))
    r2 = [x + barWidth for x in r1]

    # Make the plot
    plt.bar(r1, y3, color="#1A9BE8", width=barWidth, edgecolor="black", label="3' end")
    plt.bar(r2, y5, color="#E83F1A", width=barWidth, edgecolor="black", label="5' end")

    # Add xticks on the middle of the group bars
    plt.xlabel("Base at read terminal", fontsize=15)
    plt.ylabel("Read fraction", fontsize=15)
    plt.xticks([r + barWidth for r in range(len(y3))], xaxis)
    # Create legend & Show graphic
    plt.legend()

    plt.title("Terminal read fractions ({})".format(global_config["prefix"]))
    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "basestat",
            global_config["img_id"],
        )
    )

    plt.close()


def calculate_asymmetry_scores():
    """
    Asymmetry score per gene =

           log2(sum(densities_2nd_half)/sum(densities_1st_half))
    """

    def gethalf(vector):
        pivot = int(np.floor(len(vector) / 2))
        return (vector[:pivot], vector[pivot:])

    # Dictionary containing genename --> asymmetry score
    asymmetry = {}
    tmp_asc = []  # temporarily store scores in a list for box plot

    for gene_name, densities in ribosome_dentities["P"].items():

        if len(densities) < 100:
            notify(
                "asymmetry score: {} too short. Skipped.".format(
                    gene_name, level="warn"
                )
            )
            continue

        # Ignore the first and last 50nt of the gene
        density_5p, density_3p = gethalf(densities[50 : len(densities) - 50])
        rho_sum_5p = sum(density_5p)
        rho_sum_3p = sum(density_3p)

        if (rho_sum_3p == 0) or (rho_sum_5p == 0):
            notify(
                "asymmetry score: No density profile for {}. Skipped.".format(
                    gene_name
                ),
                level="warn",
            )
            continue

        asymmetry_score = np.log2(rho_sum_3p / rho_sum_5p)
        asymmetry[gene_name] = asymmetry_score
        tmp_asc.append(asymmetry_score)

    save_file(asymmetry, "asymmetry_scores")

    colors = ["dodgerblue"]
    box = plt.boxplot(tmp_asc, vert=False, notch=True, patch_artist=True)

    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)

    # plt.xlabel("log2(sum(#rpfs 3' half)/sum(#rpfs 5' half)")
    plt.xlabel("Asymmetry score", fontsize=15)
    plt.xticks(fontsize=15)
    plt.title("Asymmetry score distribution ({})".format(global_config["prefix"]))
    plt.axvline(x=0, color="red", linestyle="--")
    plt.gcf().set_size_inches(6, 4)
    plt.tight_layout()

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "asymmetry_score",
            global_config["img_id"],
        )
    )

    plt.close()


def consistency_score_per_transcript():
    """
    Scatterplot of readcounts mapping to the second half
    of a gene vs the first.
    NOTE: Taking the whole gene, not trimming anything
    """

    bamfile = global_config["bamfile"]
    five_prime_counts = []
    three_prime_counts = []

    consistency_dict = {}  # temporarily hold the data

    for gene_name, details in annotation.items():

        count_5p, count_3p = 0, 0

        start = details["start"]
        stop = details["stop"]
        strand = details["strand"]
        refer = details["refname"]

        pivot = int(np.ceil((stop - start) / 2)) + start

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

        if (count_5p == 0) or (count_3p == 0):
            continue

        five_prime_counts.append(count_5p)
        three_prime_counts.append(count_3p)

        consistency_dict[gene_name] = {
            "5_prime_half": count_5p,
            "3_prime_half": count_3p,
        }

    five_prime_counts = np.log2(np.array(five_prime_counts))  # y
    three_prime_counts = np.log2(np.array(three_prime_counts))  # x
    pearson = np.corrcoef(three_prime_counts, five_prime_counts)[1, 0]

    m, c = np.polyfit(
        three_prime_counts, five_prime_counts, 1  # 1st degree polynomial fit
    )

    consistency_dict["metadata"] = {
        "pearson_coeff": pearson,
        "slope": m,
        "intercept": c,
    }

    save_file(consistency_dict, "consistency data")

    regression_function = np.poly1d(m, c)

    plt.plot(
        three_prime_counts,  # x
        five_prime_counts,  # y
        "bo",
        three_prime_counts,  # x
        regression_function(three_prime_counts),  # fx
        "-k",
        markersize=1,
        alpha=0.5,
    )

    plt.title(
        "y = ({})x + {}, Pearson's R = {}".format(
            round(m, 2), round(c, 2), round(pearson, 2)
        )
    )

    plt.xlabel("log2(sum(read count 3' half))")
    plt.ylabel("log2(sum(read count 5' half))")
    plt.suptitle("{}".format(global_config["prefix"]))

    plt.savefig(
        "{}{}_{}_{}.svg".format(
            global_config["img_dir"],
            global_config["session_id"],
            "consistency_plot",
            global_config["img_id"],
        )
    )

    plt.close()


def generate_metagene_vector():
    """
    Computes the metagene by first normalizing and then taking a
    column wise sum of the endmap vectors for each gene
    """

    up = global_config["endmap_upstream"]
    dn = global_config["endmap_downstream"]
    lengths = global_config["readlengths"]

    tmp_list = []

    for length in lengths:
        for gene_name in annotation.keys():

            try:
                # chop out unwanted values
                vector = endmap_vectors["P"][gene_name][length][: (up + dn)]
            except:
                # Gene doesn't have RPFs of a particular length
                continue

            tmp_list.append(normalize_counts(vector))

        if len(tmp_list) == 0:
            # error, needs to be handled downstream
            metagene_per_readlength[length] = [-1]
        else:
            metagene_per_readlength[length] = np.apply_along_axis(
                sum, 0, np.array(tmp_list)
            )

        del tmp_list[:]


def retrieve_coding_sequence(fname, ftype="fasta"):
    """
    Loads a file mapping genes to their CDS sequences
    """
    if len(transcripts_dict) == 0:
        cds = list(SeqIO.parse(fname, ftype))
        for i in range(len(cds)):
            transcripts_dict[cds[i].id] = str(cds[i].seq)


def calculate_densities_over_genes():
    """
    Compute the ribosome E/P/A site densities over all genes
    """

    up = global_config["endmap_upstream"]
    dn = global_config["endmap_downstream"]
    lengths = global_config["readlengths"]

    skipped = []
    notify("{} gene(s) have endmap_vectors".format(len(endmap_vectors["P"].keys())))

    for gene, endmaps_per_readlength in endmap_vectors["P"].items():
        # WARNING: These values are not normalized by any means.
        # Normalize for length before any calculations can be done.
        tmp = []  # store the endmaps per readlength.
        for length, endmap in endmaps_per_readlength.items():
            if length not in lengths:
                continue  # only calculate for the particular set of RPFs
            tmp.append(endmap)

        if len(tmp) == 0:
            notify(
                "{}: no coverage data, can't determine densities".format(gene),
                level="warn",
            )
            continue

        tmp = np.array(tmp)
        tmp = np.apply_along_axis(sum, 0, tmp)
        # trimming out the flanking regions we added
        tmp = tmp[up : len(tmp) - dn]
        ribosome_dentities["P"][gene] = tmp.tolist().copy()
        tmp = tmp.tolist()
        del tmp[:]

    for gene, endmaps_per_readlength in endmap_vectors["E"].items():

        tmp = []  # store the endmaps per readlength.
        for length, endmap in endmaps_per_readlength.items():
            if length not in lengths:
                continue  # only calculate for the particular set of RPFs
            tmp.append(endmap)

        if len(tmp) == 0:
            notify(
                "{}: no coverage data, can't determine densities".format(gene),
                level="warn",
            )
            continue

        tmp = np.array(tmp)
        tmp = np.apply_along_axis(sum, 0, tmp)
        # trimming out the flanking regions we added
        tmp = tmp[up : len(tmp) - dn]
        ribosome_dentities["E"][gene] = tmp.tolist().copy()
        tmp = tmp.tolist()
        del tmp[:]

    for gene, endmaps_per_readlength in endmap_vectors["A"].items():

        tmp = []  # store the endmaps per readlength.
        for length, endmap in endmaps_per_readlength.items():
            if length not in lengths:
                continue  # only calculate for the particular set of RPFs
            tmp.append(endmap)

        if len(tmp) == 0:
            notify(
                "{}: no coverage data, can't determine densities".format(gene),
                level="warn",
            )
            continue

        tmp = np.array(tmp)
        tmp = np.apply_along_axis(sum, 0, tmp)
        # trimming out the flanking regions we added
        tmp = tmp[up : len(tmp) - dn]
        ribosome_dentities["A"][gene] = tmp.tolist().copy()
        tmp = tmp.tolist()
        del tmp[:]

    # Taking P as representarive for E and A but this function should be rewritten
    notify("calculated densities for {} gene(s)".format(len(ribosome_dentities["P"])))
    save_file(unique(skipped), "skipped_genes_no_coverage_density.json")


def main():
    print("----------------------------------------------------------------")
    print("               {}\n".format(versionstr))
    print("      Copyright (C) 2020, Molecular Genetics Laboratory")
    print("         National Institute of Immunology, New Delhi\n")
    print(" This program is distributed in the hope that it will be useful,")
    print(" but WITHOUT ANY WARRANTY; without even the implied warranty of")
    print("  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the")
    print("         GNU General Public License for more details.")
    print("----------------------------------------------------------------\n")

    try:
        set_config(load_json(sys.argv[1]))
    except:

        fh = open("bard_config_template.json", "w")
        ft = open("bard_config_help.txt", "w")
        fh.write(CONFIG_TMP)
        ft.write(CONFIG_HELP)
        fh.close()
        ft.close()

        print(
            "\033[33m       ABORT: Malformed or non-existant configuration file\033[m\n"
        )
        print("  Usage: ")
        print("         ~$ python bard.py config.json\n")
        print("  Try using the templates if you need help with the configuration\n")
        print("  Bye!\n")

        raise SystemExit()

    print("\n")
    script_init()

    check_gene_list()

    notify("Determining terminal base fractions")
    read_terminal_stats()

    notify("Found {} genes in annotation file".format(len(annotation)))

    notify("Plotting read length histogram")
    read_length_histogram()

    notify(
        "Filtering low coverage genes (cutoff={}, metric={})".format(
            global_config["coverage_cutoff"], global_config["coverage_metric"]
        )
    )
    get_max_coverage_genes()
    notify("Reduced number of genes to {}".format((len(high_coverage_genes))))

    notify("Scanning for aligned read terminals (no offset)")
    map_gene_to_endmaps(apply_offset=False)

    notify("Generating metagene vector")
    generate_metagene_vector()

    notify("Determining offsets")
    offset_delta_for_metagene()

    notify("Plotting metagene per read length (without offsets)")
    plot_metagene_per_readlength()

    notify("Mapping genes to endmaps (with offsets)")
    map_gene_to_endmaps(
        apply_offset=True, only_reads_with_codons={"E": [], "P": [], "A": []}
    )

    notify("Calculating ribosome densities")
    calculate_densities_over_genes()

    notify("Calculating pause scores")
    calculate_pause_score(site="P")
    calculate_pause_score(site="E")
    calculate_pause_score(site="A")

    notify("Plotting pause scores")
    plot_pauses_combined()

    notify("Plotting termination peak")
    plot_stop_peak()

    notify("Calculating asymmetry scores")
    calculate_asymmetry_scores()

    # notify("Consistency score per transcript is disabled")
    # consistency_score_per_transcript()

    notify("Metagene profile over codons is disabled")
    # metagene_over_codon(aa_to_codon["A"], mode="combined") # serine
    # metagene_over_codon([aa_to_codon["S"], aa_to_codon["A"]], mode="composite")
    try:
        cm_mode = global_config["codon_metagene_mode"]
        codonlist = global_config["codon_metagene_list"]
        metagene_over_codon(expand_to_codon(codonlist), mode=cm_mode)
    except:
        notify("No metagene for codons specified, skipping", level="warn")

    notify("Determining reading frame")
    plot_reading_frame(reading_frame(ignore_read_by_nt={5: [], 3: []}))

    notify("Generating metagene")
    generate_metagene_vector()

    save_file(ribosome_dentities, "ribosome_dentities")

    notify("Plotting metagene per read length (with offsets)")
    plot_metagene_per_readlength()

    notify("Plotting initiation peak")
    plot_initiation_peak(peak=True, peak_range=[5, 25])

    notify("Coverage profile is disabled")
    generate_coverage_profiles(disabled=True)

    notify("Cleaning up")
    global_config["bamfile"].close()
    del global_config["bamfile"]  # can't be saved as plain text
    # Save the configuration for reference
    save_json(global_config, global_config["conf_file"])

    notify("Done!")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        global_config["bamfile"].close()
        notify("Keyboard interrupt. Aborting.", level="crit", fatal=True)
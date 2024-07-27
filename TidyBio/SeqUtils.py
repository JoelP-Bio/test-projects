import sys
import pandas as pd
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv(r"C:\Users\JJCuco\PycharmProjects\pythonProject1\TidyBio\amino_acids_codontable.csv")
df.rename(columns={"Unnamed: 2":"Sym"},inplace=True)
aa3_to1_dict = dict(zip(df['Symbols'], df['Sym']))
full_aa_codon_dict = dict(zip(df['Amino acids'], df['Codons']))
full_aato1_dict = dict(zip(df['Amino acids'],df['Sym']))
CodonTable = {
            # 'M' - START, '*' - STOP
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TGT": "C", "TGC": "C",
            "GAT": "D", "GAC": "D",
            "GAA": "E", "GAG": "E",
            "TTT": "F", "TTC": "F",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            "CAT": "H", "CAC": "H",
            "ATA": "I", "ATT": "I", "ATC": "I",
            "AAA": "K", "AAG": "K",
            "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATG": "M",
            "AAT": "N", "AAC": "N",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TGG": "W",
            "TAT": "Y", "TAC": "Y",
            "TAA": "*", "TAG": "*", "TGA": "*"
        }
class Sequence(object):
    # Create a Valid Sequence for DNA, RNA

    # example: seq1 = Sequence('ATGC')

    def __init__(self, seq):
        super(Sequence, self).__init__()
        self.seq = seq
        if not isinstance(self.__validate_seq(seq),str):
            raise TypeError("The sequence data given to a Sequence"
                            " object should be a string (not another Sequence object)"
                            " nor a Non-valid Nucleotide [A, T, G, C, U]")

    def __repr__(self):
        return "Sequence(seq='{}')".format(self.seq)

    def __str__(self):
        return self.seq

    def __validate_seq(self, seq):
        real_seq, base_nucleotide = self.__validate_seq_type(seq)
        for base in real_seq:
            if base not in base_nucleotide:
                return False
        return real_seq

    def __validate_seq_type(self, seq):
        real_seq = seq.upper()
        if "U" in real_seq:
            for base in real_seq:
                if "T" in real_seq:
                    raise "NucleotideError: {} not an RNA nucleotide ['A,U,G,C']".format(base)
                else:
                    base_nucleotide = ["A", "U", "G", "C"]
                    return real_seq, base_nucleotide
        else:
            for base in real_seq:
                if "U" in real_seq:
                    raise "NucleotideError: {} not a DNA nucleotide ['A,T,G,C']".format(base)
                else:
                    base_nucleotide = ["A", "T", "G", "C"]
                    return real_seq, base_nucleotide

    def __len__(self):
        return len(self.seq)

    def __contains__(self, sub_char):
        return sub_char in str(self)

    def __getitem__(self, key):
        if isinstance(key, slice):
            # Get the start, stop, and step from the slice
            vals = [self[ii] for ii in range(*key.indices(len(self)))]
            return ''.join(vals)
        elif isinstance(key, int):
            if key < 0:  # Handle negative indices (indexes from the right)
                key += len(self)
            if key < 0 or key >= len(self):
                raise IndexError("The index {} is out of range.".format(key))
            return self.seq[key]  # Get the data from elsewhere
        else:
            raise TypeError("Invalid argument type.")

    # Basic Fxn
    # Count, find, index
    def count(self, subseq, start=0, end=sys.maxsize):
        # Returns the count of a nucleotide in a sequence
        return str(self).count(subseq, start, end)

    def find(self, subseq, start=0, end=sys.maxsize):
        # Returns the position of a nucleotide in a sequence
        return str(self).find(subseq, start, end)

    def rfind(self, subseq, start=0, end=sys.maxsize):
        # Returns the position of a nucleotide in a sequence from the right
        return str(self).rfind(subseq, start, end)

    def index(self, subseq, start=0, end=sys.maxsize):
        # Returns the Index/Position of a nucleotide in a sequence
        return str(self).index(subseq, start, end)

    def rindex(self, subseq, start=0, end=sys.maxsize):
        # Returns the Index/Position of a nucleotide in a sequence from the right
        return str(self).rindex(subseq, start, end)

    # Main Fxn
    def get_symbol_frequency(self):
        # Returns the Frequency of a Nucleotide in a Sequence
        if "U" in self.seq:
            base_dict = {"A": 0, "U": 0, "G": 0, "C": 0}
            for base in self.seq:
                if self.__validate_seq(base):
                    if "T" in self.seq:
                        raise "NucleotideError: {} not an RNA nucleotide ['A,U,G,C']".format(base)
                    else:
                        base_dict[base] += 1
        else:
            base_dict = {"A": 0, "T": 0, "G": 0, "C": 0}
            for base in self.seq:
                if self.__validate_seq(base):
                    if "U" in self.seq:
                        raise "NucleotideError: {} not a DNA nucleotide ['A,T,G,C']".format(base)
                    else:
                        base_dict[base] += 1
        return base_dict
    def get_symbol_percentage(self, chars):
        # Returns the Percentage of a Nucleotide in a Sequence
        freq_list = []
        res_dict = {}
        char_list = chars
        for i, char in enumerate(char_list):
            freq_list.append(self.count(char) / len(self))
            res_dict[char] = freq_list[i]
        return res_dict

    def gc_content(self):
        # Returns the gc_content of a sequence as a percentage
        result = float(str(self).count('G') + str(self).count('C')) / len(self) * 100
        return result

    def at_content(self):
        # Returns the at_content of a sequence as a percentage
        result = float(str(self).count('A') + str(self).count('T')) / len(self) * 100
        return result

    def complement(self):
        # Returns the complementary strand of a sequence
        comp_pairs = []
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        for a in self:
            if a in base_pairs.keys():
                comp_pairs.append(base_pairs[a])
        return "".join(comp_pairs)

    def reverse_complement(self):
        # Returns the reverse complementary strand of a sequence
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        comp_pairs = [base_pairs[a] for a in seq if a in base_pairs.keys()]
        reverse_pairs = "".join(comp_pairs)[::-1]
        return reverse_pairs
                
    def transcribe(self):
        # Transcribes Sequence into mRNA
        mrna_result = self.seq.replace("T", "U")
        return mrna_result
                
    def translate(self, start_pos=0):
        # Translates Sequence into Protein/Amino Acids
        amino_acids_list =[CodonTable[self.seq[pos:pos+3]] for pos in range(start_pos, len(self.seq)-2, 3)]
        return "".join(amino_acids_list)

def reverse_transcribe(mRNA):
    # Reverse transcribes mRNA to DNA
    dna_result = mRNA.replace("U", "T")
    return Sequence(dna_result)

def reverse_translate(protein):
    # Reverse translates Protein/Amino Acids to mRNA
    aa_list = []
    for i in range(len(protein)):
        aa_list.append(protein[i])
    codon_list = [get_key(aa, CodonTable) for aa in aa_list]
    return Sequence("".join(codon_list))
                
def get_key(val, my_dict):
    # Returns the key of a value 
    for key, value in my_dict.items():
        if val == value or val in value:
            return key

def get_value(val, my_dict):
    # Returns the value of a key
    for key, value in my_dict.items():
        if val == key:
            return value
                    
def count_kmers(seq, k=3):
    # Empty Dictionary
    counts = {}
    num_kmers = len(seq) - k + 1
    for i in range(num_kmers):
        # Slicing the seq to get the kmer
        kmer = seq[i:i+k]  # increment indices
        # makes sure kmers are not being iterated again
        if kmer not in counts:
            counts[kmer] = 0  # set to 0
            # increment
        counts[kmer] += 1  # adds 1 to the 0
    return counts
            
def get_kmers(seq, k=3):
    # Empty Dictionary
    counts = {}
    kmer_list = []
    num_kmers = len(seq) - k + 1
    for i in range(num_kmers):
        # Slicing the seq to get the kmer
        kmer = seq[i:i+k]  # increment indices
        # makes sure kmers are not being iterated again
        kmer_list.append(kmer)
    return kmer_list

def kmer_dist(seq1, seq2):
    seq1_k = get_kmers(seq1)
    seq2_k = get_kmers(seq2)
    seq1_set = set(seq1_k)
    seq2_set = set(seq2_k)
    union_seq = seq1_set.union(seq2_set)
    # intersection shows identical elements in either set
    intersection_seq = seq1_set.intersection(seq2_set)
    sym_difference = len(union_seq) - len(intersection_seq)
    # dissimilarites show differing elements in either set
    dissimilarities = seq1_set ^ seq2_set
    print(dissimilarities)
    print(sym_difference)
    distance = sym_difference / len(union_seq)
    return distance
  
def get_codons(seq, k=3):
    # returns codons
    codon_list = []
    for i in range(0, len(seq), k):
        codon_list.append(str(seq)[i:i+k])
    return codon_list
            
def convert1to3(seq):
    # returns 3-letter protein seq
    term_list = []
    for i in seq:
        res = get_key(i, aa3_to1_dict)
        term_list.append(res)
    return "".join(term_list)
            
def convert3to1(seq):
    # returns 1-letter protein seq
    term_list = []
    for i in get_codons(seq, k=3):
        res = get_value(i, aa3_to1_dict)
        term_list.append(res)
    return ''.join(term_list)
            
def hamming_distance(lhs,rhs):
    # returns hamming distance based on mismatches
    return len([(x,y) for x,y in zip(lhs,rhs) if x != y])
            
def occurrence(main_seq, sub_seq):
    # shows the occurrence of a sub sequence with its positions
    indices = []
    start = 0
    # while loop for multiple instances of the sub sequence within the main
    while True:
        start = main_seq.find(sub_seq, start)
        if start > 0:
            end = start + len(sub_seq)
            indices.append([start, end])
        else:
            break
        start += 1
    return indices
            
def get_acid_name(seq):
    # returns full acid name
    term_list = []
    for i in get_codons(seq):
        res = get_key(i, full_aa_codon_dict)
        term_list.append(res)
    return "".join(term_list)
            
def codon_frequency(seq, aminoacid):
    # returns frequency of each codon
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if CodonTable[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalScore = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalScore, 2)
    return freqDict
            
def delta(x, y):
    # matches are valued as 0
    return 0 if x == y else 1
            
def M(seq1, seq2, i, j, k):
    # returns sum of each coordinate (x, y)
    return sum(delta(x, y) for x, y in zip(seq1[i:i+k], seq2[j:j+k]))
            
def makeMatrix(seq1, seq2, k):
    n = len(seq1)
    m = len(seq2)
    # returns nested list of summation values for each coordinate for pltMatrix() to read
    return [[M(seq1, seq2, i, j, k) for j in range(m)] for i in range(n)]
            
def plotMatrix(M, t, seq1, seq2, nonblank=chr(0x25A0), blank = ' '):
    # plots matrix
    seq1 = str(seq1)
    seq2 = str(seq2)
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label, row in zip(seq1, M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)
                
def dotplot(seq1, seq2, k = 1, t = 1):
    # finishes dotplot
    M = makeMatrix(seq1, seq2, k)
    plotMatrix(M, t, seq1, seq2)
            
def dotplotx(seq1, seq2):
    # dotplot utilizing numpy and matplotlib
    plt.imshow(np.array(makeMatrix(seq1, seq2, 1)))
    # on x-axis list all sequences of seq 2
    plt.xticks(np.arange(len(list(seq2))), list(seq2))
    # on y-axis list all sequences of seq 1
    plt.yticks(np.arange(len(list(seq1))), list(seq1))
    plt.show()

def aa_type_counter(protein):
    zero = '0' * 20
    new_acid_dict = dict(zip(full_aato1_dict.values(),[int(i) for i in list(zero)]))
    for aa in protein:
        for acid in new_acid_dict.keys():
            if aa == acid:
                new_acid_dict[acid] += 1
    return new_acid_dict
            
def aa_most_common(protein, n):
    # returns most common proteins
    freq = Counter(protein).most_common(n)
    return freq

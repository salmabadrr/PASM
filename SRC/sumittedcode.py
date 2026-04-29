import math
import os
import numpy as np  # dealing with numpy ararys instead of lists of lists is more memory efficient
from math import log
import itertools
import math
import random
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt
#normally the data is retrived from NCBI as one fasta file, This code works on data as a folder of individual fasta files, so here is the function to convert the fasta file to a fodler
def file_to_folder(file_path, folder_name):
    path = "D:/Downloads"
    os.chdir(path)
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    os.chdir(f"{path}/{folder_name}")
    sequence = ''
    seqfile = None
    file = open(file_path)
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if seqfile is not None:
                seqfile.write(sequence + "\n")
                seqfile.close()
            sequence = ''
            x = line.split(' ')
            if not os.path.exists(x[0][1:]):
                seqfile = open(f'{x[0][1:]}.FASTA', 'w')
                seqfile.write(f"{line}\n")
        else:
            sequence += line
    #write last sequence
    if seqfile is not None:
        seqfile.write(sequence + "\n")
        seqfile.close()
    return len(os.listdir(f"{path}/{folder_name}"))
# file_path = "D:/Downloads/FASTA (2).fa"
# print(file_to_folder(file_path,'NA_data_fulllength2021'))

# ---------------------------------------------------------------------------------------------------------------------------
# Needleman-Wunsch Alignment
class Global_Align():
    def __init__(self, seq1, seq2):
        self.init_mat = []
        self.seq1 = seq1
        self.seq2 = seq2
        self.n = len(self.seq1)
        self.m = len(self.seq2)
        self.match = +1
        self.mismatch = -1
        self.gap = -2

    def intialize_matrix(self):  # create a matrix with n+1 rows and m+1 cols
        self.init_mat = np.zeros((self.n + 1, self.m + 1), dtype=int)
        return self.init_mat

    def Fill_matrix(self):
        # intialize first row and col 'only gaps'
        for i in range(1, self.n + 1):  # iterate over first col and fill
            self.init_mat[i][0] = i * self.gap
        for j in range(1, self.m + 1):  # iterate over first row and fill
            self.init_mat[0][j] = j * self.gap
        # fill matrix
        for i in range(1, self.n + 1):
            for j in range(1, self.m + 1):
                if self.seq1[i - 1] == self.seq2[j - 1]:  # match
                    self.init_mat[i][j] = max(self.init_mat[i][j - 1] + self.gap,
                                              self.init_mat[i - 1][j] + self.gap,
                                              self.init_mat[i - 1][j - 1] + self.match)
                else:  # mismatch
                    self.init_mat[i][j] = max(self.init_mat[i][j - 1] + self.gap,
                                              self.init_mat[i - 1][j] + self.gap,
                                              self.init_mat[i - 1][j - 1] + self.mismatch)

        return self.init_mat

    def Trace_back(self):# bottom right cell
        self.seq1_align = []
        self.seq2_align = []
        i = self.n
        j = self.m
        while (i > 0 or j > 0):  # both sequences have characters left
            if i > 0 and j > 0 and self.init_mat[i][j] == self.init_mat[i - 1][j - 1] + (
                    self.match if self.seq1[i - 1] == self.seq2[
                        j - 1] else self.mismatch):
                self.seq1_align.append(self.seq1[i - 1])
                self.seq2_align.append(self.seq2[j - 1])
                # we moved diagonally a cell up in the seqeunce
                i -= 1
                j -= 1
            # check the horizontal value: a gap in sequence 1
            elif j > 0 and self.init_mat[i][j] == self.init_mat[i][j - 1] + self.gap:
                self.seq1_align.append('-')
                self.seq2_align.append(self.seq2[j - 1])
                j -= 1
            # check vertical valuee: a gap in sequence 2
            elif i > 0 and self.init_mat[i][j] == self.init_mat[i - 1][j] + self.gap:
                self.seq1_align.append(self.seq1[i - 1])
                self.seq2_align.append('-')
                i -= 1

            # fill the remaining gaps
        while i > 0:  # gaps in sequence 2
            self.seq1_align.append(self.seq1[i - 1])
            self.seq2_align.append('-')
            i -= 1
        while j > 0:  # gaps in sequence 1
            self.seq1_align.append('-')
            self.seq2_align.append(self.seq2[j - 1])
            j-=1
        # the lists are filled backward, so reverse them
        self.seq1_align.reverse()
        self.seq2_align.reverse()
        # create match string
        self.match_string = ''
        for k in range(len(self.seq1_align)):
            if self.seq1_align[k] == self.seq2_align[k]:
                self.match_string += "|"
            elif self.seq1_align[k] != self.seq2_align[k]:
                if (self.seq1_align[k] == "-" or self.seq2_align[k] == "-"):
                    self.match_string += " "
                else:
                    self.match_string += "*"
        return self.seq1_align, self.match_string, self.seq2_align

    def Align_Output_file(self, name1, name2): # afunction if we wanto to save the pairwise alignment output
        os.chdir("D:/Junior/Programming-BMS 321/Project/pythonProject/Last")
        file_name = f'{name1} vs {name2}.FASTA'
        if not os.path.exists(file_name):
            with open(file_name, 'w') as output:
                output.write(f">{name1}\n{''.join(self.seq1_align)}\n")
                output.write(f">{name2}\n{''.join(self.seq2_align)}\n")

    def transition_transversion(self):
        self.p = 0  # transition
        self.q = 0  # transversion
        self.Transition = [('A', 'G'), ('G', 'A'),
                           ('T', 'C'), ('C', 'T'), ]

        self.Transversion = [('A', 'C'), ('A', 'T'),
                             ('G', 'C'), ('G', 'T'),
                             ('C', 'A'), ('C', 'G'),
                             ('T', 'A'), ('T', 'G')]
        # If the two nucleotides don't match, it is either transition or transversion
        total = 0
        for u in self.match_string:
            if u == '|' or u == "*":
                total += 1
        for k in range(len(self.seq1_align)):
            if self.seq1_align[k] != '-' and self.seq2_align[k] != '-':
                if self.seq1_align[k] != self.seq2_align[k]:
                    pair = (self.seq1_align[k], self.seq2_align[k])
                    if pair in self.Transition:
                        self.p += 1
                    elif pair in self.Transversion:
                        self.q += 1
        self.p = self.p / total
        self.q = self.q / total
        return self.p, self.q

    def calculate_alignment(self):
        match_count = 0
        mismatch_count = 0
        gap_count = 0
        for i in self.match_string:
            if i == "|":
                match_count += 1
            elif i == '*':
                mismatch_count += 1
            elif i == ' ':
                gap_count += 1
        self.alignment_score = (match_count * self.match) + (mismatch_count * self.mismatch) + (
                    gap_count * self.gap)
        self.percent_identity = (match_count / len(self.match_string)) * 100
        self.Distance = 1 - (match_count / len(self.match_string))
        return round(self.percent_identity, 2), round(self.Distance, 5), self.alignment_score

    #
    def evolutionary_distance_k80(self):
        term1 = -0.5 * math.log(1 - (2 * self.p) - self.q)
        term2 = -0.25 * math.log(1 - (2 * self.q))
        Kimura_distance = term1 - term2
        return round(Kimura_distance, 5)

def read_files(path): #the function reads a fasta file and returns a sequence
    sequence = ''
    file = open(path)
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith(">") == False:
            sequence += line
        else:
            continue
    return sequence

def Construct_Distance_Matrix():
    dataset_path = "D:/Junior/Programming-BMS 321/Project/pythonProject/Last/NA_data_fulllength2021"
    file_paths = os.listdir(dataset_path)
    import itertools
    import numpy as np
    n = len(file_paths)
    dist_matrix = np.zeros((n, n))
    # # using combinations instead of nested loops
    for SEQ1, SEQ2 in itertools.combinations(file_paths, 2):
        path1 = f"{dataset_path}/{SEQ1}"
        path2 = f"{dataset_path}/{SEQ2}"
        seq1 = read_files(path1)
        seq2 = read_files(path2)

        X = Global_Align(seq1, seq2)
        X.intialize_matrix()
        X.Fill_matrix()
        X.Trace_back()
        a, b, c = X.calculate_alignment()
        X.transition_transversion()
        distance = X.evolutionary_distance_k80()
        print(f"{SEQ1} vs {SEQ2}: {a}%\nEvolutioanry distnace:{distance}")
        i = file_paths.index(SEQ1)
        j = file_paths.index(SEQ2)
        dist_matrix[i, j] = distance
        dist_matrix[j, i] = distance

    return dist_matrix

# ----------------------------------------------------------UPGMA classes---------------------------------------------------------
class Node:
    def __init__(self, name,  left= None, right= None, height = 0):
        self.left = left
        self.right = right
        self.height = height
        self.name =name
    def __str__(self):
        if self.left is None and self.right is None:
            return self.name ##just a leaf or taxa
        return f"({self.left}, {self.right})" ##children

class UPGMA:
    def __init__(self, dist_matrix, labels):
        self.matrix = np.array(dist_matrix, dtype=float) ##to easy manipulate distance matrix

        self.clusters = []
        for name in labels:
            self.clusters.append(Node(name)) ##example [A2018, A2019, A2020]

        self.sizes = {}
        for s in range(len(labels)):
            self.sizes[s]=1

    def find_lowest_value(self):
        min_value = float("inf")
        x, y = -1, -1
        n = len(self.matrix)

        ##get the lowest value among the upper triangle of the distance matrix
        for i in range(n):
            for j in range(i + 1, n):
                if self.matrix[i][j] < min_value:
                    min_value = self.matrix[i][j]
                    x, y = i, j

        return x, y, min_value

##a, b are # integrs: cluster[0]-->node_a --> Node("2018")
    # aka it retrieves the Node object stored at index a in the list
    def merge_clusters(self, a, b, dist):
        node_a = self.clusters[a]
        node_b = self.clusters[b]

        new_height = dist / 2 ##ultrametric tree assumes that all leaves have the same height to the root
        new_node = Node(
            name=node_a.name + node_b.name,
            left=node_a,
            right=node_b,
            height=new_height
        ) ##(AB, left=A, right=B, height=dist/2)

        size_a = self.sizes[a]
        size_b = self.sizes[b]
        new_size = size_a + size_b ##no. of taxa in a cluster

        new_clusters = []
        new_sizes = {}
        index = 0

        for i in range(len(self.clusters)):
            if i != a and i != b: ##a, b are the old integrs
                new_clusters.append(self.clusters[i])
                new_sizes[index] = self.sizes[i]
                index += 1

        new_clusters.append(new_node)
        new_sizes[index] = new_size

        self.clusters = new_clusters
        self.sizes = new_sizes

    def update_matrix(self, a, b):
        n = len(self.matrix)
        size_a = self.sizes[a]
        size_b = self.sizes[b]
        new_size = size_a + size_b

        #indices we keep
        keep = []
        for i in range(n):
            if i != a and i != b:
                keep.append(i)

        #build new matrix with AB istead of A, B, initialize with zeros
        new_n = n-1
        new_matrix = np.zeros((new_n, new_n))

        #copy old distances
        for i_new in range(len(keep)):
            for j_new in range(len(keep)):
                i_old = keep[i_new] #i old which was 2=0
                j_old = keep[j_new]
                new_matrix[i_new][j_new] = self.matrix[i_old][j_old]

        #compute distances to new cluster
        k = 0
        for i_old in keep:
            dist_ka = self.matrix[i_old][a]
            dist_kb = self.matrix[i_old][b]
            new_dist = (dist_ka * size_a + dist_kb * size_b) / new_size

            new_matrix[k][new_n - 1] = new_dist
            new_matrix[new_n - 1][k] = new_dist
            k += 1

        self.matrix = new_matrix

    def build_tree(self):
        while len(self.clusters) > 1:
            a, b, dist = self.find_lowest_value()

            #update el.distance matrix first
            self.update_matrix(a, b)

            #merge tree nodes to clusters
            self.merge_clusters(a, b, dist)

        return self.clusters[0] ##root
    #============================================
    def to_tuple(self, node):
        if node.left is None and node.right is None:
            return node.name
        else:
            return (self.to_tuple(node.left), self.to_tuple(node.right))

# ------------------------------------class of MSA stats and bootstrap--------------------------------------------
class Statistics():
    def __init__(self, msafile):
        self.msaf = msafile
        self.sequences = self.readmsaf()
        self.aligcheck()

    def readmsaf(self):
        sequences = []
        with open(self.msaf) as file:
            seq = ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if seq:
                        sequences.append(seq)
                        seq = ''
                else:
                    seq += line
            if seq:
                sequences.append(seq)
        return sequences

    def aligcheck(self):
        seqlengths = set()
        for i in self.sequences:
            seqlengths.add(len(i))
        if len(seqlengths) != 1:
            print("not aligned")

    def colums(self, position):  # extract each column in the msa
        col = []
        for i in self.sequences:
            base = i[position]
            col.append(base)
        return col

    # alignment statistics
    def stat(self):
        length = len(self.sequences[0])
        taxa = len(self.sequences)
        A = 0
        T = 0
        C = 0
        G = 0
        for seq in self.sequences:
            for i in seq:
                if i == 'A':
                    A += 1
                elif i == 'T':
                    T += 1
                elif i == 'G':
                    G += 1
                elif i == 'C':
                    C += 1
        total_bases = A + T + G + C
        AT_content = ((A + T) / total_bases) * 100
        CG_content = ((C + G) / total_bases) * 100
        return {
            'length': length,
            'taxa': taxa,
            'AT_content': AT_content,
            'CG_content': CG_content
        }

    ##Conservation
    def shannon(self, column):
        nucleotidecounts = {}
        for i in column:
            if i != '-':
                if i in nucleotidecounts:
                    nucleotidecounts[i] += 1
                else:
                    nucleotidecounts[i] = 1
        total = sum(nucleotidecounts.values())
        H = 0.0
        for i in nucleotidecounts:
            p = nucleotidecounts[i] / total
            H -= p * math.log2(p)
        return H

    def conservation(self, entropy):
        Hm = 2
        C = 1 - (entropy / Hm)
        return C

    def allconservation(self):
        results = []
        for i in range(len(self.sequences[0])):
            column = self.colums(i)
            H = self.shannon(column)
            C = self.conservation(H)
            results.append(C)
        overall = sum(results) / len(results)
        return overall

    #bootstrap(bs)
    def bsalignment(self):
        numcol = len(self.sequences[0])
        newseq = []
        for i in self.sequences:
            newseq.append('')
        for j in range(numcol):
            randomcol = random.randint(0, numcol - 1)
            for seqidx in range(len(self.sequences)):
                nucleotide = self.sequences[seqidx][randomcol]
                newseq[seqidx] += nucleotide
        return newseq

    def splits(self, salmatree: tuple):
        splits = []
        if salmatree:
            left, right = salmatree
            if type(left) != tuple:
                left = (left,)
            if type(right) != tuple:
                right = (right,)
            splits.append((left, right))
            if type(left) == tuple and len(left) > 1:
                splits.extend(self.splits(left))
            if type(right) == tuple and len(right) > 1:
                splits.extend(self.splits(right))
        return splits

    def leaf(self, group):
        results = []
        for i in group:
            if type(i) == tuple:
                results.extend(self.leaf(i))
            else:
                results.append(i)
        return results

    def bstree(self, n=1000):
        branchs = {}
        for i in range(n):
            bssequences = self.bsalignment()
            splitslist = self.splits(tree)
            for b in splitslist:
                leaves = []
                for group in b:
                    leaves.extend(self.leaf(group))
                bset = frozenset(leaves)
                if bset in branchs:
                    branchs[bset] += 1
                else:
                    branchs[bset] = 1
        for b in branchs:
            branchs[b] = (branchs[b] / n) * 100
        return branchs

    #identify the mutation
    def translation(self, seq):
        protein = ''
        codons = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            protein += codons.get(codon, 'x')
        return protein

    def proteins(self):
        proteins = []
        for seq in self.sequences:
            proteins.append(self.translation(seq))
        return proteins

    def strain(self):
        name = []
        with open(self.msaf) as f:
            for line in f:
                if line.startswith('>'):
                    name.append(line.strip()[1:])
        return name

    def mutation(self):
        proteins = self.proteins()
        strains = self.strain()
        ref = proteins[0]
        strainm = {}
        for s in strains:
            if s == strains[0]:
                continue
            strainm[s] = []
        for pos in range(len(ref)):
            aaset = set()
            for i in range(1, len(proteins)):
                aaset.add(proteins[i][pos])
            if len(aaset) == 1 and list(aaset)[0] != ref[pos]:
                continue
            for i in range(1, len(proteins)):
                aa = proteins[i][pos]
                if aa != ref[pos]:
                    mut = ref[pos] + str(pos + 1) + aa
                    strainm[strains[i]].append(mut)
        return strainm

    def filteration(self):
        mutations = self.mutation()
        filtered = {}
        standardaa = "ACDEFGHIKLMNPQRSTVWY"
        for strain, mut in mutations.items():
            filtered[strain] = []
            for m in mut:
                ref = m[0]
                # pos = m[1:-1]
                newaa = m[-1]
                if newaa in standardaa:
                    if ref in standardaa:
                        filtered[strain].append(m)
        return filtered

    # visualization
    def tuple_to_newick(self, tree):
        if type(tree) == tuple:
            left, right = tree
            return f"({self.tuple_to_newick(left)},{self.tuple_to_newick(right)})"
        else:
            return str(tree)

    def visualization(self, tree):
        newick = self.tuple_to_newick(tree) + ";"
        tree = Phylo.read(StringIO(newick), "newick")
        bsdict = self.bstree()
        for clade in tree.find_clades():
            if not clade.is_terminal():
                leaves = []
                for leaf in clade.get_terminals():
                    leaves.append(leaf.name)
                leafset = frozenset(leaves)
                if leafset in bsdict:
                    clade.confidence = bsdict[leafset]
        Phylo.draw(tree)

    def report(self):
        print("report() started") ##check
        stats = self.stat()
        cons = self.allconservation()
        bs = self.bstree()
        with open('statistics.txt', 'w') as f:
            f.write('Alignment Statistics\n')
            f.write('...................................................\n')
            f.write('Sequence length: ' + str(stats['length']) + '\n')
            f.write('Number of taxa: ' + str(stats['taxa']) + '\n')
            f.write('AT content (%): ' + str(stats['AT_content']) + '\n')
            f.write('CG content (%): ' + str(stats['CG_content']) + '\n\n')
            f.write("Conservation Score\n")
            f.write('...................................................\n')
            f.write(str(cons) + '\n\n')
            f.write('Bootstrap Score for Each Branch\n')
            f.write('...................................................\n')
            for b in bs:
                f.write('Branch: ')
                for l in b:
                    f.write(l + ' ')
                f.write(':' + str(bs[b]) + '%' + '\n')


    def mutationsfile(self):
        mutations = self.filteration()
        with open('Mutation.txt', 'w') as f:
            for strain, mut in mutations.items():
                f.write('Strain: ' + strain + '\n')
                f.write('Mutations : \n')
                for m in sorted(set(mut)):
                    f.write(m + '\n')
                f.write('\n')



# ----------------------------------------driver code---------------------------------------------
# -------------------------------call_distance matrix function and save sequence names in a list -------------------------------------
distance_matrix = Construct_Distance_Matrix()
taxa_labels=[]
for i in os.listdir("D:/Junior/Programming-BMS 321/Project/pythonProject/Last/NA_data_fulllength2021"):
    i=i[0:-6]
    taxa_labels.append(i)
# -------------------------calling the UPGMA class and building the tree-----------------------------------------------------------------------------------------
upgma = UPGMA(distance_matrix, taxa_labels)
Node_obj=upgma.build_tree()
tree=upgma.to_tuple(Node_obj)
# print(tree)
# -------------------------calling the statistics class-----------------------------------------------------------------------------------------
stats=Statistics("D:/Junior/Programming-BMS 321/Project/pythonProject/Last/NAdata_MSAoutput2021.FASTA")
stats.report()
stats.mutationsfile()
print(stats.bsalignment())
stats.readmsaf()
stats.aligcheck()
stats.allconservation()
stats.visualization(tree)


from typing import final
from Bio import SeqIO
from Bio import AlignIO 
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

# NOTES: 
# max k-mer length is 100 (length of provided reads)
# hash table key is k-mer, value is all isoforms with that k-mer

# EXAMPLE TRANSCRIPT:
# >ENST00000529614
# GAGGCGGCTGTGGCCGCGTCGGTGTCCGCGTCGAGGAGCCGGGGCAGGGCACGATGGCGGACTGGGCTCGGGCTCAGAGC
# CCGGGCGCTGTGGAAGAGATTCTAGACCGGGAGAACAAGCGAATGGCTGACAGCCTGGCCTCCAAAGTCACCAGGCTCAA
# ATCGGACTCGGATTTCACAAGCATGACCAGCCTGCTTACAGGGAGCGTGAAGCGCTTTTCCACAATGGCAAGGTCCGGAC
# AAGACAACCGGAAGCTTCTATGTGGCATGGCCGTGGGTCTAATTGTGGCCTTCTTCATCCTCTCCTACTTCTTGTCCAGG
# GCAAGGACGTGAGCCAGTGGGAGCTGGTGTCTGTGGGTGCCAAGGGCAGCCAGGGTCTTCCCTGCCTGGTGTTTTGGGCT
# CCAGAGGACTTACCTACAAAATACTCCTTTGCAATGATAATTGTGGGTCAGGAATCTTCTTCCTGTGTGGCAGGAGGCTG
# CGGCTGCCTGTGACCTGATGAGCTCATGTTGGCTGGTCCCATGTGTGAAAGGGACTCTCTCGGGGAAACCAAGGCCCAGC

# EXAMPLE SET OF READS: 
# >read48/ENST00000410108;mate1:234-333;mate2:334-432
# TTTCAGGACAGCCAGGAGGGGGCGCACATCCGCCGAGAAACTGTGAGCAAGAGCGTCTGTGCTGAACCATGGCGCCACCAGAGGGCGCGCGATCCCGCCC
# >read49/ENST00000410108;mate1:89-188;mate2:189-287
# GCTCTTGCTCACAGTTTCTCGGCGGATGTGCGCCCCCTCCTGGCTGTCCTGAAATACCATGCCATCCAGGTACCGGTTCTGATCCTCTGCATCCCTATCG
# >read50/ENST00000410108;mate1:227-326;mate2:316-414
# CCAAATTAACACGACCCCCGTGCTGCCCTGAGGAAGTTGAAGCTCCTCGCTGCTTCTGGCACTTCAGCGGGAAGTTGGTTGGGGCGGGATCGCGCGCCCT

# EXAMPLE OUTPUT:
# counts | number of items in equivalence class | isoforms in equivalence class |
# 30     | 0                                    | NA                            |
# 5      | 1                                    | ENST000003679                 |
# 10379  | 2                                    | ENST000003679,ENST000009216   |

# ISSUES FACED:
# Adding in isoforms as a list so we can add multiple isoforms to one k-mer
# Counting the number of duplicates for each equivalence class to get the total counts
# Making sure to test the reverse strand

# k-mer length (max is 100)
k = input("Enter a k-mer length: ")
k = int(k)
if k > 100:
    print("Too long! The reads have a length of 100.")
    k = input("Enter a new k-mer length: ")
    k = int(k)
else:
    if k <= 0:
        print("Invalid input. Please enter a positive number.")
        k = input("Enter a new k-mer length: ")
        k = int(k)
print(k)

# Create hash table that we will use to store k-mers and isoforms
hash_table = {}
# hash_table = defaultdict(list)

# Use SeqIO from Biopython to read in the transcriptome and reads fasta files
# Nucleotide sequence is saved under .seq 
# Isoform is saved under .id
reads = list(SeqIO.parse("reads.fasta", "fasta"))
transcriptome = list(SeqIO.parse("transcriptome.fasta", "fasta"))

# Testing purposes
# SeqIO saves both the sequence and id (transcript/read name) for any item in a FASTA file
# print(transcriptome[0].seq)
# print(transcriptome[5].seq)
# print(transcriptome[5].id)

# FUNCTION: ADD_KMER
# Checks if a k-mer already exists in the hash table
# If no --> add the kmer and isoform as a new entry
# if yes --> check if the isoform is included in the list of isoforms for that k-mer
# if it is not --> add it to the list of isoforms for that k-mer 
# if it is --> does not do anything
def add_kmer(hashTable, kmer, isoform):
    if kmer not in hashTable:
        hashTable[kmer] = [isoform]
    else:
        if (isoform not in hashTable[kmer]):
            hashTable[kmer].append(isoform)
        
# Loop to add all the k-mers from every isoform in the transcriptome
# to the hash table using the add_kmer function from before 
# Modified from the Quiz 1 solutions code 
for t in transcriptome:
    for j in range(0, len(str(t.seq)) - k + 1):
        kmer=str(t.seq[j:j+k])
        # Testing purposes
        # print(kmer)
        add_kmer(hash_table, kmer, str(t.id))

# Testing purposes 
print(hash_table)

matches = []
equivalenceClass = None

# Loop through all of the reads 
for r in reads:
    # Loop through each read
    for j in range(0, len(r.seq) - k + 1):
        kmer=str(r.seq[j:j+k])
        if (equivalenceClass == "N/A"):
            exit
        if kmer not in hash_table:
            equivalenceClass = "N/A"
            exit
        # If we are on the first k-mer of a read set that k-mer's isoform as the current equivalence class
        if j == 0:
            equivalenceClass = hash_table.get(kmer)
        # Otherwise find the overlap between the current k-mer's isoforms and the current saved equivalence class
        else:
            isoforms = hash_table.get(kmer)
            # If no overlap possible --> exit
            if (equivalenceClass is None) or (isoforms is None): 
                equivalenceClass = "N/A"
                exit
            # Otherwise set equivalence class as overlap between previous EC and new isoforms
            else: 
                equivalenceClass = list(set(equivalenceClass) & set(isoforms))
    # If we found no equivalence class matches --> check reverse strand 
    if equivalenceClass == "N/A":
        for j in range(0, len(r.seq) - k + 1):
            original_strand = r.seq[j:j+k]
            reverse_strand = original_strand.reverse_complement()
            kmer=str(reverse_strand)
            # testing purposes
            # print(kmer)
            if (equivalenceClass == "N/A"):
                exit
            if kmer not in hash_table:
                equivalenceClass = "N/A"
                exit
            # If we are on the first k-mer of a read set that k-mer's isoform as the current equivalence class
            if j == 0:
                equivalenceClass = hash_table.get(kmer)
            # Otherwise find the overlap between the current k-mer's isoforms and the current saved equivalence class
            else:
                isoforms = hash_table.get(kmer)
                # If no overlap possible --> exit
                if (equivalenceClass is None) or (isoforms is None): 
                    equivalenceClass = "N/A"
                    exit
                # Otherwise set equivalence class as overlap between previous EC and new isoforms
                else: 
                    equivalenceClass = list(set(equivalenceClass) & set(isoforms))
    # Housekeeping
    if equivalenceClass is None:
        equivalenceClass = "N/A"
    # Housekeeping
    if len(equivalenceClass) == 0:
        equivalenceClass = "N/A"
    # Convert "N/A" from a string to a list so that the entries of the matches list are consistent
    # since the equivalence classes are saved as a list
    if equivalenceClass == "N/A":
        equivalenceClass = equivalenceClass.split()
    matches.append(equivalenceClass)

# Testing purposes 
print("MATCHES: ")
print(matches)

# Testing purposes
matches_array = np.array(matches)
print("MATCHES_ARRAY: ")
print(matches_array)

# Testing purposes
print("FINAL VALUES (ARRAY): ")
final_values, final_counts = np.unique(matches_array, return_counts=True)
print(final_values)
print("FINAL COUNTS (ARRAY): ")
print(final_counts)

# Count the total number of reads
sum = 0
for i in final_counts:
    sum = sum + i
print("TOTAL COUNTS: ")
print(sum) 

# Count the total number of matched reads
# N/A is saved as last entry so summing everything else gives us the total matches
sum = 0
for i in final_counts[:-1]:
    sum = sum + i
print("TOTAL MATCHES: ")
print(sum)

# For the number of isoforms in each equivalence class, we just count the number of 
# items in the equivalence class list
final_isoform_counts = []
for i in final_values:
    if i == 'N/A':
        n = 0
    else:
        n = len(i)
    final_isoform_counts.append(n)

# Convert into a numpy array and use pandas to build the final table
final_isoform_counts = np.array(final_isoform_counts)
matches_df = pd.DataFrame(data=[final_counts, final_isoform_counts, final_values]).T 
matches_df.columns = ["counts", "number of items in equivalence class" , "isoforms in equivalence class"]
print("MATCHES DATAFRAME: ") 
print(matches_df)

# matches_df.hist(column="counts")
# matches_df.hist(column="number of items in equivalence class")

# Save results to a csv file 
matches_df.to_csv(r'C:\Python310\results.csv', header=True, index=False)

# Testing purposes
# for record in reads:
#     print(record.seq + " " + record.id)

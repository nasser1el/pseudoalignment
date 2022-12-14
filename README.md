# Pseudoalignment

Given some RNA-seq data in FASTA format, the goal is to find the vector of equivalence class counts. The implementation does the following:
+ Given a gene annotation, an index should be produced which can be read in by your pseudoalignment procedure.
+ Your pseudoalignment procedure should take in:
  + (1) the previously mentioned index
  + (2) RNA-seq data and output equivalence class counts
  + (3) a k-mer length

Equivalence class counts are provided in the following format:

|Counts | Number of Items in Equivalence Class | Isoforms in Equivalence Class |
|-------|--------------------------------------|-------------------------------|
|30     | 0                                    | NA                            |
|5      | 1                                    | ENST000003679                 |
|10379  | 2                                    | ENST000003679,ENST000009216   |

<h1>Pseudoalignment Implementation</h1>

Pseudoalignment is a procedure used to find the equivalence classes for each set of reads generated from RNA-seq data using a gene annotation. For this project, the chromosome 11 transcriptome was provided as the transcriptome and a file of reads was also provided, both in FASTA format. In performing the pseudoalignment procedure, I implemented the naive version of the algorithm, using a hash table built from the transcriptome file. I then used the hash table to find the equivalence classes for each of the reads provided. In the end, about 82% of the reads matched to an equivalence class when using a k-mer length of 3. 

In implementing pseudoalignment, the first step was building the hash table, or index, from the provided gene annotation. To do this, I first utilized SeqIO from the Biopython library to read in the provided FASTA files. When parsing FASTA files, SeqIO is able to store both the id and sequence of a read. Using this feature, I was able to build a list of all of the sequences from the chromosome 11 transcriptome and the reads file as well as the corresponding ids. I then created a function that added k-mers to a hash table as keys with the corresponding list of isoforms as the values. This function first checked if the k-mer was not already in the hash table and, if so, added it along with its isoform. Otherwise, the function would check if the isoform of that k-mer was not already included in the list of isoforms. If it wasn???t, the isoform was appended to the list of isoforms. The hash table was built by iterating through every k-mer of length k of every transcript in the transcriptome and using the add_k-mer function, described above. 

![Flowchart of Logic!](/results/Flowchart%20of%20Logic.png)

To find the equivalence classes for each of the reads, I iterated through every k-mer of length k of every read in the file. For each of the k-mers, I first checked if the k-mer existed in the hash table as a key. If not, the equivalence class was set to N/A and we break out of the loop for that read. Next, I checked if the k-mer was the first k-mer of the read being mapped. If so, I set the current equivalence class, denoted by the variable equivalenceClass, as the corresponding isoforms for that k-mer. For all other k-mers in the read, I set a variable, isoforms, as the list of isoforms for that k-mer. I then updated the current equivalenceClass by finding the overlap between isoforms and the previously set variable, equivalenceClass. If at any point, the equivalenceClass is equal to ???N/A??? or None, we exit out of the loop since there???s no overlap between all of the isoforms found for all of the k-mers in the read. If we???ve traversed through the whole read and still found no equivalence class, I checked the reverse read. To do this, I again iterated through ever k-mer in the read, but this time using SeqIO???s reverse.complement() command to find the k-mer drawn from the opposite strand of DNA. Afterwards, I followed 
the same process outlined previously to find the equivalence class for the ???reverse??? k-mer. After going through an entire read, the final equivalenceClass result was appended to a list of matches. The list was then converted to a numpy array and the duplicate entries in the list, representing different reads with the same equivalence class, were counted and combined to produce the final counts and final list of matched equivalence classes with all of their isoforms. This data, along with the number of isoforms in each equivalence class, was exported as a .csv file and included in the submission for this project.

<h1>Issues in the Implementation of Pseudoalignment</h1>

In implementing pseudoalignment, I faced several problems including figuring out how to add multiple isoforms, as values, to a single k-mer (or key) in the hash table, how to count the total number of matched reads for each equivalence class, and testing the reverse strand for each read. Initially, the main problem when building the hash table was finding a proper way of storing multiple isoforms for one k-mer. I tried using append, but I ran into errors with invalid object types. The solution was to first store all of the transcripts in a list and make sure to convert the sequences of each transcript to a string since the source of the error was that SeqIO uses its own object type. The next major problem I faced was, after matching all of the reads to the hash table, outputting my results. I stored the final equivalence classes for each of the reads in a list, appending a new equivalence class each time we went through an entire read. However, this resulted in a list of lists (since each equivalence class is a list of isoforms) and this proved difficult to combine and count duplicate entries. To solve this, I converted the list to a numpy array and used numpy???s unique command to store all of the unique equivalence classes and their corresponding number of matches in two separate arrays which I then combined using pandas. However, after looking through my results, I had a low amount of matches, hovering around 40-50%. After talking to colleagues, I realized I had forgotten to account for the fact that the reads could be drawn from both strands of DNA. This led me to discover the reverse.complement() command from SeqIO which I then used to test the reverse strand for all unmatched reads. Finally, after fixing all of the above outlined errors, I was able to get sensible results.

<h1>Notable Statistics from Pseudoalignment Tests</h1>

Using the provided data on my pseudoalignment implementation as well as a k-mer length of 3, I found that 1,050,806 of the 1,282,526 total reads (~82%) ended up matching to one of the 10,180 unique equivalence classes. It should be noted that using a high k-mer length allowed the program to run more quickly compared to a lower k-mer length due to the lower number of k-mers that the program had to iterate through. However, I also adjusted the code to accept a valid input from the user for a k-mer length (defined as k=1 to 100), so it can be used on several different lengths of k. Some final interesting statistics were that the greatest number of isoforms in one equivalence class was 54, the average number of isoforms in an equivalence class was 5.02, the highest number of reads mapping to one equivalence class (excluding N/A) was 50,420, and the average number of reads mapping to each equivalence class was 126.

![Counts for Different Sizes of Equivalence Classes!](/results/Counts%20for%20Different%20Sizes%20of%20Equivalence%20Classes.png)
![Distribution of Counts!](/results/Distribution%20of%20Counts.png)
![Distribution of the Number of Isoforms in an Equivalence Class!](/results/Distribution%20of%20the%20Number%20of%20Isoforms%20in%20an%20Equivalence%20Class.png)

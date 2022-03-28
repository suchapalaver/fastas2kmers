# fastas2kmers

Written in Python, this k-mer counter parses the contents of an unzipped fasta file according to length k, which can be specified by the user. It outputs a table for each sequence of all kmers of length k, their reverse complement, and their total number in the sequence.

Download fastas2kmers.py to the directory with the fasta file you want to parse.     
From that directory run:

```
python fastas2kmers.py <fasta_file> <k>
```
e.g.: 
```
python fastas2kmers.py sample.fa 31
```

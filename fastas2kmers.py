#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 13:50:48 2021
@author: jojo
"""
import sys, time

def indexfasta(filename):
    infile = open(filename, 'rb') # opens and reads the fasta file in binary
    chunksize = 1024*1024 # it reads in these chunks
    filepos = 0
    headstart = list()
    headend = list()
    while True: # Exit loop when, chunk by chunk, we've gone through the whole file
        content = infile.read(chunksize)
        if len(content) == 0:
            break
        chunkpos = 0 # chunks away!
        while chunkpos != -1: # exit this loop when we're at the file's end
            chunkpos = content.find(b'>', chunkpos) # chunk from 1st identifier after current chunkpos
            if chunkpos != -1: # i.e. when there are no more left, find will return -1
                headstart.append(chunkpos + filepos) # headstart from '>' in record
                chunkpos += 1 # chunking beyond previous '>' to get to next record
        for i in range(len(headend), len(headstart)): # how many records we're looking for
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()
    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart)-1, 0, -1):
        if headend[i] == headend[i-1]:
            del headstart[i]
            del headend[i]            
    headstart.append(filepos)
    fastaindex = list()
    with open(filename, 'rb') as fh:
        for i in range(len(headend)):
            seq_start = headend[i]+1
            seq_end = headstart[i+1] - 1
            fh.seek(headstart[i])
            identifier = str(fh.read((headend[i]) - headstart[i])).strip("b'\n ")
            fastaindex.append((identifier, (seq_start, seq_end, seq_end-seq_start)))
    return fastaindex

def indexsequence(seq):
    pointer = 0
    seqindex = list()
    while len(seq) > pointer:
        potenstart = [seq.find(b'a', pointer), seq.find(b't', pointer), seq.find(b'c', pointer), seq.find(b'g', pointer)]
        realstart = min(potenstart)
        if realstart == -1:
            # happens rarely, so slow code is ok, apparently
            potenstart = [i for i in potenstart if i > -1]
            if len(potenstart) == 0:
                break
            realstart = min(potenstart)
        realend = seq.find(b'N', realstart)
        if realend == -1:
            realend = len(seq)
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex

def find_kmers(fasta, i):
    infile = open(fasta, 'rb')
    infile.seek(i[1][0])     
    identifier = i[0]
    seq = infile.read(i[1][1] - i[1][0]+1).translate(transtable, b'\r\n\t ')
    infile.close()
    # Index sequence
    seqindex = indexsequence(seq)
    subdict = dict()
    seqdict = dict()
    for start, stop in seqindex:
        for i in range(start, stop-kmer_len+1):
            kmer = str(seq[i:i+kmer_len]).strip("b'").upper()
            kmer = (kmer, kmer[::-1].translate(revcomptable))
            if kmer not in subdict:
                subdict[kmer] = 1
            else:
                subdict[kmer] += 1
        seqdict[identifier] = subdict
        yield seqdict

if __name__ == '__main__':
    
    file = str(sys.argv[1]) 
    kmer_len = int(sys.argv[2])
    
    transtable = bytes.maketrans(b'ATCGMRYKVHDBWmrykvhdbxnsw', b'atcgNNNNNNNNNNNNNNNNNNNNN')
    revcomptable = bytes.maketrans(b'ACGT', b'TGCA')
    
    filesprinted = []
    
    print()
    print()
    print()
    print(f"\tWorking file: {file}")
    print(f"\tSelected k-mer length: {kmer_len}")
    print()
    
    start = time.time()
    my_indexes = indexfasta(file)
    index_time = time.time()-start
    
    start = time.time()
    
    for i in my_indexes:
        for seq_dict in find_kmers(file, i):
            for key in seq_dict.keys():
                pfile = f"{key.strip('>').replace(' ', '_')}__{kmer_len}-mers.tsv"
                filesprinted.append(pfile)
                with open(f"{pfile}", "w+", newline='') as tsv:
                    print('K-mer\tR-C\tTotal', file=tsv)
                    # for kmer, freq in seq_dict[key].items(): # optional if sorting not required
                    for kmer, freq in sorted(seq_dict[key].items(), key=lambda fd: fd[1], reverse=True):
                        print('\t'.join([kmer[0], kmer[1], str(freq)]), file=tsv)
    search_time = time.time() - start

    print("")
    print("\tRESULTS FOR KMER FINDER:")
    print("\tIndexing time:", round(index_time,3))
    print("\tFinding kmers time:", round(search_time,3))
    print("\tTOTAL:",round(index_time + search_time,3))
    print()
    print(f"\tThe following {len(set(filesprinted))} files\n\thave been created:\n")
    for printedfile in sorted(set(filesprinted)):
        print(f"\t{printedfile}")
    print()
    sys.exit()
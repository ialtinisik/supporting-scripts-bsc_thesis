import peptides
import numpy as np
from copy import copy
import sys
import csv


class PeptideEnumerator():
    def __init__(self, max_len_per_chunk=6):
        self.max_len_per_chunk = max_len_per_chunk
        # self.peptide_len = peptide_len
        # self.aminoacids = ["A", "C", "D"]
        self.aminoacids = ["A", "C", "D", "E", "F", "G", "H",
                           "I", "K", "L", "M", "N", "P", "Q",
                           "R", "S", "T", "V", "W", "Y"]
        self.n_acids = len(self.aminoacids)

    def enumerate(self, peptide_len=6):
        ''' enumerates all possible sequences of peptide_len 
            by copying seqs generated for length i
            and appending self.aminoacids 
        '''
        seqs = self.aminoacids
        for _len in range(1, peptide_len):
            seqs = [
                s + aa for s in seqs for aa in self.aminoacids
            ]
        return seqs

    def __call__(self, peptide_len):
        ''' splits generation into chunks to keep memory in check
            max_len_per_chunk defines the cutoff, and shards sequences
            into subsets with prefix [A*, C*, ..., Y*] or [AA*, AC*, ...] '''
        max_len = self.max_len_per_chunk
        if peptide_len <= max_len:
            yield "", self.enumerate(peptide_len)
            # return self.enumerate(peptide_len)
        elif peptide_len >= 10:
            raise NotImplementedError
        else:
            extra_len = (peptide_len - max_len)
            print(f'peptide length {
                  peptide_len} greater than max len {max_len}')
            print(f'splitting into {self.n_acids**extra_len} chunks of {
                  self.n_acids**max_len}')

            if extra_len == 1:
                prefixes = self.aminoacids
            elif extra_len > 1:
                prefixes = self.enumerate(peptide_len=extra_len)

            seqs_hexa = self.enumerate(peptide_len=max_len)
            for prefix in prefixes:
                yield f'_{prefix}', [prefix + s for s in seqs_hexa]



if __name__ == "__main__":
    from src.design_heuristics import filter_synthesizable

    # for peptide_len in [2, 3, 4, 5]:
    total = 0
    pepgen = PeptideEnumerator(max_len_per_chunk=5)
    for peptide_len in [6]:
        pep_chunks = pepgen(peptide_len)
        for affix, seqs in pep_chunks:
            seqs = [filter_synthesizable(seq) for seq in seqs]


            count = sum([1 if s != None else 0 for s in seqs])
            total += count
            print(affix, 'count', count, '\t', 'frac', count/len(seqs))

            with open(f"./data/seqs/filtered_sequences_len{peptide_len}{affix}.txt", 
                      'w') as f:
                writer = csv.writer(f)
                for seq in seqs:
                    if seq is not None:
                        writer.writerow([seq])

        print(f'total count len{peptide_len}: {total}, frac {total/20**peptide_len}')

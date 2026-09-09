import gzip
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

n = 100_000
per_read = []   
grid     = []  

with gzip.open("data/ERR458495.fastq.gz", "rt") as handle:
    for i, record in enumerate(SeqIO.parse(handle, "fastq")):
        if i >= n:
            break
        q = record.letter_annotations["phred_quality"]
        per_read.append(sum(q) / len(q))
        grid.append(q)

grid = np.array(grid)              
per_position = grid.mean(axis=0)    
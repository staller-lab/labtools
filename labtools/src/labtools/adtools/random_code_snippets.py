import time
import os
from adtools.finder import pull_AD
from adtools.seqlib import read_fastq

def race(file):
    """ Two methods of counting reads compete."""

    print(os.path.getsize(file))
    tic = time.perf_counter()
    seq_counter(file)
    toc = time.perf_counter()
    print(toc-tic)
    tic = time.perf_counter()
    seq_counterX(file)
    toc = time.perf_counter()
    print(toc-tic)

# deprecated 9/27/2022 
# keeping it as a snippet because it is useable in downsampling
# until I add in the subset method to the new fastq method
def seq_counter_old(fastq, design_to_use = None, barcoded = False, subset = None, **kwargs):
    """ Count occurences of ADs or AD/bc pairs."""

    seqCounts = {}
    for line in read_fastq(fastq, subset):
        AD,bc = pull_AD(line[1], barcoded, **kwargs)
        
        if barcoded and AD != None:
            AD = (AD, bc)
        if AD not in seqCounts and AD != None: seqCounts[AD] = 1
        elif AD != None: seqCounts[AD] += 1
    counts = pd.Series(seqCounts)

    if design_to_use:
        design = pd.read_csv(design_to_use)
        if barcoded:
            counts = counts.where(counts.index.droplevel(1).isin(design.ArrayDNA)).dropna()
        else:
            counts = counts.where(counts.index.isin(design.ArrayDNA)).dropna()
    return counts
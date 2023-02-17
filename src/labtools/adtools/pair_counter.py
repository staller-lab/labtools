from adtools.finder import pull_AD
from adtools.seqlib import read_fastq, read_fastq_big
import pandas as pd

def seq_counter(fastq, design_to_use = None, barcoded = False, **kwargs):
    """Count occurences of ADs or AD/bc pairs."""

    seqCounts = {}
    for line in read_fastq_big(fastq):
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

def sort_normalizer(pair_counts, bin_counts):
    """Normalize by reads per sample, reads per tile and reads per bin."""
    
    df = pd.DataFrame(pair_counts)
    df.fillna(0, inplace=True)
    # 10 is the read minimum, should make this user defined
    df = df.loc[:, (df > 10).any(axis=0)]
    df = df.transpose()
    numreads = df.sum(axis = 1)
    reads = df.copy(deep = True)
    #df = df_in.copy(deep=True)
    # row i column j
    for j in df:
        df[j] = (df[j]/df[j].sum())/bin_counts[j]
    for i, pair in enumerate(df.index):
        df.iloc[i] = df.iloc[i]/df.iloc[i].sum()
    
    return df, numreads, reads

def calculate_activity(df_in, bin_values, min_max = False):
    df = df_in.copy(deep=True)
    activities = df_in.dot(bin_values)
    if min_max:
        activities = minmax_scale(activities)
    
    df.loc[:,"Activity"] = activities
    return df

def main():
    pass

if __name__ == '__main__':
    main
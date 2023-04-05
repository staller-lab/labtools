from labtools.adtools.finder import pull_AD, pull_barcode
from labtools.adtools.seqlib import read_fastq, read_fastq_big
import pandas as pd

def seq_counter(fastq, design_to_use = None, barcoded = False, only_bcs = False, **kwargs):
    """Counts occurences of ADs or AD-barcode pairs in a fastq file.
    
    Parameters
    ----------
    fastq : str 
        Path to fastq or fastq.gz file.
    design_to_use : str, default None
        Path to csv file containing ArrayDNA column.
    barcoded : bool, default False
        Whether to count ADs with different barcodes separately.
    only_bcs : bool, default False
        Whether or not to only look for barcodes in the sequence.
    
    Returns
    ----------
    counts : pandas.core.series.Series
        Pandas series where indices are AD or AD/barcode sequences and values are counts.
    
    Examples
    ----------
    >>> seq_counter("../exampledata/mini.fastq")
    GGTTCTTCTAAATTGAGATGTGATAATAATGCTGCTGCTCATGTTAAATTGGATTCATTTCCAGCTGGTGTTAGATTTGATACATCTGATGAAGAATTGTTGGAACATTTGGCTGCTAAA    1
    GAAGAATTGTTTTTACATTTGTCTGCTAAGATTGGTAGATCTTCTAGGAAACCACATCCATTCTTGGATGAATTTATTCATACTTTGGTTGAAGAAGATGGTATTTGTAGAACTCATCCA    3
    dtype: int64
    """

    if type(fastq) == list:
        seqCounts = {}
        for file in fastq:
            for line in read_fastq_big(file, **kwargs):
                AD,bc = pull_AD(line[1], barcoded, **kwargs)
                
                if barcoded and AD != None:
                    AD = (AD, bc)
                if AD not in seqCounts and AD != None: seqCounts[AD] = 1
                elif AD != None: seqCounts[AD] += 1
        counts = pd.Series(seqCounts)
    elif only_bcs == True and design_to_use == None:
        seqCounts = {}
        for line in read_fastq_big(fastq, **kwargs):
            bc = pull_barcode(line[1],**kwargs)
            
            if bc not in seqCounts and bc != None: seqCounts[bc] = 1
            elif bc != None: seqCounts[bc] += 1
        counts = pd.Series(seqCounts)
    elif only_bcs == True and design_to_use != None:
        raise(TypeError, "cannot use a design file with only barcodes")
    else:
        seqCounts = {}
        for line in read_fastq_big(fastq, **kwargs):
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

def create_map(ad_bcs):
    """Converts output of seq_counter with AD,bc pairs to a dict  map.
    
    If the barcode is found with two different ADs, it is not included in 
    the dictionary.

    Parameters
    ----------
    ad_bcs : pd.Series 
        output counts from seq_counter with barcoded = True.

    Returns
    ----------
    bc_dict : dict
        Dictionary with barcodes as keys and 1 AD as value.
    """
    bc_dict = {}
    for line in zip(ad_bcs.index, ad_bcs):
        ad = line[0][0]
        bc = line[0][1]
        count = line[1]
        if bc not in bc_dict:
            bc_dict[bc] = ad
        elif bc in bc_dict and bc_dict[bc] == ad:
            pass
        elif bc in bc_dict and bc_dict[bc] != ad:
            pass
            #bc_dict[bc] = "multiple_matches"
    return(bc_dict)

def convert_bcs_from_map(bcs, bc_dict):
    """Takes bc only data and uses a barcode dictionary to return AD counts.
    
    If the barcode is found with two different ADs, it is not included in 
    the dictionary.

    Parameters
    ----------
    bcs : pd.Series 
        output counts from seq_counter with only_bcs = True.
    bc_dict : dict
        Dictionary with barcodes as keys and 1 AD as value from create_map().

    Returns
    ----------
    converted : pd.Series
        Pandas series where indices are AD sequences and values are counts.
    """
    ads = []
    for bc in bcs.index:
        ad = None
        if bc in bc_dict:
            ad = bc_dict[bc]
        ads.append(ad)
    ad_col = pd.Series(ads, name = "AD")
    x = pd.DataFrame(bcs).reset_index()
    x["AD"] = ads
    convereted = x[[0, "AD"]].groupby("AD").sum()[0]
    return convereted


def sort_normalizer(pair_counts, bin_counts):
    """Normalize by reads per sample, reads per tile and reads per bin.
    
    Parameters
    ----------
    pair_counts : list of pandas.core.series.Series 
        List of pandas series where indices are AD or AD/barcode sequences and values are counts.
    bin_counts : list
        List of number of cells per bin in the same order as the pair counts.
    
    Returns
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the normalized read counts.
    numreads : 
        unknown
    reads :
        unknown
    
    Examples
    ----------
    >>> sort_normalizer([count1, count2], [1000,1000])
    """
    df = pd.DataFrame(pair_counts)
    df.fillna(0, inplace=True)
    # 10 is the read minimum, should make this user defined
    df = df.loc[:, (df.sum() > 10)]
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
    """Calculate the activity of a normalized sort df.
    
    Parameters
    ----------
    df_in : pandas.DataFrame
        Dataframe output of sort_normalizer()
    bin_values : list
        List of mean or median fluorescence values per bin in the same order as the pair counts.
    min_max : bool, default False
        Whether to normalize the activity using min 0 max 1.
    
    Returns
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the activity values per sequence or sequence-barcode pair.
    """
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
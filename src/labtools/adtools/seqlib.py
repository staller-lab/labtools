import gzip
from random import sample
from sklearn.utils.random import sample_without_replacement

# fasta reader
def read_fasta(filename):
    """Generator for reading entries in a fasta file.
    
    Yields 2 lines of a fasta file at a time (name, seq).
    
    Parameters
    ----------
    filename : str 
        Path to fasta or fasta.gz file.
    
    Yields
    ----------
    (name, seq) : (str, str)
        Name of sequence, biological sequence.
    
    Examples
    ----------
    >>> for line in read_fasta("example.fasta"):
    ...     name = line[0]
    ...     seq = line[1]
    ...     print(name, seq)
    Geraldine
    ACGTGCTGAGGCTGCGCTAGCAT
    Gustavo
    CTGATGCTAGATGCTGATA
    """
    name = None
    seqs = []

    fp = None 
    if filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
    else: fp = open(filename)

    for line in fp.readlines():
        line = line.rstrip()
        if line.startswith('>'):
            if len(seqs) > 0:
                seq = ''.join(seqs)
                yield(name, seq)
                name = line[1:]
                seqs = []
            else:
                name = line[1]
        else:
            seqs.append(line)
    yield(name, ''.join(seqs))
    fp.close()

# fastq reader
# this method opens the entire fastq into memory
# which is great... as long as the fastq file isn't > your RAM lol
# yes this happened so please see the read_fastq_big method below for that
def read_fastq(filename, subset=None):
    """Generator for reading entries in a fastq file.
    
    Yields 4 lines of a fastq file at a time (name, seq, +, error).
    
    Parameters
    ----------
    filename : str 
        Path to fastq or fastq.gz file.
    subset : int, optional
        Number of reads to randomly sample from the fastq file.
    
    Yields
    ----------
    (name, seq, qual) : (str, str, str)
        tuple of str containing name, seq and quality for entry.
    
    Examples
    ----------
    >>> for line in read_fastq("example.fasta"):
    ...     name = line[0]
    ...     seq = line[1]
    ...     qual = line[2]
    ...     print(name, seq)
    Geraldine
    ACGTGCTGAGGCTGCGCTAGCAT
    Gustavo
    CTGATGCTAGATGCTGATA
    """
    # add a warning for large fastq files and suggest using read_fastq_big
    name = None
    seqs = []

    fp = None 
    if filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
    else: fp = open(filename)
    lines = fp.readlines()

    all_reads = range(0, len(lines), 4)
    if subset:
        subset_indices = sample_without_replacement(len(all_reads), subset)
        subset_reads = [all_reads[i] for i in subset_indices]
    else: subset_reads = all_reads

    for num in subset_reads: 
        name = lines[num]
        seq = lines[num+1]
        opt = lines[num+2]
        qual = lines[num+3]
        yield(name.rstrip(), seq.rstrip(), qual.rstrip())
    fp.close()

def read_fastq_big(filename):
    """Generator for fastq file without opening into memory.
    
    Yields 4 lines of a fastq file at a time (name, seq, +, error).
    Useful in situations where the fastq file is large and opening into RAM
    would crash computer. Does not currently support subsetting.
    
    Parameters
    ----------
    filename : str 
        Path to fastq or fastq.gz file.
    
    Yields
    ----------
    (name, seq, qual) : (str, str, str)
        tuple of str containing name, seq and quality for entry.
    
    Examples
    ----------
    >>> for line in read_fastq_big("example.fasta"):
    ...     name = line[0]
    ...     seq = line[1]
    ...     qual = line[2]
    ...     print(name, seq)
    Geraldine
    ACGTGCTGAGGCTGCGCTAGCAT
    Gustavo
    CTGATGCTAGATGCTGATA
    """
    name = None
    seqs = []

    if filename.endswith('.gz'): 
        opener = gzip.open(filename, 'rt')
    else: opener = open(filename)

    with opener as file:
        linenum = 0
        for line in file:
            linenum += 1
            if linenum == 5:
                linenum = 1
            if linenum == 1: name = line
            elif linenum == 2: seq = line
            elif linenum == 3: opt = line
            elif linenum == 4: 
                qual = line
                yield(name.rstrip(), seq.rstrip(), qual.rstrip())
    opener.close()

def get_numreads(filename):
    """Returns number of reads in a fastq file.
    
    Parameters
    ----------
    filename : str 
        Path to fastq or fastq.gz file.
    
    Returns
    ----------
    numreads : int
        Number of reads in the fastq file.
    
    Examples
    ----------
    >>> get_numreads("example.fastq")
    124
    """

    fp = None
    if filename.endswith(".gz"): fp = gzip.open(filename, 'rt')
    else: fp = open(filename)
    lines = fp.readlines()

    numreads = len(range(0, len(lines), 4))

    return numreads
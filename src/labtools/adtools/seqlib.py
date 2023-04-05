import gzip
from random import sample
import subprocess
from sklearn.utils.random import sample_without_replacement
import csv

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
                name = line[1:]
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

def read_fastq_big(filename, subset = None, **kwargs):
    """Generator for fastq file without opening into memory.
    
    Yields 4 lines of a fastq file at a time (name, seq, +, error).
    Useful in situations where the fastq file is large and opening into RAM
    would crash computer. Supports subsetting with sklearn.sample_without_replacement().
    
    Parameters
    ----------
    filename : str 
        Path to fastq or fastq.gz file.
    subset: int
        Number of reads to randomly subsample from file.
    
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

    numreads = get_numreads(filename)
    # every 4 lines is a read in fastq
    all_reads = range(0, numreads*4, 4)
    if subset:
        # a list of indices to subsample
        subset_indices = sample_without_replacement(len(all_reads), subset)
        # get the actual read number from the index (aka the index*4)
        subset_reads = [all_reads[i] for i in subset_indices]
    else: subset_reads = all_reads
        
    with opener as file:
        linenum = 0
        readnum = 0
        for line in file:
            linenum += 1
            if linenum == 5:
                linenum = 1
            if linenum == 1: name = line
            elif linenum == 2: seq = line
            elif linenum == 3: opt = line
            elif linenum == 4: 
                qual = line
                if (readnum+1) in subset_reads:
                    yield(name.rstrip(), seq.rstrip(), qual.rstrip())
            readnum += 1
    opener.close()

def get_numreads(filename):
    """Returns number of reads in a fastq or fastq.gz file.
    
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
    if filename.endswith(".gz"): 
        bashCommand = f"zgrep -c '@' {filename}"
    else:
        bashCommand = f"grep -c '@' {filename}"
    
    process = subprocess.run(bashCommand, shell=True, capture_output=True, text=True)
    numreads = int(process.stdout)
        
    return numreads

def get_numreads_old(filename):
    """Returns number of reads in a fastq file.
    
    Parameters
    ----------
    filename : str 
        Path to fastq file.
    
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

def write_bc_dict(bc_dict, name):
    """Writes bc_dict to a csv.
        
    Parameters
    ----------
    bc_dict : dict 
        Dictionary output from counter.create_map().
    name : str
        Filename for output csv. Ex "Library1_dictionary"
    """
    with open(f'{name}', 'w') as f:
        w = csv.DictWriter(f, bc_dict.keys())
        w.writeheader()
        w.writerow(bc_dict)

def read_bc_dict(filename):
    """Reads bc_dict from a csv.
        
    Parameters
    ----------
    filename : str
        Path to csv containing a single dictionary.

    Returns
    ----------
    bc_dict : dict
        Dictionary.     
    """
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for d in reader:
            bc_dict = d
    return(bc_dict)
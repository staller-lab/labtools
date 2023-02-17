import numpy as np
import random
import pandas as pd
import math
import more_itertools
from collections import Counter

def windowed_shuffle(seq, window_size=5, step_size=1, n_shuffles=1, short_last_window = False):
    """Return a number of uniquely shuffled sequences.
    
    If n_shuffles is higher than the number of possible unique shuffles in that window, 
    only the unique shuffles are returned. Returns sequences with shuffled letters in lowercase.
    
    Parameters
    ----------
    seq : str
        Sequence to shuffle.
    window_size : int, default 5
        Size of shuffle region.
    step_size : int, default 1
        Distance between window starts (default 1).
    n_shuffles : int, default 1
        Number of shuffles to make for each region (default 1).
    short_last_window : bool, default False
        Whether last window can be less than window_size to include all letters (default False).
    
    Returns
    ----------
    list
        A list of shuffled sequences.
    
    Examples
    ----------
    >>> shuffles_list, names_list = shuffle.windowed_shuffle("ABCDEFGHIJKLMNOPQRSTUVWYXZ")
    ...shuffles_list
    ['decbaFGHIJKLMNOPQRSTUVWYXZ',
     'AfdecbGHIJKLMNOPQRSTUVWYXZ',
     'ABfcdegHIJKLMNOPQRSTUVWYXZ',
     'ABCfdhgeIJKLMNOPQRSTUVWYXZ',
     'ABCDhgeifJKLMNOPQRSTUVWYXZ',
     'ABCDEifgjhKLMNOPQRSTUVWYXZ',
     'ABCDEFhjgkiLMNOPQRSTUVWYXZ',
     'ABCDEFGkjilhMNOPQRSTUVWYXZ',
     'ABCDEFGHikmljNOPQRSTUVWYXZ',
     'ABCDEFGHIjmlknOPQRSTUVWYXZ',
     'ABCDEFGHIJomklnPQRSTUVWYXZ',
     'ABCDEFGHIJKpomlnQRSTUVWYXZ',
     'ABCDEFGHIJKLpoqnmRSTUVWYXZ',
     'ABCDEFGHIJKLMronqpSTUVWYXZ',
     'ABCDEFGHIJKLMNsqrpoTUVWYXZ',
     'ABCDEFGHIJKLMNOrsqtpUVWYXZ',
     'ABCDEFGHIJKLMNOPsutrqVWYXZ',
     'ABCDEFGHIJKLMNOPQstruvWYXZ',
     'ABCDEFGHIJKLMNOPQRvsutwYXZ',
     'ABCDEFGHIJKLMNOPQRSuvwytXZ',
     'ABCDEFGHIJKLMNOPQRSTxuvywZ',
     'ABCDEFGHIJKLMNOPQRSTUwxzyv']
    """
    sequence_list = []
    name_list = []
    step_size = int(step_size)
    
    # determine the number of windows
    if short_last_window:
        rounder = math.ceil
    else: 
        rounder = math.floor
    if window_size == step_size:
        n_windows = rounder(len(seq)/step_size)
    else:
        n_windows = rounder((len(seq)-window_size+step_size)/step_size)
    
    for window in range(0,n_windows):
        # determine first and last indices of the shuffle window
        start = 0+step_size*window
        end = start+window_size
        # divide the sequence into regions to shuffle or not shuffle
        pre_shuffle = seq[0:start]
        shuffle_region = seq[start:end]
        post_shuffle = seq[end:]
        
        # calculate the maximum number of unique shuffles for the given sequence
        # see explanation at https://getcalc.com/statistics-permutation-calculator.htm
        # this prevents the later while loop from looking for new shuffles when none exist
        max_unique_shuffles = math.factorial(len(shuffle_region)) \
                            /(np.prod([math.factorial(num) for num in list(Counter(shuffle_region).values())]))-1
        n = n_shuffles
        if max_unique_shuffles < n_shuffles:
            n = max_unique_shuffles
            
        # create shuffles
        shuffle_list = []
        for j in np.arange(n):
            shuffle = ''.join(more_itertools.random_permutation(shuffle_region))
            # keep generating shuffles until a new one is obtained
            while shuffle in shuffle_list or shuffle == shuffle_region:
                shuffle = ''.join(more_itertools.random_permutation(shuffle_region))
            # keep track of generated shuffles
            shuffle_list.append(shuffle)
            # append the full sequences and names
            shuffled_seq = pre_shuffle + shuffle.lower() + post_shuffle
            sequence_list.append(shuffled_seq)
            name = f'shuffle_size{window_size}_ss{step_size}_pos{window}_{int(j+1)}'
            name_list.append(name)
    return sequence_list, name_list
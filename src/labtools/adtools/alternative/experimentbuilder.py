import sys
from typing import Sequence, Type
from ADtools.seqlib import read_fastq
from ADtools.finder import pull_AD

# this was an exercise in OOP but it is not really useful (yet)

# define a sequence parent class
class _Sequence:
    """Parent class of AD and barcode."""

    def __init__(self, seq):
        self.seq = seq
        self.count = 1
    
    def update(self):
        """Report that the sequence was found again."""

        self.count += 1

    def __eq__(self, other):
        """Custom operator to compare equality of sequence objects."""

        return (self.seq == other.seq and self.count == other.count)

    def __str__(self):
        """When printing the object, show the sequence."""
        #print(self.count)
        return f"{(self.seq, self.count)}"

    def __lt__(self, other):
        if type(other) == type(self):
            return self.seq < other.seq
        else: 
            raise TypeError(f"'<' not supported in instance of "\
            "{type(self)} and {type(other)}")

    def __gt__(self, other):
        if type(other) == type(self):
            return self.seq > other.seq
        else: 
            raise TypeError(f"'>' not supported in instance of "\
            "{type(self)} and {type(other)}")

    def has_same_seq(self, other):
        if type(other) == type(self):
            return self.seq == other.seq
        else: raise TypeError("only works with another _Sequence")

    def is_seq_in_list(self, _list):
        if self.seq in [barcode.seq for barcode in _list]:
            return True
        else: return False

    def merge(self, other):
        if self.has_same_seq(other):
            print("I am going to merge")
        else: print("Not merging")

class SeqList(list):
    """List of sequence objects."""

    def has_same_sequences(self, other):
        if type(other) == SeqList:
            return self.sort() == other.sort()
        else: raise TypeError("only works with another SeqList")
    
    def __str__(self):
        return f"{[(item.seq, item.count) for item in self]}"

# define a barcode
class Barcode(_Sequence):
    """A short DNA sequence used to identify a read."""

    def __init__(self, seq):
        super().__init__(seq)
        self.ADs = SeqList()

    def link_AD(self, AD):
        """Add the barcode object to the list of associated barcodes."""

        if not AD.is_seq_in_list(self.ADs):
            self.ADs.append(AD)
        else:
            current_AD = next(ad for ad in self.ADs if ad.seq == AD.seq)
            current_AD.update()

# define an AD
class AD(_Sequence):
    """A designed DNA tile attached to a transcription factor."""
    
    def __init__(self, seq):
        super().__init__(seq)
        self.AA_seq = None # implement
        self.barcodes = SeqList()

    def link_barcode(self, barcode):
        """Add the barcode object to the list of associated barcodes."""

        if not barcode.is_seq_in_list(self.barcodes):
            self.barcodes.append(barcode)
        else:
            current_barcode = next(bc for bc in self.barcodes if bc.seq == barcode.seq)
            current_barcode.update()

# pure chaos I have no idea what is happening
class Sample:
    """An individual sequencing sample (bin)."""

    def __init__(self, file, name, position, barcoded = True, event_count = 100000):
        self.file = file
        self.name = name
        self.position = position
        self.ADs = SeqList()
        self.barcoded = barcoded
        self.event_count = event_count

    def parse(self):
        """Count tiles and assign barcodes in the sample."""
        
        none_code = AD("uninterpretable")
        self.ADs.append(none_code)

        for line in read_fastq(self.file):
            read = line[1]
            tile, barcode = pull_AD(read, self.barcoded)
            if tile == None:
                none_code.update()
                continue
            if not any(tile for AD in self.ADs if AD.seq == tile):
                current_AD = AD(tile)
                self.ADs.append(current_AD)
            else:
                current_AD = next(AD for AD in self.ADs if AD.seq == tile)
                current_AD.update()
            current_barcode = Barcode(barcode)
            current_barcode.link_AD(current_AD)
            current_AD.link_barcode(Barcode(barcode))           


class Experiment:
    """A group of samples from the same sort"""

    def __init__(self, name, num_bins):
        self.name = name
        self.num_bins = num_bins
        self.samples = []
        self.event_counts = {}

    def add_sample(self, sample):
        """Add a sample to the experiment."""

        if type(sample) == Sample:
            self.samples.append(sample)
        else: raise TypeError("Sample must be a Sample object.")

        self.event_counts[sample.name] = sample.event_count
    
    def normalize(self, data_in, bincounts = None):
        """Normalize the read data."""

        if len(self.samples) != self.num_bins:
            raise ValueError(f"Not enough samples. Expected {self.num_bins}.")
            
        data=data_in.copy(deep=True)
        for col in data.columns:
            data[col] = data[col]/data[col].sum()
        if bincounts:
            count = 0
            for bincount in bincounts:
                data[data.columns[count]] = data[data.columns[count]]/bincount
                count += 1
        for row in data.index:
            data.iloc[row] = data.iloc[row]/data.iloc[row].sum()
        return data



    
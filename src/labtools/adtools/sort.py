from labtools.adtools.pair_counter import *
from labtools.adtools.seqlib import get_numreads
from labtools.adtools.pair_counter import seq_counter
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

class Sort():
    """A sort object.
    
    Attributes
    ----------
    data_files : list
        List of str paths to fastq files in order of bins.
    bin_counts : list
        List of number of cells per bin in order of data files.
    bin_values : list
        List of mean or median fluorescence values per bin in order of data files.
    design_file : str
        Path to design csv file containing "ArrayDNA" header with DNA sequences to search
    """

    def __init__(self, data_files, bin_counts, bin_values, design_file = None):
        self.data_files = data_files
        self.bin_counts = bin_counts
        self.bin_values = bin_values
        self.design_file = design_file

    def process(self, csv=False, **kwargs):
        """Calculate the activity for each tile.
        
        Parameters
        ----------
        csv : bool, default False
            Use a csv input for predetermined AD counts (for example, bowtie output).
        **kwargs : dict, optional
            Extra arguments to pull_AD (for example changing the anchor sequences).
            
        Returns
        ----------
        processed_sort : pandas.DataFrame
        numreads : 
        reads :
       
        Examples
        ----------
        >>>Sort.process()
        """

        sort_list = []
        
        # super weird and needs to be double checked
        # converts output from bowtie bam to csv into same format
        # as the seq_counter() fxn output for downstream processing
        if csv == True: 
            for sample in self.data_files:
                parsed_sample = pd.read_csv(sample, index_col=1, header = None).squeeze('columns')
                sort_list.append(parsed_sample)
        else: 
            for sample in self.data_files:
                parsed_sample = seq_counter(sample, design_to_use = self.design_file, **kwargs)
                sort_list.append(parsed_sample)
        
        normed_sort, numreads, reads = sort_normalizer(sort_list, self.bin_counts)

        processed_sort = calculate_activity(normed_sort, self.bin_values)

        return processed_sort, numreads, reads

    def downsample(self, subset_size):
        """Perform downsampling on raw reads, then analyze."""

        sort_list = []
        for sample in self.data_files:
            numreads = get_numreads(sample)
            subsample_list = []
            
            for num in range(0, numreads, subset_size):
                subsample = seq_counter(sample, design_to_use = self.design_file, subset = num)
                subsample_list.append(subsample)
            sort_list.append(subsample_list)
        
        all_subsamples = []

        min_reads = min([len(sort) for sort in sort_list])
        for subsample in range(0, min_reads):
            group = [sort_list[i][subsample] for i in range(0, len(sort_list))]
            df, _ = sort_normalizer(group, self.bin_counts)
            final = calculate_activity(df, self.bin_values)
            all_subsamples.append(final)

        for i in range(1, len(all_subsamples)):
            all_subsamples[0][f"{i}"] = all_subsamples[i].Activity

        downsampling_df = all_subsamples[0].iloc[:,len(sort_list):]

        RMSE_list = []
        for j in downsampling_df:
            if j == "Activity":
                RMSE = downsampling_df.Activity
            else:
                RMSE = ((downsampling_df.Activity - downsampling_df[j]) ** 2 ) ** 0.5
            RMSE.rename(j, inplace = True)
            RMSE_list.append(RMSE)

        errors = pd.DataFrame(RMSE_list)
        plot_error = errors.transpose().iloc[:,1:]

        fig, axes = plt.subplots(1,1, figsize=(18, 10), sharex = True)
        axes.set_title("Subsampling RMSE", size=14)
        #axes.set_ylim(0,75)
        axes.set_ylabel("RMSE")
        axes.set_xlabel(f"Number of {subset_size} reads")
        sns.boxplot(data = plot_error)

        return errors.transpose()
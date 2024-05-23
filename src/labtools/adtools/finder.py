import re

def pull_AD(read, barcoded = False, ad_preceder = "GCTAGC", 
bc_preceder = "GGGCCCG", bc_anteceder = "GGAGAGAA", ad_length = 120, bclength = 11, loss_table=None, **kwargs):
    """Find the activation domain tile in a read.
    
    Takes a read sequence and uses customizable anchor sequences to locate a 
    variable sequence (AD/seq of interest) in the read. Includes support for barcodes.
    
    Parameters
    ----------
    read : str 
        The biological read of interest.
    barcoded : bool, default False
        Whether or not the sequence includes a barcode in addition to the AD/seq of interest.
    ad_preceder : str, default "GCTAGC"
        The anchor sequence directly before the AD.
    bc_preceder : str, default "GGGCCCG"
        The anchor sequence directly before the barcode.
    bc_anteceder : str, default "GGAGAGAA"
        The anchor sequence directly after the barcode.
    ad_length : int, default 120
        The length of the AD/seq of interest.
    bc_length : int, default 11
        The length of the barcode sequence if used.
    
    Returns
    ----------
    AD : str
        The sequence of interest, if located. Else None.
    barcode : str
        The barcode, if used and located. Else None.
    
    Examples
    ----------
    >>> pull_AD("ACTTTTATVGCTAGCATGGCTGGTAGATCTTGGTTGATTGATTCTAATAGAATTGCTACTAAGATTATGTCTGCTTCTGCTTCTTCTGATCCAAGACAAGTTGTTTGGAAATCTAATCCATCTAGACATTGTCCAGCTGATCGATGCTAGTAGAGAGAGA")
    ATGGCTGGTAGATCTTGGTTGATTGATTCTAATAGAATTGCTACTAAGATTATGTCTGCTTCTGCTTCTTCTGATCCAAGACAAGTTGTTTGGAAATCTAATCCATCTAGACATTGTCCA
    """
    searched_read = re.split(ad_preceder, read, maxsplit=1)
    AD = None
    barcode = None

    if len(searched_read) == 2:
        roi = searched_read[1]
        if barcoded:
            searched_read = re.split(bc_preceder, roi[ad_length:], maxsplit=1)
            if len(searched_read) == 2:
                barcode = searched_read[1][:bclength]
            else:
                searched_read = re.split(bc_anteceder, roi[ad_length:], maxsplit=1)
                if len(searched_read) == 2:
                    barcode = searched_read[0][-bclength:]
                # LT: Increment loss_table if barcode is not found by either flanking sequence
                else:
                    loss_table['bc_flanks'] += 1
            if barcode is not None and len(barcode != bclength:
                # LT: Increment loss_table if barcode is not the correct length
                loss_table['bc_length'] += 1
            if barcode == None or len(barcode) != bclength:
                loss_table['bc_length'] += 1
                barcode = None
                # LT: Increment loss_table to show total number of reads without barcodes
                if barcode == None:
                loss_table['total_bc_not_found'] += 1
            AD = roi[:ad_length]
        else: AD = roi[:ad_length]
    # LT: Increment loss_table if the ad_preceder sequence was not found
    else:
    loss_table['ad_preceder'] += 1
        
    return AD, barcode

def pull_barcode(read, bc_preceder = "GGGCCCG", bc_anteceder = "GGAGAGAA", bclength = 11, **kwargs):
    """Find the barcode in a read.
    
    Takes a read sequence and uses customizable anchor sequences to locate a 
    variable sequence (barcode) in the read.
    
    Parameters
    ----------
    read : str 
        The biological read of interest.
    bc_preceder : str, default "GGGCCCG"
        The anchor sequence directly before the barcode.
    bc_anteceder : str, default "GGAGAGAA"
        The anchor sequence directly after the barcode.
    bc_length : int, default 11
        The length of the barcode sequence if used.
    
    Returns
    ----------
    barcode : str
        The barcode, if used and located. Else None.
    
    Examples
    ----------
    >>> pull_barcode("ACTTTTATVGCTAGCATGGCTGGTAGATCTTGGTTGATTGATTCTAATAGAATTGCTACTAAGATTATGTCTGCTTCTGCTTCTTCTGATCCAAGACAAGTTGTTTGGAAATCTAATCCATCTAGACATTGTCCAGCTGATCGATGCTAGTAGAGAGAGA")
    ATGGCTGGTAGATCTTGGTTGATTGATTCTAATAGAATTGCTACTAAGATTATGTCTGCTTCTGCTTCTTCTGATCCAAGACAAGTTGTTTGGAAATCTAATCCATCTAGACATTGTCCA
    """

    barcode = None
    searched_read = re.split(bc_preceder, read, maxsplit=1)

    if len(searched_read) == 2:
        barcode = searched_read[1][:bclength]

    else:
        searched_read = re.split(bc_anteceder, read, maxsplit=1)
        if len(searched_read) == 2:
            barcode = searched_read[0][-bclength:]
        # LT: Increment loss_table if barcode is not found by either flanking sequence
        else:
            loss_table['bc_flanks'] += 1
    
    return barcode

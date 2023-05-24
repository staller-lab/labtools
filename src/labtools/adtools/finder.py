import re

def pull_AD(read, barcoded = False, ad_preceder = "GCTAGC", 
bc_preceder = "GGGCCCG", bc_anteceder = "GGAGAGAA", ad_length = 120, bclength = 11, **kwargs):
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
                # there might be an error here pls check
                if len(searched_read) == 2:
                    barcode = roi[0][-bclength:]
            if barcode == None or len(barcode) != bclength:
                barcode = None
            else: AD = roi[:ad_length]
        else: AD = roi[:ad_length]
        
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
    
    return barcode
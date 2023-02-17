

def get_ADs(ADs):
    """Take a set of AD objects and return the sequences."""

    for AD in ADs:
        print(AD.count, AD.seq)
        for code in AD.barcodes:
            print(code.seq, code.count)
"""Computer thre reverse complement of a DNA string

Original from http://onlamp.com/pub/a/python/2002/10/17/biopython.html
"""


import string
import sys


sys.stderr.write("use bio_helper instead of revcomp module\n")


# just use [::-1]
#def reverse(string): 
#    """
#    Return a string in reverse order.
#    """
#    
#    letters = list(string) 
#    letters.reverse() 
#    return ''.join(letters) 
     
  
def dna_comp(dnaseq): 
    """
    Return the reverse complement of the dna string.
    """
    
    dna_transtable = string.maketrans("ATCGatcg", "TAGCtagc")
    return dnaseq.translate(dna_transtable)


def dna_revcomp(dnaseq): 
    """
    Return the reverse complement of the dna string.
    """
    return dna_comp(dnaseq)[::-1]


def test():
    """
    Test for correct answer
    """

    dnaseq = 'CCGGAAGAGcttacttag'
    dnarevcomp = dna_revcomp(dnaseq) 
    if dnarevcomp != 'ctaagtaagCTCTTCCGG':
        raise ValueError, "Oops...revcomp gave wrong answer"

    print "Reverse complement of %s is %s." % (
        dnaseq, dnarevcomp)

    
if __name__ == "__main__":
    test()

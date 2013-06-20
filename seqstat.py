#!/usr/bin/env python
"""seqstat.py (a Python implementation of SQUID's seqstat and alistat)
Prints statistics of a sequence file (aligned or unaligned).
"""


#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
from itertools import izip
from itertools import chain
from math import sqrt
from operator import itemgetter

#--- third-party imports
#
from Bio import SeqIO
from Bio import Seq

#--- project specific imports
#
import bioutils


                                                        
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"


# global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


NUCLEIC_ACID_LETTERS = 'ACGTUN'


def argminmax(values, what=None):
    """Return index and value of min/max from values.

    From http://stackoverflow.com/questions/2474015/getting-the-index-of-the-returned-max-or-min-item-using-max-min-on-a-list
    (Matt Anderson's answer)

    Same as numpy.argmin and argmax

    Just to make us independent of numpy        

    Test:
    import numpy
    l = numpy.random.randint(0, 100, 10)
    (numpy.argmin(l), numpy.min(l)) == seqstat.argminmax(l, 'min')
    (numpy.argmax(l), numpy.max(l)) == seqstat.argminmax(l, 'max')
    """
    
    if what == 'min':
        return (min(enumerate(values), key=itemgetter(1)))
    elif what == 'max':
        return (max(enumerate(values), key=itemgetter(1)))
    else:
        raise ValueError, (
            "Only no how to do arg'min' and arg'max' but arg'%s'% what")

    

def meanstd(values):
    """Return mean and standard deviation of values

    Same as: return numpy.mean(values), numpy.stdv(values)

    Just to make us independent of numpy
    
    Test:
    import numpy
    values = numpy.random.randint(0, 100, 1000)
    meanstdv(values) == (numpy.mean(values), numpy.std(values))
    """
    
    size = len(values)
    assert size != 0
    
    mean = sum(values)/float(size)
    std = sqrt(sum([(x - mean)**2 for x in values])/float(size))
    
    return (mean, std)





def guess_if_nucleic_acid(seq, thresh = 0.90,
                          nucleic_acid_letters=NUCLEIC_ACID_LETTERS):
    """Guess if the given sequence is a nucleic acid.

    It's considered nucleic acid if more than 90% of the sequence is
    nucleic_acid_letters (uppercase). The threshold is configurable
    via the thresh parameter.

    Based on http://www.mailinglistarchive.com/html/biopython@biopython.org/2009-01/msg00009.html
    (Author: Brad Chapman)
    """
    
    assert isinstance(seq, Seq.Seq) or \
      isinstance(seq, type("")) or \
      isinstance(seq, type(u""))

    # could use Seq.ungap if Seq.Seq
    seq = bioutils.ungap(str(seq).upper())
    
    nuc_alpha_count = 0
    # don't use collections Counter (Python 2.7 only)
    for letter in nucleic_acid_letters:
        nuc_alpha_count += seq.count(letter)
        #print "DEBUG", len(seq), letter, seq.count(letter)

    if len(seq) == 0:
        return False
    elif float(nuc_alpha_count) / float(len(seq)) >= thresh:
        return False
    else:
        return True

    

def pairwise_identity(s1, s2):
    """Return fractional pairwise identity between two aligned
    strings, which is defined here as the number of identical residues
    (case sensitive), divived by the smaller of the two ungapped
    sequences.

    Uppercase your sequence for case insensitivity. For mixed RNA/DNA
    you might want to replace T's with U's vice versa.
    
    Based on ideas from
    http://code.activestate.com/recipes/499304-hamming-distance/
    """
    
    assert len(s1) == len(s2)
    idents = sum(c1 == c2
                 for c1, c2 in izip(s1, s2) 
                 if not bioutils.isgap(c1) and not bioutils.isgap(c2))
    min_ungapped_len = min(len(bioutils.ungap(s1)), len(bioutils.ungap(s2)))
    return idents / float(min_ungapped_len)




def comp_pairwise_ident_matrix(seqrecs):
    """Returns a fake matrix (symmetric 2d list) of pairwise
    identities. Valid index range is [i][j], where i>=j, j>=0 and
    i<nseqs. values for i=j are None!
    """
    nseqs = len(seqrecs)

    # intentionally a list, not a matrix, because numpy doesn't know
    # about symmetric arrays
    mx = []
    for i in xrange(nseqs):
        jdists = []
        for j in xrange(0, i):
            s1 = str(seqrecs[i].seq).upper()
            s2 = str(seqrecs[j].seq).upper()
            pwid = pairwise_identity(s1, s2)
            jdists.append(pwid)

            if False:
                # tmp hack dna dist
                dist = sum(c1 != c2
                         for c1, c2 in izip(s1, s2)
                         if not bioutils.isgap(c1) and not bioutils.isgap(c2))
                print "TMP: dist %s vs %s: %d" % (seqrecs[i].id, seqrecs[j].id, dist)

        jdists.append(None) # self comparison not defined
        mx.append(jdists)
    return mx
    


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--verbose",
                      dest="verbose",
                      action="store_true",
                      help=SUPPRESS_HELP) #"be verbose")
    parser.add_option("", "--debug",
                      dest="debug",
                      action="store_true", 
                      help=SUPPRESS_HELP) #"debugging")
    parser.add_option("-a", "--all",
                      action="store_true", 
                      dest="info_for_all",
                      help="Report information for all sequences")
    parser.add_option("-f", "--infmt",
                      dest="informat",
                      help="Input format (must be supported by Biopython)")
    return parser



def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)
        
    if len(args) != 1:
        parser.error("Need sequence file as input argument")
        sys.exit(1)

        
    fseq = args[0]
    if fseq == "-":
        fhandle = sys.stdin
    else:
        fhandle = open(fseq, "rU")


    fmt = opts.informat
    if not fmt:
        fmt = bioutils.guess_seqformat(fseq)


    seqrecs = []
    seqlens = []
    seqlens_ungapped = []
    
    # read all into memory: makes id computation easier. the file
    # might come from stdin so we can't read twice
    for seqrec in SeqIO.parse(fhandle, fmt):
        seqrecs.append(seqrec)
        seqlens.append(len(seqrec.seq))
        seqlens_ungapped.append(len(bioutils.ungap(str(seqrec.seq))))
    if fhandle != sys.stdin:
        fhandle.close()
            
    nseqs = len(seqlens)        
    if nseqs == 0:
        LOG.warn('No sequences found. Try changing the format (just tried: %s)' % fmt)
        sys.exit(0)

    
    aligned = False
    if nseqs > 1 and len(set(seqlens)) == 1:
        # add and len(set(seqlens_ungapped)) != 1 to make sure
        # unaligend sequence length are identical
        aligned = True
        aln_len = seqlens[0] # any will do as we know they're aligned
        pw_id_mx = comp_pairwise_ident_matrix(seqrecs)

    if not aligned and seqlens != seqlens_ungapped:
        LOG.warn("Found gaps, but sequences do not seem to be aligned."
                 " Stats will be for ungapped seqs.")
         
    # guess type from first entry
    if guess_if_nucleic_acid(seqrecs[0].seq):
        seqtype = 'protein' 
    else:
        seqtype = 'nucleic'
    print "Type (of 1st seq):   %s" % (seqtype)

    print "Number of sequences: %d" % (nseqs)        
    print "Smallest:            %d" % (
        min(seqlens_ungapped))
    print "Largest:             %d" % (
        max(seqlens_ungapped))
    print "Average length:      %.1f" % (
        sum(seqlens_ungapped)/float(len(seqlens_ungapped)))
    #print "Format:              %s" % (fmt)
    
    print "Aligned:             %s" % ("yes" if aligned else "no")
    if aligned:
        # make sure to ignore self-comparison None's
        flat_pw_id_mx = [x for x in chain.from_iterable(pw_id_mx) if x]
        print "Alignment length:    %d" % (aln_len)        
        (mean, std) = meanstd(flat_pw_id_mx)
        print "Average identity:    %0.2f" % (
            mean)
        print "Standard deviation:  %0.2f" % (
            std)
        print "Most related pair:   %0.2f" % (
            max(flat_pw_id_mx))
        print "Most unrelated pair: %0.2f" % (
            min(flat_pw_id_mx))
    
    if opts.info_for_all:
        # spacer
        print ""
        
        header = "# Name\tLength"
        if aligned:
            header += "\thigh-id to\tlow-id to"
        print header
        
        for (i, seqrec) in enumerate(seqrecs):
            line = "%s\t%d" % (
                seqrec.id, seqlens_ungapped[i])
            
            if aligned:
                # construct list of pairwise ids from fake matrix. 
                pw_ids = pw_id_mx[i]
                pw_ids.extend([pw_id_mx[j][i] 
                               for j in xrange(i+1, nseqs)])
                assert len(pw_ids) == nseqs, (
                    "len(pw_ids)=%d, but expected %d" % (len(pw_ids), nseqs))

                # Find min and max and corresponding partner index,
                # but take care to ignore self-comparison value 'None'
                pw_ids[i] = -1.0
                (pw_id_max_idx, pw_id_max_val) = argminmax(pw_ids, 'max')
                pw_ids[i] = 1.1
                (pw_id_min_idx, pw_id_min_val) = argminmax(pw_ids, 'min')
                pw_ids[i] = None # reset even though not strictly necessary

                line += "\t%.4f\t%s\t%.4f\t%s" % (
                    pw_id_max_val, seqrecs[pw_id_max_idx].id,
                    pw_id_min_val, seqrecs[pw_id_min_idx].id)
            print line

    print "%d names are unique and %d sequences are unique (including gaps)." % (
        len(set([s.id for s in seqrecs])),
        len(set([str(s.seq) for s in seqrecs])))

            

if __name__ == "__main__":
    if sys.version_info < (2 , 6) or sys.version_info > (2 , 8):
        sys.stderr.write(
            "WARNING: Python version %s untested\n" % sys.version_info)

    #biopython_version = tuple([int(x) for x in Bio.__version__.split('.')])
    #if biopython_version < (1 , 55) or biopython_version > (1 , 59):
    #    sys.stderr.write(
    #        "WARNING: Biopython version %s untested\n" % (Bio.__version__))

    main()
    LOG.info("Successful exit")

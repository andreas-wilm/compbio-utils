#!/usr/bin/env python

"""Extract and prints annotations from reference Genbank file,
including their positions and corresponding positions for a query
sequence in given pairwise alignment"""


#--- standard library imports
#
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser, SUPPRESS_HELP
import difflib

#from collections import namedtuple
#Annotation = namedtuple('Annotation', ['start', 'end', 'type', 'descr'])

#--- third-party imports
#
import Bio
from Bio import SeqIO


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
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

class PosMap(object):
    """Position map class

    NOTE: all unit-offset!
    """
    
    def __init__(self, seqrecs=None):
        """
        """

        self.seq_ids = []
        self.pos_map = dict()

        if seqrecs:
            self.generate(seqrecs)
        
    @staticmethod
    def isgap(res):
        """Return true if given residue is a gap character
        """
        return (res in ['-', '~', '.'])
    

    
    def generate(self, seqrecs):
        """Computes a position map, which is a dict with aligned
        positions as main key. Sequence ids are 2nd dim key and their
        corresponding unaligned position is the value
        
        NOTE: if a residue is aligned to a gap, the previous position
        is used
    
        NOTE: the format is terribly inefficient. should be spit out
        as blocks/ranges asfor liftover chains (troublesome for >2
        though)")
        """
       
        self.pos_map = dict()
        self.seq_ids = [s.id for s in seqrecs]
        
        aln_len = len(seqrecs[0].seq)
        for s in seqrecs:
            assert len(s.seq) == aln_len, (
                "Looks like your seqs are not aligned")

        # all offset one
        cur_unaligned_pos = len(seqrecs) * [0]
        for aln_pos in xrange(aln_len):
            for s in xrange(len(seqrecs)):
                res = seqrecs[s][aln_pos]
                if not self.isgap(res):
                    cur_unaligned_pos[s] += 1
            self.pos_map[aln_pos+1] = dict()
            for (p, s) in zip(cur_unaligned_pos, [s.id for s in seqrecs]):
                self.pos_map[aln_pos+1][s] = p

                
    
    def output(self, fh=sys.stdout):
        """Print position map
        """

        print "aln-pos\t%s" % ('\t'.join(self.seq_ids))
        for aln_pos in sorted(self.pos_map.keys()):
            line = "%d" % aln_pos
            for s in self.seq_ids:
                line += " %d" % (self.pos_map[aln_pos][s])
            fh.write("%s\n" % (line))
            

            
    def parse(self, pos_map_file):
        """Parse position map from file
        """

        LOG.critical("Untested function")
        fh = open(pos_map_file, 'r')
        line = fh.readline()
        header = line.rstrip().split('\t')
        assert header[0] == 'aln-pos', (
            "Was expecting first field to be aln-pos, but it's '%s'" % (
                header[0]))
    
        self.pos_map = dict()
        for line in fh:
            # note: offset untouched, i.e. as in file (unit-offset)
            positions = [int(x) for x in line.rstrip().split('\t')]
            assert len(positions) == len(header)
    
            aln_pos = positions[0]
            self.pos_map[aln_pos] = dict()
            for (p, s) in zip(positions[1:], header[1:]):
                self.pos_map[aln_pos][s] = p
        fh.close()
    
    

    def convert(self, src=None, query=None):
        """Mangles input pos_map and returns a dict containing unaligned
        position matching between src and query ids.

        If query is None, then aligned positions for src are returned.

        Likewise if src is None, then unaligned positions for aligned pos are returned
        """
    
        # FIXME Would this break if we had gap vs gap alignment and
        # some residues following?

                
        if query and src:
            d = dict([(v[src], v[query]) 
                      for (k, v) in self.pos_map.iteritems()])
        elif src:
            d = dict([(v[src], k) 
                      for (k, v) in self.pos_map.iteritems()])
            #for (k, v) in d.iteritems():
            #    assert k<=v
        elif query:            
            d = dict([(k, v[query]) 
                      for (k, v) in self.pos_map.iteritems()])            
        else:
            raise ValueError
            
        return d


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
    parser.add_option("-r", "--ref-gb",
                      dest="ref_gb",
                      help="Reference Genbank file")
    parser.add_option("-p", "--pw-aln",
                      dest="pw_aln",
                      help="Pairwise alignment including the refseq given with --ref-gb")
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
        
    if len(args) != 0:
        parser.error("Unrecognized args found")
        sys.exit(1)

    for (f, d) in [(opts.pw_aln, "Pairwise alignment"),
                   (opts.ref_gb, "Reference Genbank")]:
        if not f:
            parser.error("Missing %s argument" % d)
        if not os.path.exists(f):
            LOG.fatal("%s file '%s' does not exist" % (d, f))
            sys.exit(1)


    refseq = list(SeqIO.parse(opts.ref_gb, "genbank"))
    assert len(refseq)==1
    refseq = refseq[0]


    pw_aln = list(SeqIO.parse(opts.pw_aln, bioutils.guess_seqformat(opts.pw_aln)))
    assert len(pw_aln)==2, (
        "Was expecting two sequences, but parsed %d from %s" % (
            len(pw_aln), opts.pw_aln))


    # determine ref id
    #
    # seqids in alignment should match genbank id but might not
    matches = difflib.get_close_matches(refseq.id, [s.id for s in pw_aln])
    assert len(matches), (
        "Couldn't find a sensible match between sequence ids in alignment and genbank")
    aln_ref_id = matches[0]
    if aln_ref_id != refseq.id:
        LOG.warn("Assuming %s (from alignment) is the same as %s (from Genbank)" % (
            aln_ref_id, refseq.id))
    LOG.info("%s is the ref id" % (aln_ref_id))

    # determine query id
    assert len(pw_aln) == 2
    for sid in [s.id for s in pw_aln]:
        if sid != aln_ref_id:
            query_id = sid
    LOG.info("%s is the query id" % (query_id))
    
    pos_map = PosMap(pw_aln)
    #pos_map.output()
    ref_to_query_map = pos_map.convert(aln_ref_id, query_id)

    print "#QUERY-POS (%s)\tREF-POS (%s)\tSTRAND\tTYPE\tQUALIFIERS" % (
        query_id, aln_ref_id)
    for feat in refseq.features:
        query_pos_str =  "%d-%d" % (
            ref_to_query_map[feat.location.start.position+1],
            ref_to_query_map[feat.location.end.position])
        orig_pos_str = "%d-%d" % (
            feat.location.start.position+1,
            feat.location.end.position)
        strand_str = "%s" % feat.strand
        type_str = "%s" % feat.type
        qualifiers_str = '; '.join("%s %s" % (k, ', '.join(v))
                                    for (k, v) in feat.qualifiers.iteritems() 
                                    if k != 'translation')
        print '\t'.join([query_pos_str, orig_pos_str, 
                         strand_str, type_str, qualifiers_str])
    LOG.info("Note: feature positions might overlap")
        
            
    
if __name__ == "__main__":
    if sys.version_info < (2 , 7):
        sys.stderr.write("WARNING: only tested Python 2.7 so far\n")
    elif sys.version_info > (2 , 8):
        sys.stderr.write("WARNING: only tested Python 2.7 so far\n")
        
    biopython_version = tuple([int(x) for x in Bio.__version__.split('.')])
    if biopython_version < (1 , 55):
        sys.stderr.write("WARNING: only tested Biopython 1.55 so far\n")
    elif biopython_version > (1 , 55):
        sys.stderr.write("WARNING: only tested Biopython 1.55 so far\n")
        
    main()
    LOG.info("Successful exit")

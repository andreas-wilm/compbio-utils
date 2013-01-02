#!/usr/bin/env python
"""Take contigs and orient them accordint to a reference using Mummer.
Optionally fill in the missing bits with the reference itself instead
of Ns """

#--- standard library imports
#
import sys
import os
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import subprocess
import tempfile
import shutil

# invocation of ipython on exceptions
if False:
    import sys, pdb
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(
        mode='Verbose', color_scheme='Linux', call_pdb=1)


#--- third-party imports
#
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


#--- project specific imports
#
# /

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""



# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')



def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: " + __doc__ + "\n" \
                  "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="debugging")
    parser.add_option("-c", "--c",
                      dest="fcontigs",
                       help="Input file containing contigs to join (fasta)")
    parser.add_option("-n", "--dont-fill-with-ref",
                      dest="dont_fill_with_ref", 
                      action="store_true",
                      help="Don't fill gaps with reference, but keep Ns isntead")
    parser.add_option("-r", "--ref",
                      dest="fref",
                      help="Input file containing reference sequence to fill gaps between contigs (fasta)")
    parser.add_option("-o", "--output",
                      dest="fout",
                      help="output file (fasta)")
    parser.add_option("", "--keep-tmp-files",
                      dest="keep_temp_files",
                      action="store_true",
                      default=False,
                      help="Optional: Don't delete temp (nucmer etc) files")
    parser.add_option("", "--tmp-dir",
                      dest="tmp_dir", # type="string|int|float"
                      help="Optional: directory to save temp files in")

    return parser




def main():
    """
    The main function
    """
    
    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)
        
    if opts.verbose:
        LOG.setLevel(logging.INFO)

    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    if not opts.fref:
        parser.error("Reference input file argument missing.")
        sys.exit(1)
    if not opts.fcontigs:
        parser.error("Contig input file argument missing.")
        sys.exit(1)

    for fname in [opts.fref, opts.fcontigs]:
        if not os.path.exists(fname):
            LOG.fatal(
                "file '%s' does not exist.\n" % fname)
            sys.exit(1)

    if not opts.fout:
        parser.error("Fasta output file argument missing.")
        sys.exit(1)
    if os.path.exists(opts.fout):
        LOG.fatal(
            "Refusing to overwrite existing file '%s'." % opts.fout)
        sys.exit(1)
        

    tmp_files = []

    # run mummer's nucmer
    #
    out_prefix = tempfile.NamedTemporaryFile(
        delete=False, dir=opts.tmp_dir).name
    fdelta = out_prefix + ".delta"
    cmd_args = ['nucmer', opts.fref, opts.fcontigs,  
                '-p', out_prefix]
    process = subprocess.Popen(cmd_args, 
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (p_stdout, p_stderr) =  process.communicate()
    assert os.path.exists(fdelta), (
        "Couldn't find expected output ('%s') for command: '%s'" % (
            fdelta, ' '.join(cmd_args)))
    tmp_files.append(fdelta)
    LOG.info("Delta written to %s" % fdelta)


    # run mummer's show-tiling which creates a pseudo molecule from
    # the contigs, filled with Ns
    #
    fpseudo = fdelta + "_pseudo.fa"
    cmd_args = ['show-tiling', fdelta, '-p', fpseudo]
    process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (p_stdout, p_stderr) =  process.communicate()
    assert os.path.exists(fpseudo), (
        "Couldn't find expected output ('%s') for command: '%s'" % (
            fpseudo, ' '.join(cmd_args)))
    tmp_files.append(fpseudo)
    LOG.info("Pseudo written to %s" % fpseudo)

    
    if opts.dont_fill_with_ref:
        LOG.info(
            "Early exit as requested: keeping nucmer's stitched version using N's instead of ref and copying to %s" % (
                opts.fout))
        shutil.copyfile(fpseudo, opts.fout)
        if not opts.keep_temp_files:
            for f in tmp_files:
                os.unlink(f)
        sys.exit(0)
        
    
    # concatenate the pseudo file and the reference
    #
    faln_in = fpseudo + "_plus_ref.fa"
    fh = open(faln_in, 'w')
    for f in [opts.fref, fpseudo]:
        seqrecs = list(SeqIO.parse(f, "fasta"))
        assert len(seqrecs)==1, (
            "Expected exactly one sequence in '%s'" % f)
        SeqIO.write(seqrecs, fh, "fasta")
    fh.close()
    tmp_files.append(faln_in)
    

    # align pseudo file (gapped contigs) and reference
    #
    faln_out = fpseudo + "_plus_ref_aln.fa"                                                                                  
    #cmd_args = ['muscle', '-maxiters', '1', '-in', faln_in, '-out', faln_out]
    #process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE,
    #                           stderr=subprocess.PIPE)
    cmd = 'mafft --auto --anysymbol %s > %s' % (faln_in,  faln_out)
    process = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    (p_stdout, p_stderr) =  process.communicate()
    assert os.path.exists(faln_out), (
        "Couldn't find expected output ('%s') for command: '%s'" % (
            faln_out, ' '.join(cmd_args)))
    tmp_files.append(faln_out)
    LOG.info("Alignment written to %s" % faln_out)

    fh = open(faln_out, "r")
    seqrecs = list(SeqIO.parse(fh, "fasta"))
    #cleaned_stderr = ''.join(process.stderr.readlines()).replace("WARNING : Unknown character n\r", "")
    fh.close()
    assert len(seqrecs)==2, (
        "Expected exactly two sequences in alignment %s but got %d." % (
            faln_out, len(seqrecs)))


    # now walk through the pairwise alignment. keep the pseudo contig
    # and use the reference only in case of gaps (marked as Ns) or
    # other ambiguity characters
    #
    stitched_seq = []
    # both seqrecs are aligned
    (ref_seqrec, pseudo_seqrec) = sorted(seqrecs, 
                                         key = lambda s: s.id.startswith("pseudo_"))
    for pos in range(len(pseudo_seqrec.seq)):
        pseudo_nt = pseudo_seqrec.seq[pos]
        ref_nt = ref_seqrec.seq[pos]

        # N: "real" N in contig or to be replaced with reference
        # internal gap: deletion with respect to reference
        # terminal gap: missing fragment
        # delay gap handling therefore
        # tilde in contigs will be kept as indel
        if pseudo_nt == 'N' or pseudo_nt == "-":
            stitched_seq.append(ref_nt)
        elif pseudo_nt.upper() not in 'ACGTU~' and ref_nt.upper() in 'ACGTU':
            stitched_seq.append(ref_nt)
        else:
            stitched_seq.append(pseudo_nt)

    # delayed: gap handling: substitute terminal ones with refseq and
    # replace internal ones

    stitched_seq = ''.join(stitched_seq)

    tmp_len = len(stitched_seq)
    stitched_seq = stitched_seq.lstrip("-")
    num_gaps_start = tmp_len - len(stitched_seq)

    tmp_len = len(stitched_seq)
    stitched_seq = stitched_seq.rstrip("-")
    num_gaps_end = tmp_len - len(stitched_seq)

    stitched_seq  = "%s%s%s" % (ref_seqrec.seq[:num_gaps_start],
                                stitched_seq,
                                ref_seqrec.seq[:-num_gaps_end])

    stitched_seq = stitched_seq.replace("-", "")
    stitched_seq = stitched_seq.replace("~", "")
    
    # output
    #
    out_seqrec = SeqRecord(Seq(stitched_seq, IUPAC.ambiguous_dna),
                           id="stitched_seq", 
                           description="Mosaic of %s and contigs" % (ref_seqrec.id))
    SeqIO.write(out_seqrec, opts.fout, "fasta")
    LOG.info("Stitched/Mosaic sequence with fake id %s written to %s" % (
            out_seqrec.id, opts.fout))

    if not opts.keep_temp_files:
        for f in tmp_files:
            os.unlink(f)

if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
    LOG.warn("Best to check the output through alignment")
    LOG.warn("Tilde in contigs are kept as insertions")

#!/usr/bin/env python
"""Wrapper for benchmarking alignment programs on Prefab
"""


#--- standard library imports
#
import sys
import os
from optparse import OptionParser;# NOTE: deprecated since python 2.7
import logging
#import tempfile
import shutil

#--- third-party imports
#

#--- project specific imports
#


__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


#MY_DIR = os.path.dirname(sys.argv[0])
MY_DIR = os.path.abspath(os.path.dirname(sys.argv[0]))
MY_NAME = os.path.basename(sys.argv[0])
DEFAULT_PREFAB_FILE_LIST = os.path.join(MY_DIR, "prefab4_file_list.txt")
#QSCORE_BIN = os.path.join(MY_DIR, "qscore", "qscore")
QSCORE_BIN = "qscore"
PREFAB_DIR = os.path.join(MY_DIR, "prefab4")
PREFAB_IN_DIR = os.path.join(PREFAB_DIR, "in")
PREFAB_REF_DIR = os.path.join(PREFAB_DIR, "ref")
SREFORMAT_BIN = "sreformat"



def binary_in_path(binary):
    """Tests if binary is detected by 'which'

    Arguments:
    - `binary`: the binary to test for

    Results:
    - Return True of found, False otherwise
    """
    
    cmd = "which %s >/dev/null 2>/dev/null" % (binary)
    result = os.system(cmd)
    if result:
        return False
    else:
        return True


def parse_prefab_listing(filename):
    """Parses the prefab file listing

    Arguments:
    - `filename`: the listing filename

    Results:

    - Returns a list with dictionaries for each entry in filename. each dictionary has the keys
    basename, numseq and pwid.
    """

    ret_listing = list()
    
    fid = open(filename, "rU")
    for line in fid:
        # ignore comments
        if line.startswith("#"):
            continue
        assert 3 == len(line.split()), (
            "Expected format: 'file numseq pwid' in '%s'" % filename)
            
        (basename, numseq, pwid) = line.split()
        numseq = int(numseq)
        pwid = int(pwid)
        assert numseq > 1, (
            "Parsed number of sequences (%d) is not >1" % numseq)
        assert pwid >= 1 and pwid < 100, (
            "Parsed pairwise identity (%d) should be >0 and <100" % pwid)

        ret_listing.append(
            dict(basename=basename, numseq=numseq, pwid=pwid)
            )
    
    fid.close()

    return ret_listing
    

    
def cmdline_parser():
    """Creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog: A PREFAB wrapper\n" \
            "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("", "--debug",
                      action="store_true",
                      dest="debug",
                      help="switch debugging on")
    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_option("-o", "--outdir",
                      dest="outdir",
                      help="Directory where files will be stored")

    # limits
    parser.add_option("", "--max-nseq",
                      dest="max_nseq",
                      default=sys.maxint,
                      type=int,
                      help="Skip files with more than this number of sequences")
    parser.add_option("", "--min-nseq",
                      default=2,
                      type=int,
                      dest="min_nseq",
                      help="Skip files with less than this number of sequences")
    parser.add_option("", "--max-id",
                      default=100,
                      type=int,
                      dest="max_pwid",
                      help="Skip files whose pairwise is higher than this value")
    parser.add_option("", "--min-id",
                      default=1,
                      dest="min_pwid",
                      type=int,
                      help="Skip files whose pairwise is lower than this value")

    parser.add_option("-l", "--filelist",
                      dest="filelist",
                      default=DEFAULT_PREFAB_FILE_LIST,
                      help="Alternative file listing (default is %s)" % DEFAULT_PREFAB_FILE_LIST)
    parser.add_option("-c", "--cmd",
                      dest="cmd",
                      help="%s" % ''.join([
                          "command to run.",
                          " @IN@ and @OUT@ will be replaced accordingly.",
                          " If @OUT@ is missing,",
                          " alignment is expected to be printed to stdout by command."
                          ])
                      )
    return parser



def qscore_aln(faln, fref, fscore):
    """Reformats faln, runs qscore and appends result to fscore


    Arguments:
    - `faln`: test alignment file  
    - `fref`: reference alignment file
    - `fscore`: output score file
    
    Results:
    - Returns False on error, True otherwise
    - fscore will be created
    """

    faln_fasta = faln + "_sreformat.fa"

    cmd = "%s a2m %s >> %s" % (SREFORMAT_BIN, faln, faln_fasta)
    result = os.system(cmd)
    if result:
        LOG.critical("sreformat failed with exit code '%s'. Command was '%s'" % (
            result, cmd))
        return False

    # could us -seqdiffwarn but would have to redirect stderr
    # FIXME: do we need -ignoretestcase?
    cmd = "%s -test %s -ref %s > %s" % (
        QSCORE_BIN, faln_fasta, fref, fscore)
    
    result = os.system(cmd)
    if result:
        LOG.critical("qscore failed with exit code '%s'. Command was '%s')" % (
            result, cmd))
        return False

    # remove reformatted copy of alignment file to save space
    os.remove(faln_fasta)

    return True



def run_aln_cmd(cmd, workdir, infile_src):
    """

    Arguments:
    - `cmd`: command to run, including @IN@ and @OUT@ markup
    - `workdir`: dir to run it
    - `infile_src`: full path to input sequence file
    
    Results:
    - Runs cmd in workdir with @IN@ and @OUT@ replaced
    - Returns output filename
    """

    assert "@IN@" in cmd, (
        "Can't find @IN@ markup in command.")


    # default filenames and commaneds
    #
    aln_raw_out = os.path.basename(infile_src) + ".out"
    readme = "README"
    time_cmd = "/usr/bin/time -p"
    stdout_log = "stdout.log"
    stderr_log = "stderr.log"


    # work on a copy of infile which will be deleted later
    #
    origdir = os.getcwd()
    shutil.copy(infile_src, workdir)
    os.chdir(workdir)
    infile = os.path.basename(infile_src)

 
    # replace @IN@ and @OUT@ in cmd
    #
    cmd = cmd.replace("@IN@", infile)
    if "@OUT@" in cmd:
        cmd = cmd.replace("@OUT@", aln_raw_out)

    # create README file
    #
    fid = open(readme, 'w')
    fid.write("This file was created by '%s'\n" % \
                  ' '.join(sys.argv))
    fid.write("Will now try to execute: %s\n" % cmd)
    fid.write("stdout and stderr can be found in %s and/or %s\n" % (
        stdout_log, stderr_log))
    fid.close()


    # run aligner
    #
    # if @OUT@ was not part of command, assume that the alignment goes
    # to stdout instead
    LOG.debug('Running "%s" in %s' % (cmd, workdir))
    if aln_raw_out in cmd:
        result = os.system("%s %s > %s 2>%s" % (
            time_cmd, cmd, stdout_log, stderr_log))
    else:
        result = os.system("%s %s > %s 2>%s" % (
            time_cmd, cmd, aln_raw_out, stderr_log))

    # check result and output
    #
    err_msg = None
    if result:
        err_msg = "'%s' failed with error code '%s'." % (cmd, result)        
    if not os.path.exists(aln_raw_out) or os.path.getsize(aln_raw_out) == 0:
        err_msg = "'%s' produced no output." % (cmd)
        
    if err_msg:
        LOG.critical("%s Check log-files in '%s'" % \
                      (err_msg, workdir))
        fid = open(readme, 'a')
        fid.write(err_msg)
        fid.close()

        os.chdir(origdir)
        return None

    # remove copy of infile
    os.remove(infile)

    os.chdir(origdir)
    return os.path.join(workdir, aln_raw_out)
        


def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()
    
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    if opts.verbose:
        LOG.setLevel(logging.INFO)

    if not opts.cmd:
        parser.error("Missing command argument")
        sys.exit(1)
    if '@IN@' not in opts.cmd:
        parser.error("Can't find @IN@ markup in command '%s'" % (
            opts.cmd))

    if not opts.outdir:
        parser.error("Missing output directory argument")
        sys.exit(1)
    #if os.path.isdir(opts.outdir):
    #    parser.error("Cowardly refusing to use already existant output directory '%s'" % (
    #        opts.outdir))
    #    sys.exit(1)


    # startup checks
    #
    for dirname in [PREFAB_DIR, PREFAB_IN_DIR, PREFAB_REF_DIR]:
        if not os.path.isdir(dirname):
            LOG.fatal("Missing directory: '%s'" % dirname)
            sys.exit(1)
    for fname in [DEFAULT_PREFAB_FILE_LIST]:
        if not os.path.exists(fname):
            LOG.fatal("Missing file: '%s'" % fname)
            sys.exit(1)
    for binary in [SREFORMAT_BIN, QSCORE_BIN]:
        if not binary_in_path(binary):
            LOG.fatal("Couldn't find binary '%s'" % (binary))
            sys.exit(1)
            

    try:
        os.mkdir(opts.outdir)
    except OSError:
        # FIXME ok if basedir exists, as long as subdirs don't
        pass

                     
    prefab_listing = parse_prefab_listing(opts.filelist)
    org_num_entries = len(prefab_listing)
    #print "pdb after parse_prefab_listing"; import pdb; pdb.set_trace()


    # filter as requested
    #
    
    prefab_listing = filter(lambda x: x['numseq'] >= opts.min_nseq,
                            prefab_listing)
    prefab_listing = filter(lambda x: x['numseq'] <= opts.max_nseq,
                            prefab_listing)
    prefab_listing = filter(lambda x: x['pwid'] >= opts.min_pwid,
                            prefab_listing)
    prefab_listing = filter(lambda x: x['pwid'] <= opts.max_pwid,
                            prefab_listing)


    LOG.info("Will try to align %d out of %d entries" % (
            len(prefab_listing), org_num_entries))
    for (ctr, entry) in enumerate(prefab_listing):
        
        fin = os.path.join(PREFAB_IN_DIR, entry['basename'])

        entry_outdir = os.path.join(opts.outdir, entry['basename'])
        try:
            os.mkdir(entry_outdir)
        except OSError:
            LOG.critical(
                "Cowardly refusing to overwrite existant dir '%s'" % (
                entry_outdir))
            continue
        
        LOG.info("Step %d of %d: Aligning %s (#seq=%d, pwid=%d%%) in %s" % (
            ctr+1, len(prefab_listing),
            entry['basename'],
            entry['numseq'], entry['pwid'],
            entry_outdir))
        
        fout = run_aln_cmd(opts.cmd, entry_outdir, fin)
        if not fout:
            # error message already printed
            continue

        fref = os.path.join(PREFAB_REF_DIR, entry['basename'])
        fscore = fout + "_qscore.txt"
        if qscore_aln(fout, fref, fscore):
            LOG.info("qscore result stored in %s" % fscore)
            
if __name__ == "__main__":
    main()

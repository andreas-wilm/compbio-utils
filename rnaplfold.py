

from subprocess import Popen, PIPE, STDOUT
import os
import tempfile

def rnaplfold(seq, prob_threshold=0.0, rnaplfold_binary="RNAplfold"):
    """
    
    Arguments:
    - `seq`: FIXME
    - `prob_threshold`: FIXME
    """
    
    # There are two ways to call RNAplfold:
    #
    # echo seq | RNAplfold; produces plfold_dp.ps
    # or
    # RNAplfold < seq.file; produces seqid__dp.ps
    # 
    # We use the first to avoid issues with messy seq-ids

    RNAPLFOLD_DEFAULT_DPPS_NAME = "plfold_dp.ps"
    
    # RNAplfold produces a file in current dir, so move to temp dir
    # before calling
    orig_workdir = os.getcwd()
    os.chdir(tempfile.tempdir)

    p = Popen([rnaplfold_binary], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    # gives OSError if not installed
    (stdoutdata, stderrdata) = p.communicate(input=seq + "\n")[0]    
    # FIXME test stderrdata and stdoutdata?
    # FIXME test if file exists

    fid_dpps = open(RNAPLFOLD_DEFAULT_DPPS_NAME, 'r')
    for line in fid_dpps:
        # base pairs are stored in the following form:
        # % i  j  sqrt(p(i,j)) ubox
        line_split = line.split()
        if len(line_split) == 4 and line_split[3] == "ubox":
            print "DEBUG: got a bp (FIXME: unit offset?) at %s" % ':'.join(line.split()[0:2])
            
            #(bp_i, bp_j) = [int(x) for x in line_split[0:2]
            #bp_prob = float(line_split[2])**2

            # FIXME store as list of lists
            # FIXME zero-/unit-offset?            
            # FIXME apply prob_threshold

    fid_dpps.close()

    os.remove(RNAPLFOLD_DEFAULT_DPPS_NAME)

    os.chdir(orig_workdir)

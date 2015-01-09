#!/usr/bin/env python
"""FIXME:add-doc
"""


import sys
import gzip
import os
from itertools import chain
from collections import Counter
#import pickle
import cPickle as pickle

# FIXME hardcoded lofreq path for vcf import
#sys.path.insert(0, '/home/wilma/lofreq/lofreq2-gis.git/src/lofreq_python/')
#from lofreq_star import vcf
import vcf

def load_mutect_ext(mutect_ext_file):
    """FIXME:add-doc
    """
    
    mutect_rej = dict()
    mutect_ext_pkl = mutect_ext_file + ".pkl"

    # load from pickled file for speedup
    if os.path.exists(mutect_ext_pkl):
        assert os.path.getctime(mutect_ext_file) < os.path.getctime(mutect_ext_pkl)
        sys.stderr.write("Loading data from pickle file instead\n")
        fh = open(mutect_ext_pkl, 'rb')
        mutect_rej = pickle.load(fh)
        fh.close()
        return mutect_rej

    if mutect_ext_file == "-":
        fh = sys.stdin
    elif mutect_ext_file[-3:] == ".gz":
        fh = gzip.open(mutect_ext_file)
    else:
        fh = open(mutect_ext_file)
        
    for line in fh:
        line = line.strip()
        if line.startswith('#') or len(line)==0:
            continue
        lsplit = line.split('\t')
        if lsplit[-1] != 'REJECT':
            #print "DEBUG: ignoring %s" % line
            continue
        key = "%s,%s" % (lsplit[0], lsplit[1])
        reasons = lsplit[-2].split(',')
        mutect_rej[key] = reasons
        
    if fh != sys.stdin:
        fh.close()

        # try dumping pkl file
        # should not get here if one existed so there's no risk overwriting it
        sys.stderr.write("Creating pickle file\n")
        assert not os.path.exists(mutect_ext_pkl)
        fh = open(mutect_ext_pkl, 'wb')
        pickle.dump(mutect_rej, fh)
        fh.close()

    return mutect_rej


def main():
    """main function
    """
    try:
        mutect_ext_file = sys.argv[1]
    except IndexError:
        sys.stderr.write("%s: report mutect rejection categories\n" % os.path.basename(sys.argv[0]))
        sys.stderr.write("usage: %s mutect-ext[.gz] [vcf[.gz]]\n" % os.path.basename(sys.argv[0]))
        sys.stderr.write("If a vcf-file is given then rejection categories for variants listed in this file are reportd\n")
        sys.exit(1)
        
    try:
        vcf_file = sys.argv[2]
    except IndexError:
        vcf_file = None

    mutect_rej = load_mutect_ext(mutect_ext_file)    
    print "INFO: Loaded %d rejected variants from Mutect file %s" % (
        len(mutect_rej), mutect_ext_file)

    cat_count = Counter(chain.from_iterable(mutect_rej.values()))


    if not vcf_file:
        print "#rejection category\tfreq"
        total_count = sum(cat_count.values())
        for (cat, count) in sorted(cat_count.items(), key=lambda x:x[1]):
            print "%s\t%.6f" % (cat, count/float(total_count))
    
    else:
        if vcf_file[-3:] == ".gz":
            vcf_fh = gzip.open(vcf_file)
        else:
            vcf_fh = open(vcf_file)
        vcfreader = vcf.VCFReader(vcf_fh)
        lofreq_uniq_vars =  [v for v in vcfreader]
        vcf_fh.close()
        print "INFO: Loaded %d variants from %s" % (len(lofreq_uniq_vars), vcf_file)
        print "INFO: Reporting corresponding Mutect rejection categories (NA = not found)" 
        
        rcounts = dict()
        s = 0
        for v in lofreq_uniq_vars:
            key = "%s,%d" % (v.CHROM, v.POS)
            
            if not mutect_rej.has_key(key):
                rcounts['NA'] = rcounts.get('NA', 0) + 1
                s += 1
                #print "DEBUG: %s not listed in mutect-ext" % (key)
                continue
            
            #print "DEBUG: %s %s" % (key, ', '.join(mutect_rej[key]))
            for reason in mutect_rej[key]:
                #rcounts[reason] = rcounts.get(reason, 0) + 1/float(len(mutect_rej[key]))
                rcounts[reason] = rcounts.get(reason, 0) + 1
                s += 1
                
        print "#rejection category\tfreq"
        for (k, v) in sorted(rcounts.items(), key = lambda x: x[1]):
            #print k, v/float(s)
            print "%s\t%.6f" % (k, v/float(len(lofreq_uniq_vars)))
            
        


    
if __name__ == "__main__":
    main()
    #LOG.info("Successful program exit")

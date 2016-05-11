#!/bin/env python3
import os
import sys
import bisect
from collections import namedtuple
from collections import OrderedDict

Region = namedtuple('Region', ['chrom', 'start', 'end', 'len'])


# Largest non-N containing region is ~4% of total length (without blacklisted chromosomes):
#Region(chr='4', start=75452279, end=191044276, len=115591997)
MAX_SIZE = 0.04


def region_len_sum(rs):
    if len(rs) == 0:
        return 0
    else:
        return sum([r.len for r in rs])


def parse_fai(fai):
    """parse fai file and return chromosome and length as dict
    """
    chroms = OrderedDict()
    with open(fai) as fh:
        for line in fh:
            ls = line.rstrip().split()
            assert len(ls) >= 5
            c, l = ls[0], int(ls[1])
            assert c not in chroms
            chroms[c] = l
    return chroms


def parse_bed(bed):
    """parse bed file and return parse regions as dict with chroms as key
    and region namedtuples as list
    """
    regions = OrderedDict()
    with open(bed) as fh:
        for line in fh:
            ls = line.rstrip().split()
            assert len(ls) >= 3
            chrom, start, end = ls[0], int(ls[-2]), int(ls[-1])
            assert end > start
            reg = Region._make([chrom, start, end, end-start])
            if reg.chrom not in regions:
                regions[reg.chrom] = []
            regions[reg.chrom].append(reg)
    return regions


def main(bed, fai, outprefix):
    """main function
    """

    full_chrom_lens = parse_fai(fai)
    print("Parsed {} full chroms from {}".format(len(full_chrom_lens), fai))

    regions = parse_bed(bed)
    print("Parsed {} regions from {}".format(len(regions), bed))
    for chrom in regions:
        assert chrom in full_chrom_lens

    blacklist = ['NC_007605', 'hs37d5', 'hs37d5']
    blacklist.extend([r for r in regions if r.startswith("GL")])
    print("Blacklisted chromosomes: {}".format(blacklist))

    print("Deleting regions in blacklisted chromosomes")
    for chrom in blacklist:
        if chrom in regions:
            del regions[chrom]

    total_non_n_size = 0
    for chrom, reg_list in regions.items():
        total_non_n_size += sum([r.end-r.start for r in reg_list])
    print("Total non-N length (without blacklisted regions) is {}".format(total_non_n_size))


    non_n_size = dict()
    for chrom, reg_list in regions.items():
        non_n_size[chrom] = region_len_sum(reg_list)
        print("Non_n_size[{}]\t=\t{}".format(chrom, non_n_size[chrom]))


    out_bed_ctr = 0
    nosplit_list = []
    for chrom, reg_list in regions.items():

        # first the ones that need no split
        chrom_non_n_perc = non_n_size[chrom]/float(total_non_n_size)
        if chrom_non_n_perc <= MAX_SIZE:
            nosplit_list.append(chrom)
            print("No need to split {} ({:.2f}% of non-N regions)".format(chrom, chrom_non_n_perc))
            with open(outprefix + ".{}.bed".format(out_bed_ctr), 'w') as fh:
                fh.write("{}\t{}\t{}\n".format(
                    chrom, 0, full_chrom_lens[chrom]))
            out_bed_ctr += 1

        else:

            # then divide the rest into to regions such that the non-N split is even
            # and divide at N region
            running_perc = []
            cur_len = 0
            for reg in reg_list:
                cur_len += reg.len
                running_perc.append(cur_len/float(non_n_size[chrom]))
            assert cur_len == non_n_size[chrom]
            split_at = bisect.bisect_left(running_perc, 0.5)# split at 50%
            # since we split at 50% we only need two bed entries:
            # split before index split_at
            with open(outprefix + ".{}.bed".format(out_bed_ctr), 'w') as fh:
                fh.write("{}\t{}\t{}\n".format(
                    chrom, 0, reg_list[split_at-1].end))
            out_bed_ctr += 1
            with open(outprefix + ".{}.bed".format(out_bed_ctr), 'w') as fh:
                fh.write("{}\t{}\t{}\n".format(
                    chrom, reg_list[split_at].start, full_chrom_lens[chrom]))
            out_bed_ctr += 1


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} bed fai outprefix\n".format(os.path.basename(sys.argv[0])))
        sys.exit(1)
    BED, FAI, OUTPREFIX = sys.argv[1:]
    for f in [BED, FAI]:
        assert os.path.exists(f)
    main(BED, FAI, OUTPREFIX)

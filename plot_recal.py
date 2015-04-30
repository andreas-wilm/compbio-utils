#!/usr/bin/env python

import os
import sys

from collections import namedtuple

import matplotlib.pyplot as plt
import numpy
    
ObsEmpQualData = namedtuple('ObsEmpQualData', ['ReadGroup', 'QualityScore', 'EventType', 'EmpiricalQuality', 'Observations', 'Errors'])


def parsed_obs_emp_qual(recal_table):
    """tested with gatk 2.7-4-g6f46d11
    """

    parsed_data = []
    with open(recal_table) as fh:
        in_data = False
        for line in fh:
            line_split = line.split()
            if in_data:
                if line == '\n':
                   break
                parsed_data.append(ObsEmpQualData._make(line_split))
            elif line_split == list(ObsEmpQualData._fields):
                in_data = True
    
    return parsed_data


def plot(recal_before, recal_after, title, outplot):

    fig = plt.figure(figsize=(12, 6))
    fig.suptitle(title)

    nplot = 0
    for (recal_title, recal_data) in [("Before", recal_before), ("After", recal_after)]:
        nplot += 1
        ax = fig.add_subplot(1, 2, nplot)

        for (event_type, color) in [("I", "red"), ("D", "blue")]:
            (obsq, empq, nobs) = zip(*[(int(x.QualityScore), float(x.EmpiricalQuality), int(x.Observations)) 
                                       for x in recal_data if x.EventType==event_type])
            area = numpy.log2(nobs)
            area = numpy.pi  * (15 * area/max(area))**2# 0 to 15 point radiuses as in scatter_demo.py
            ax.scatter(obsq, empq, s=area, c=color, label=event_type, alpha=0.5)
        
        (xlim, ylim) = ax.get_xlim(), ax.get_ylim()
        plot_lim = (0, max(max(ylim), max(xlim)))
        ax.set_xlim(plot_lim)
        ax.set_ylim(plot_lim)
        diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
        
        loc = 2 # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend
        ax.legend(loc=loc)
        
        ax.set_title(recal_title)

        #ax.set_xlabel("Observed")
        ax.set_xlabel("Reported")
        ax.set_ylabel("Empirical")

    plt.savefig(outplot)
    

if __name__ == "__main__":

    try:
        table_before = sys.argv[1]
        table_after = sys.argv[2]
        title = sys.argv[3]
        outplot = sys.argv[4]
    except IndexError:
        sys.stderr.write("args: recal_before.table recal_after.table title plot\n")
        sys.exit(1)

    for f in [table_before, table_after]:
        assert os.path.exists(f)
    assert not os.path.exists(outplot)

    table_data_before = parsed_obs_emp_qual(table_before)
    #print "Before\n", before_recal
    tabe_data_after = parsed_obs_emp_qual(table_after)
    #print "After\n", after_recal

    plot(table_data_before, tabe_data_after, title, outplot)



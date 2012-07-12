#!/bin/bash

# quick and dirty binning of data (column 1!)

#
# window/bin size
win=$1
# optional filename. default is stdin
file=$2

test -n "$1" || exit 1


awk -v win=$win '
BEGIN {max=-2147483648; min=max*-1}
{
  binno=int($1/win); bin[binno]+=1;
  if (binno<min) {min=binno};
  if (binno>max) {max=binno};
  printf "DEBUG in=%f binno=%f min=%f max=%f\n", $1, binno, min, max
}
END {
  print "# from", "to", "count"
  for (binno=min; binno<=max; binno+=1) {
    if (bin[binno]>0) {print binno*win, binno*win+win-1, bin[binno]}
  }
}' $2

#!/bin/bash
# make fasta sequence names unique by appending a numerical suffix that's incremented when needed
# test success with grep '^>' | sort  | uniq -d
awk '{if (/^>/) {id=$0; newid=id; i=1; while (d[newid]>0) {i+=1; newid=id "-" i;} id=newid; print id; d[id]+=1} else {print}}' $@


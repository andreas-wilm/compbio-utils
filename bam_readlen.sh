#!/bin/bash

# length of first read only
samtools view $1 | head -n1 | cut -f 10 | tr -d '\n' | wc -c


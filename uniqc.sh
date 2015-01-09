#!/bin/sh

# print a word/line histogram using awk
# original from http://kuscsik.blogspot.com/2008/02/how-to-create-histogram-using-awk.html
# alternative to sort | uniq -c, e.g. if sorting takes to long
awk ' NF > 0{ counts[$0] = counts[$0] + 1; } END { for (word in counts) print word, counts[word]; }'

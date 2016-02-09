from collections import OrderedDict
import sys


infile = open("/dev/stdin", 'r')
outfile = open("/dev/stdout", 'w')

hits = OrderedDict()
for line in infile:
    id = line.strip().split("\t")[8]
    if id not in hits:
        hits[id] = [line]
    else:
        hits[id].append(line)

for id in hits:
    if len(hits[id]) >= 2:
        for item in hits[id]:
            outfile.write(item)



infile.close()
outfile.close()
#! /usr/bin/env python

import sys
import random

in_file = sys.argv[1]

num_lines = 0
with open(in_file, 'r') as f:
    for line in f:
        num_lines += 1

# Percentage of lines to keep
pct = float(sys.argv[2])
out_file = sys.argv[3]

lines = range(num_lines)
target_num = int(num_lines * pct / 100)
print('Sampling %d of %d lines' % (target_num, num_lines))
random.shuffle(lines)
include = [0] * num_lines
for i in range(target_num):
    include[lines[i]] = 1

with open(in_file, 'r') as f_in:
    with open(out_file, 'w') as f_out:
        line_num = 0
        for line in f_in:
            if include[line_num]:
                f_out.write(line)
            line_num += 1

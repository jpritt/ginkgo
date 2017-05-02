#! /usr/bin/env python

import sys
import random

in_file = sys.argv[1]
num_lines = 222847

# Percentage of lines to keep
pct = float(sys.argv[2])
out_file = sys.argv[3]

lines = range(num_lines)
target_num = int(num_lines * pct / 100)
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

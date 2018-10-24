import sys

CN_col = int(sys.argv[3])
with open(sys.argv[1], 'r') as f_in:
  with open(sys.argv[2], 'w') as f_out:
    last_line = None
    for line in f_in:
      new_line = line.rstrip().split('\t')
      if last_line and not (last_line[0] == new_line[0] and last_line[CN_col] == new_line[CN_col]):
        f_out.write('\t'.join(last_line[:2]) + '\t' + last_line[CN_col] + '\n')
      last_line = new_line
    f_out.write('\t'.join(last_line[:2]) + '\t' + last_line[CN_col] + '\n')

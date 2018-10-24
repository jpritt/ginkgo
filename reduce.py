import sys

with open(sys.argv[1], 'r') as f_in:
  with open(sys.argv[2], 'w') as f_out:
    last_chrom = None
    last_pos = 0
    last_cn = 0
    for line in f_in:
      row = line.rstrip().split('\t')
      if last_chrom and (not(last_chrom == row[0]) or not(row[2] == last_cn)):
        f_out.write(last_chrom + '\t' + last_pos + '\t' + last_cn + '\n')
      last_chrom, last_pos, last_cn = row[0], row[1], row[2]
    f_out.write(last_chrom + '\t' + last_pos + '\t' + last_cn + '\n')

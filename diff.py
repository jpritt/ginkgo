import sys

with open(sys.argv[1], 'r') as f1:
  with open(sys.argv[2], 'r') as f2:
    f1.readline()
    f2.readline()
    line_id = 2
    for line1 in f1:
      line2 = f2.readline()
      row1 = line1.rstrip().split('\t')
      row2 = line2.rstrip().split('\t')
      if not row1[2] == row2[2]:
        print(line_id)
      line_id += 1

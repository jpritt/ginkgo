total = 0
count = 0

filename = 'skbr3_counts.txt'
with open(filename, 'r') as f:
  for line in f:
    row = line.rstrip().split(' ')
    total += int(row[0])
    count += 1

print(total / float(count))

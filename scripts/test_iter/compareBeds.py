#! /usr/bin/env python

import sys

fileA = sys.argv[1]
formatA = int(sys.argv[2]) # 0 = Chr End Depth, 1 = Chr Start End
fileB = sys.argv[3]
formatB = int(sys.argv[4])
lengths = sys.argv[5]

def chromOffsets(f):
  offsets = dict()
  curr_offset = 0
  for line in f:
    row = line.rstrip().split(' ')
    offsets[row[0]] = curr_offset
    curr_offset += int(row[1])

  return offsets

def readBedWithDepth(f, offsets):
  bed = []
  lastPos = 0
  for line in f:
    row = line.rstrip().split('\t')
    pos = int(row[1]) + offsets[row[0]]
    bed.append((lastPos, pos, int(row[2])))
    if lastPos > pos:
      print(line)
      exit()
    lastPos = pos

  return bed

def readBedNoDepth(f, offsets):
  bed = []
  lastEnd = 0
  for line in f:
    row = line.rstrip().split('\t')
    start = int(row[1]) + offsets[row[0]]
    end = int(row[2]) + offsets[row[0]]

    if len(bed) == 0:
      if not start == lastEnd:
        bed.append((lastEnd, start, 0))
      bed.append((start, end, 1))
    else:
      last = bed[-1]
      if start == last[0] and end == last[1]:
        bed[-1] = (last[0], last[1], last[2]+1)
      else:
        if not start == lastEnd:
          bed.append((lastEnd, start, 0))
        bed.append((start, end, 1))
    lastEnd = end

  max_id = offsets['chrY'] + 59373566
  if not lastEnd == max_id:
    bed.append((lastEnd, max_id, 0))

  return bed
  
with open(lengths, 'r') as f:
  chromOffsets = chromOffsets(f)

bedA, bedB = [], []
with open(fileA, 'r') as f:
  header = f.readline()
  if formatA:
    bedA = readBedNoDepth(f, chromOffsets)
  else:
    bedA = readBedWithDepth(f, chromOffsets)

with open(fileB, 'r') as f:
  header = f.readline()
  if formatB:
    bedB = readBedNoDepth(f, chromOffsets)
  else:
    bedB = readBedWithDepth(f, chromOffsets)

for i in range(1, len(bedA)):
  if bedA[i][1] < bedA[i][0]:
    print('Error! First bed file chroms not in same order as chromosomes list!')
    print(i)
    print(bedA[i-1:i+2])
    exit()
for i in range(1, len(bedB)):
  if bedB[i][1] < bedB[i][0]:
    print('Error! Second bed file chroms not in same order as chromosomes list!')
    exit()

startA, endA, valA = bedA[0][0], bedA[0][1], bedA[0][2]
startB, endB, valB = bedB[0][0], bedB[0][1], bedB[0][2]
indexA, indexB = 0,0
numA = len(bedA)
numB = len(bedB)

correct = 0
incorrect = 0

while indexA < numA or indexB < numB:
  if indexA >= numA:
    incorrect += (endB - startB)
    startB = endB
  elif indexB >= numB:
    incorrect += (endA - startA)
    startA = endA
  elif not startA == startB:
    newStart = max(startA, startB)
    incorrect += newStart - min(startA, startB)
    startA, startB = newStart, newStart
  else:
    newEnd = min(endA, endB)
    if valA == valB:
      if newEnd < startA:
        print('Negative')
        print('  %d - %d: %d' % (startA, endA, valA))
        print('  %d - %d: %d' % (startB, endB, valB))
      correct += (newEnd - startA)
    else:
      incorrect += (newEnd - startA)

    startA, startB = newEnd, newEnd

  if startA == endA:
    indexA += 1
    if indexA < numA:
      startA, endA, valA = bedA[indexA][0], bedA[indexA][1], bedA[indexA][2]
  if startB == endB:
    indexB += 1
    if indexB < numB:
      startB, endB, valB = bedB[indexB][0], bedB[indexB][1], bedB[indexB][2]

#print('%f %% match (%d bases correct, %d incorrect)' % (float(correct) / (correct+incorrect), correct, incorrect))
#print('%s\t%f' % (sys.argv[6], float(correct) / (correct+incorrect)))
print(float(correct) / (correct+incorrect))

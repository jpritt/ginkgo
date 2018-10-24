
genomes = dict()
with open('dir_list.txt', 'r') as f:
  for line in f:
    with open('uploads/' + line.rstrip(), 'r') as config:
      for l in config:
        l = l.rstrip()
        if l[:14] == 'chosen_genome=':
          g = l[14:]
          if g[0] == "'":
            g = g[1:-1]
          if g in genomes:
            genomes[g] += 1
          else:
            genomes[g] = 1

g_counts = [(v,k) for k,v in genomes.items()]
g_counts.sort(reverse=True)
for g in g_counts:
  print(g[1] + '\t' + str(g[0]))


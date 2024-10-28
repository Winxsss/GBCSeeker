import sys
Gene = {}
with open(sys.argv[1], 'r') as f1:
    for line in f1:
        if line.startswith('Sample'):
            next
        s,g,af = line.strip().split('\t')
        Gene[g] = af
with open(sys.argv[2], 'r') as f2:
    for line in f2:
        line = line.strip()
        g = line.strip()
        if g in Gene:
            print(f'{g}\t{Gene[g]}')
        else:
            print(f'{g}\t0')

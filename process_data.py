from collections import defaultdict
import csv
from pyclts import CLTS
clts = CLTS('clts-2.2.0')

def cltsconvert(s):
    return(s.replace('͡','').replace('ɡ','g'))

values = list(csv.reader(open('values.csv','r'),delimiter=','))

languages = list(csv.reader(open('languages.csv','r'),delimiter=','))

languages = [l for l in languages if l[8] == 'Austronesian']

langnames = [l[0] for l in languages]

values = [l for l in values if l[1] in langnames]

for i,l in enumerate(values):
    values[i][3] = cltsconvert(values[i][3])

values = [l for l in values if l[3] in clts.bipa.sounds.keys()]

segs = sorted(set([l[3] for l in values]))

D = len(segs)

langvalues = defaultdict(list)

for l in values:
    langvalues[l[1]].append(l[3])

char_data = [['lang']+segs]

for lang in langnames:
    vals = ['0']*D
    for s in langvalues[lang]:
        vals[segs.index(s)] = '1'
    char_data.append([lang]+vals)

segfeats = {s:clts.bipa[cltsconvert(s)].name.split() for s in segs}

featlist = sorted(set([f for v in segfeats.values() for f in v]))

J = len(featlist)

seg_info = [['segment']+featlist]
for s in segs:
    vals = ['0']*J
    for t in featlist:
        if t in segfeats[s]:
            vals[featlist.index(t)] = '1'
    seg_info.append([s]+vals)

f = open('segmental_features.tsv','w')
for l in seg_info:
    print('\t'.join(l),file=f)

f.close()

f = open('character_data.tsv','w')
for l in char_data:
    print('\t'.join(l),file=f)

f.close()
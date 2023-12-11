import pandas as pd
selected = set(s.strip() for s in open('selected.txt'))
interesting = []
for ch in pd.read_table('./GMGC10.emapper2.annotations.topk10k.tsv.gz', index_col=0, header=None, comment='#', chunksize=100):
    ch = ch[ch.index.map(selected.__contains__)]
    ch = ch.copy()
    interesting.append(ch)
    
interesting = pd.concat(interesting)


# Tips &amp; tricks for large-scale data processing in Python


---

## 1: Use [Jug](https://jug.rtfd.io/) (_Luis_)


```python
import pandas as pd
from glob import glob

from jug import TaskGenerator

@TaskGenerator
def process(fna, ofile):
    from collections import Counter
    from fasta import fasta_iter

    r = []
    hs = []
    for h,seq in fasta_iter(fna):
        r.append(pd.Series(Counter(seq)))
        hs.append(h)
    r = pd.DataFrame(r, index=hs)
    r.to_csv(ofile, sep='\t')
    return ofile

ifiles = list(sorted(glob('demo-data/*.fna.gz')))

partials = {}
for i, fna in enumerate(ifiles):
    ofile = f'outputs/chunk{i:02}.tsv'
    ofile = process(fna, ofile)
    partials[fna] = gc_fraction(ofile)
```

---

## 2: Always use compression (_Luis_)

## BAD

```python
with open('output.txt', 'wt') as out:
    ...
```

## GOOD

```python
import gzip
with gzip.open('output.txt.gz', 'wt', compresslevel=0) as out:
    ...
```

## WHY

1. Saves some diskspace (obvious, but maybe least important)
2. Provides builtin error checking!
3. If you use `compresslevel=0`, there is almost no computational cost

---

## 3: YOUR TIP (_YOUR NAME_)


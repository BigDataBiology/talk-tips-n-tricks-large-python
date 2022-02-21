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

# 3: Change traditional thinking (_Yiqian_)
## Question purpose

**QUERY** dataset contains millions of smORFs.

Many **dbin** dataset contains very short peptides(file size diverses from 1M to 500M).

We want to find exact match of pepetides against smORFs.


## Works very well for small datasets, but becomes too slow at scale

```python
def dbsearch(dbin,query):
    for ID,seq in fasta_iter(query):
        for hash_string in dbin:
            if (len(seq) >= len(hash_string)) and (hash_string in qstr):
                ...
```

### Example:
**QUERY**:
ABCDEFGHIJKLMN,
OPQRSTUVWXYZ
...

**dbin**:
CDEF,
XYZ,
IJKLMN
...

Two layer `for` cycle to find：

If CDEF is a substring of ABCDEFGHIJKLMN and OPQRSTUVWXYZ ... 

IF XYZ is a substring of ABCDEFGHIJKLMN and OPQRSTUVWXYZ ... 
...

### Problem:
If **QUERY** has n sequences,**dbin** has m sequences.

Time complexity is O(n×m).

Running time depends on how big the `QUERY` dataset is and the `dbin` dataset is. 

It will take several minutes when `dbin` dataset is 1Mb.But it will take too long to run when `dbin` dataset is 500Mb. 

## Work well for big datasets

```python
def dbsearch(dbin,query,db_min,db_max):
    for ID,seq in fasta_iter(query):
        for i in range(db_min,db_max+1):
            if len(seq) >= i:
                for j in range(0,len(seq)-i+1):
                    qstr = seq[j:j+i]
                    if qstr in dbin:
                       ...
```

### Example:
**QUERY**:
ABCDEFGHIJKLMN,
OPQRSTUVWXYZ
...

**dbin**:
CDEF,
XYZ,
IJKLMN
...

Loop each sequence in `QUERY` and split it into substrings according to min length and max length of `dbin` sequences.

Eg.
min length of `dbin` = 3,max length of `dbin` = 6

`QUERY` seq = ABCDEFGHIJKLMN,split it into:

BCD CDE DEF EFG FGH GHI HIJ IJK JKL KLM LMN

ABCD BCDE CDEF...

ABCDE BCDEF...

ABCDEF BCDEFG...

Then find if each substring is in `dbin` `set`. If so,it means this substring(the same peptides) in `dbin` is a substring of this sequence in `QUERY`.

### Solve problem
If **QUERY** has n sequences,**dbin** has m sequences.

Time complexity is like O(n)? Because seaching in a `set` is O(1)

Running time only depends on how big the 'QUERY' dataset is.

---

## 4: Use datatable to read large table into memory (_Hui_)

## BAD

```python
import pandas as pd
df = pd.read_table('input.tsv', sep='\t')
```

## GOOD

```python
import datatable as dt
df = dt.fread('input.tsv', sep='\t')
```

## WHY

1. `datatable` uses multi-threading for file reading.
2. `datatable` provides `to_pandas` API so that you can do downstream analysis using pandas.
3.  Depends on how many cores your machine has, this could save you 50 ~ 70% of time.

## 5 Use Progress bar to see whether your loop is becoming slower and slower

```python
from tqdm import tqdm, trange


# use trange instead of range, you get a progress bar
for i in trange(10**6):
	...


# add tqdm(), you get a progress bar
for item in tqdm(items):
	...
```

## WHY

1. The progress bar looks like this: `76%|████████████████████████        | 7568/10000 [00:33<00:10, 229.00it/s]`
2. You're the first user of your script. If you have experiences strugling to wait for hours/days to get a feedback from your script, try this.
3. Sometimes your loop is becoming slower and slower as iteration and you don't know. With the progress bar you'll see it.

## 6: YOUR TIP (_YOUR NAME_)


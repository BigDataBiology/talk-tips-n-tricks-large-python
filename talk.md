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

## 3: Consider time complexity (_Yiqian_)

**QUERY** dataset contains millions of smORFs. Many **dbin** dataset contains very short peptides(file size diverses from 1M to 500M).

We want to find exact match of pepetides against smORFs.

## Works very well for small datasets, but becomes too slow at scale

```python
def dbsearch(dbin,query):
    for ID,seq in fasta_iter(query):
        for hash_string in dbin:
            if (len(seq) >= len(hash_string)) and (hash_string in qstr):
                ...
```

---

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

---

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

---

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

Loop each sequence in `QUERY` and split it into substrings according to min length and max length of `dbin` sequences
(_eg._, min length of `dbin` = 3,max length of `dbin` = 6)

`QUERY` seq = ABCDEFGHIJKLMN,split it into:

BCD CDE DEF EFG FGH GHI HIJ IJK JKL KLM LMN

ABCD BCDE CDEF...

ABCDE BCDEF...

ABCDEF BCDEFG...

Then find if each substring is in `dbin` `set`. If so,it means this substring(the same peptides) in `dbin` is a substring of this sequence in `QUERY`.

---

### Solve problem

If **QUERY** has n sequences,**dbin** has m sequences.

Time complexity is like O(n)? Because seaching in a `set` is O(1)

Running time only depends on how big the 'QUERY' dataset is.

---

## Basic costs in Python

```python
if elem in container:
   ...
```

is `O(N)` if `container` is a list, but `O(1)` if `container` is a
`set`/`dict`

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
3. Depends on how many cores your machine has, this could save you 50 ~ 70% of time.

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

---

## 6: Sparse matrix to reduce the memory usage (_Shaojun_)

## Generating the graph

```python
embedding_matrix = kneighbors_graph(...)
kmer_matrix = kneighbors_graph(...)
kmer_matrix.eliminate_zeros()
kmer_matrix.data.fill(1.)
embedding_matrix = embedding_matrix.multiply(kmer_matrix)
embedding_matrix.eliminate_zeros()
```

## Finding the edges in the graph

```python
embedding_matrix.data[embedding_matrix.data <= 1e-6] = 0
X, Y, V = sparse.find(embedding_matrix)
above_diag = Y > X
X = X[above_diag]
Y = Y[above_diag]
edges = [(x,y) for x,y in zip(X, Y)]
edge_weights = V[above_diag]
```

---

## 7: numexpr (_Shaojun_)

```python
a = np.random.rand(1e6)
b = np.random.rand(1e6)
timeit 2*a + 3*b
timeit ne.evaluate("2*a + 3*b")
```

```python
res = ne.evaluate('(log(v1) - log(v2))/2 + ( (m1 - m2)**2 + v2 ) / ( 2 * v1 ) - half',
                    {
                     'v1': v1,
                     'v2': v2,
                     'm1': m1,
                     'm2': m2,
                     'half': np.float32(0.5),
                     })
```

---

## 8: tempfile.TemporaryDirectory() (_Shaojun_)

```python
with tempfile.TemporaryDirectory() as tdir:
		output = os.path.join(tdir, 'output.fasta')
		...
```

## Set the output path of the temporary directory

```python
with tempfile.TemporaryDirectory(tdir = tmpdir) as tdir:
		...
```

or

```python
os.environ['TMPDIR'] = tmpdir
```

---

## 9: Understand a regular expression (_Anna_)

Regular expressions (regex) are characters that define a search pattern. They are used for string manipulation, such as extracting or modifying text.
I use https://regex101.com/ to understand and test easily the regular expression with each particular case.

An example, extracting the mOTUs code from the mOTUs results table:

```python
print(motus.head())

index                                               SAMEA4689043  ...  SAMEA4546741
Abiotrophia defectiva [ref_mOTU_v25_04788]                    40  ...            22
Absiella dolichum [ref_mOTU_v25_03694]                       142  ...             8
Acaricomes phytoseiuli [ref_mOTU_v25_06702]                    1  ...           287

motus['motus_code'] = motus['index'].str.extract('\[(.*?)\]')
```

---

## 10: Process by chunks in Pandas (_Luis_)

If you need to process very large tables with Pandas, they might not fit in
memory comfortably, but often you can process them by chunks using
[`pandas.read_csv`](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html)
and the `chunksize` argument

```python
for chunk in pd.read_csv('very-large-table.xz', \
                sep='\t', \
                index_col=0, \
                header=None, \
                names=['col1','col2'], \
                chunksize=1_000_000): \
    chunk = chunk.query('col1 == 1')
    ...
```

---

## 11: Generate index file and use .npy file to quickly select data (_Yiqian_)

Original file

```
GMSC10.100AA.000_000_000    d__Archaea
GMSC10.100AA.000_000_001    Unknown
GMSC10.100AA.000_000_002    d__Archaea;p__Aenigmatarchaeota;c__Aenigmatarchaeia
```

Index file

```
0   Unknown
1   d__Archaea
2   d__Archaea;p__Aenigmatarchaeota;c__Aenigmatarchaeia
3   d__Archaea;p__Aenigmatarchaeota;c__Aenigmatarchaeia;o__Aenigmatarchaeales;f__Aenigmatarchaeaceae;g__Aenigmatarchaeum;s__Aenigmatarchaeum subterraneum
4   d__Archaea;p__Aenigmatarchaeota;c__Aenigmatarchaeia;o__CG10238-14
5   d__Archaea;p__Aenigmatarchaeota;c__Aenigmatarchaeia;o__CG10238-14;f__CG10238-14
```

Index formatted file `GMSC10.100AA.habitat.idx.tsv`

```
GMSC10.100AA.000_000_000    1
GMSC10.100AA.000_000_001    0
GMSC10.100AA.000_000_002    2
```

Generate .npy file

```python
import numpy as np
df = np.loadtxt('GMSC10.100AA.habitat.idx.tsv',dtype=int, usecols=(1))
np.save('GMSC10.100AA.habitat.npy',df)
```

Use .npy file

```
import numpy as np
tax = np.load('GMSC10.100AA.habitat.npy',mmap_mode='r')
tax[2]
```

---

## 12: Generator using yield (_Shaojun_)

'yield' provides an efficient way to handle iteration and data generation, especially suitable for processing large-scale or dynamically generated datasets.

For example, using 'yield' to iterator the fasta file:

```python
def fasta_iter(fname, full_header=False):
    '''Iterate over a (possibly gzipped) FASTA file

    Parameters
    ----------
    fname : str
        Filename.
            If it ends with .gz, gzip format is assumed
            If .bz2 then bzip2 format is assumed
            if .xz, then lzma format is assumerd
    full_header : boolean (optional)
        If True, yields the full header. Otherwise (the default), only the
        first word

    Yields
    ------
    (h,seq): tuple of (str, str)
    '''
    header = None
    chunks = []
    if hasattr(fname, 'readline'):
        op = lambda f,_ : f
    elif fname.endswith('.gz'):
        import gzip
        op = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        op = bz2.open
    elif fname.endswith('.xz'):
        import lzma
        op = lzma.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

for h, seq in fasta_iter(fasta_path):
	...
```

---

## 13: pysam for manipulation of SAM and BAM files (_Vedanth_)

[`pysam`](https://pysam.readthedocs.io/en/latest/api.html) provides functionalities to read, write and manipulate SAM and BAM files.

For example,

```python
import pysam

# Read SAM/BAM file
samfile = pysam.AlignmentFile("alignment.bam", "rb")

# Iterate over alignments
for read in samfile.fetch():
    # Access alignment information
    print(read.query_name, read.reference_name, read.reference_start)

```

---

## 14: Extracting specific contigs from a metagenome assembly (_JP_)

````python
from Bio import SeqIO
import sys
fasta_file = sys.argv[1] #fasta file
number_file = sys.argv[2] #ID file

wanted = []
handle = open(number_file)
for line in handle:
    line = line.strip()
    if line != "":
        wanted.append(line)
wanted = set(wanted)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False
for seq in fasta_sequences:
    if seq.id in wanted:
        print (">" + seq.id +"\n" + seq.seq)```

- Save as .py file
- Requires Biopython to run
- Needs two input files: an assembly fasta file and a text file with the contig IDs of interest
- Run as follows:

```bash
python file_name.py input.1.fastafile.fasta input.file.2.with.contigs.ofinterest.txt
````

---

## 15: Polars is a much faster/memory efficient alternative to Pandas (_Sebastian_)

Polars is a Rust-based DataFrame library that is much faster and more memory efficient than Pandas. It is multithreaded out of the box, therefore be careful with parallization as it will most often not lead to better performance.

The syntax is little different but easy to get used to. If you have worked with tidyverse in R, you will find it very similar.

Polars has an elaborate Python API documentation that you can find [here](https://docs.pola.rs/py-polars/html/reference/index.html).

```python
# Initialize max threads for polars
os.environ['POLARS_MAX_THREADS'] = '1'
import polars as pl

df = pl.read_csv('data.csv')

df = df.filter( # filter rows in df
        pl.col("column_name") == "value"
    )\
    .with_columns( # modify and create new columns
        pl.col("column_name").cast(pl.Int64), # change column type
        (pl.col("column_name") * 2).alias("new_column_name"), # create new column
        third_column = pl.col("column_name") + pl.col("new_column_name") # create new column
    )

# conditional modification
forward_indices = pl.Series([1, 2, 3, 4, 5])
motif = "GATC"

## If position is in forward_indices and strand is "+", create
df = df.with_columns(
        pl.lit("").alias("motif"), # create empty string column
        pl.when(                   # conditionally modify the column
            (pl.col("position").is_in(forward_indices)) & (pl.col("strand") == "+")
        ).then(pl.lit(motif)).otherwise(pl.col("motif")).alias("motif")
    )

df.write_csv('output.csv')
```

It takes some getting use to but once you get the hang of it, you will find it much faster.

---

## 16: [Snakemake](https://snakemake.readthedocs.io/en/stable/) to create automated analysis workflow  (___Rizky___)

Snakemake is a workflow management system to create reproducible and scalable data analysis workflow using Python-based language.
It does take time (and patience) to build it but the resulting files (.smk, .yaml, profile) can be shared pretty easily.

Example (part of my minimap_danica.smk file),
```python
import os

###data, scripts and database directory###
parent_dir = os.path.dirname(workflow.basedir)
db_dir = os.path.join(parent_dir, "0_database")
data_dir = os.path.join(parent_dir, "0_data")

#mapping the data to the user defined reference database using minimap2
rule map_minimap2:
    input:
    #F and R should have the same name with R1 suffix for forward and R2 for reverse located in 0_data directory.
    #example: data1.fastq.gz and data2.fastq.gz and the rule should read it automatically and use the filename (accession name) to ensure reproducibility
        F = f"{data_dir}/cami/{{filename}}_R1.fastq.gz",
        R = f"{data_dir}/cami/{{filename}}_R2.fastq.gz",
    #ref is the indexed database from user
        ref = f"{db_dir}/danica.db"
    output:
        "1_assign_danica2/{filename}.sam"
	#the directory would be created by snakemake automatically

```

---


## 16: [Snakemake](https://snakemake.readthedocs.io/en/stable/) part 3  (___Rizky___)

```python
    conda:
        "envs/minimap2.yaml"
	#user can create minimap2.yaml to list all dependencies and snakemake will make its own environment based on the yaml file
    threads:
        1
	#both threads and resources can be user defined per rule, needed for certain profile
    resources:
        mem_mb = 16000,
        runtime = "96h"
    benchmark:
        "0_logs_danica2/{filename}/map_minimap2.benchmark.txt"
	#provide time and memory usage per rule
    log:
        "0_logs_danica2/{filename}/map_minimap2.log"
        #log and benchmark files are located in 0_logs folder and the subdirectory with read accession name, it will put STDERR
    shell:
        "(minimap2 -ax sr -t {threads} {input.ref} {input.F} {input.R} > {output}) 2> {log}"
	#shell command to be executed in the rule
```

---


## 16: [Snakemake](https://snakemake.readthedocs.io/en/stable/) part 3  (___Rizky___)

Checking if the rule is okay
```bash
snakemake --lint --snakefile minimap_danica.smk
```
Dry-run to see the flow
```bash
snakemake -np --snakefile minimap_danica.smk
```
General execution is snakemake and then the intended output
```bash
snakemake CAMI_1.sam
```
Execution and assigning thread usage -> 10 threads
```bash
snakemake --cores 10 CAMI_1.sam
```
Using profile (this one is snakemake adaptor to Lyra - QUT supercomputer)
```bash
snakemake --profile mqsub-lyra-v8 --snakefile minimap_danica.smk
```
After successful execution, Snakemake will write-protect the output file in the filesystem, so that it can’t be overwritten or deleted.

---


## 17: [Numba](https://numba.pydata.org/) for just in time compiling in Python (_Vedanth_)

- Numba is a just-in-time compiler for Python.
- Works best on code that uses NumPy arrays and functions, and loops.
- Most common way to use Numba is through its decorators that can be applied to functions to instruct Numba to compile them.

Example:
```bash
from numba import njit
import random

@njit
def monte_carlo_pi_calculation(nsamples):
    acc = 0
    for i in range(nsamples):
        x = random.random()
        y = random.random()
        if (x ** 2 + y ** 2) < 1.0:
            acc += 1
    return 4.0 * acc / nsamples
```

- ~400ms without njit to ~6ms with njit
- ---


## 18: BioG-HGT pipeline - BioGeochemical Horizontal Gene Transfer pipeline (JP)
- https://github.com/Shaedeycool/BioG-HGT_wd
- for identification of horizontally transferred biogeochemical genes
- Identifies horizontally transferred biogeochemical genes with mobile genetic elements found on the same contig
- Output file: the BGCs, along with the MGE and corresponding MAG
- The relative abundance of the BGC gene, MGE e-value and the sequences for both the BGC gene and MGE in the contig
[
](https://github.com/Shaedeycool/BioG-HGT_wd/blob/main/workflow_1.png)![image](https://github.com/BigDataBiology/talk-tips-n-tricks-large-python/assets/124160719/fd709d4e-6b96-4894-93e8-323e03889a9d)


## 19: Use query() for more efficient filtering of a dataframe (Anna)

## Less efficient: Iteration

- Iteration over a DataFrame involves using loops to process each row or column individually (Python-level looping)
- With large datasets, iterating row by row can be particularly slow
- Does not take advantage of NumPy's operations

```python
import pandas as pd

metadata = pd.read_csv('metadata.csv')
filtered_rows = []

for index, row in metadata.iterrows():
    if row['size'] == 'large' and row['age'] != 'puppy' and row['kg'] > 25:
        filtered_rows.append(row)

metadata_filt = pd.DataFrame(filtered_rows)
```
## Less efficient: Boolean indexing

- By applying conditions to columns in the DataFrame, it creates boolean arrays. These boolean arrays are then used to select rows from the DataFrame.
- Boolean indexing uses NumPy operations for efficient filtering. But, the conditions are evaluated element-wise for each column.
- Boolean indexing allows for more flexibility in combining conditions using logical operators (&, |, ~)

```python
import pandas as pd

metadata = pd.read_csv('metadata.csv')
metadata_filt = metadata[metadata['size'] == 'large' & metadata['age'] != 'puppy' & metadata['kg'] > 25]
```

## Better: Use query()

- More concise and readable way to express complex filtering conditions: you provide a boolean expression as a string, which represents the filtering conditions.
- When you use the query() method in pandas, pandas parses the string expression and converts it into a boolean array using NumPy operations.
- These boolean arrays represent whether each row in the DataFrame satisfies the filtering conditions or not.
- Particularly advantageous when dealing with large datasets because it delegates the filtering process to the underlying NumPy operations.

```python
import pandas as pd

metadata = pd.read_csv('metadata.csv')
metadata_filt = metadata.query("size == 'large' and age != 'puppy' and kg > 25")
```

In summary, while query() tends to be the most efficient option for filtering large datasets with complex conditions, the performance difference between query(), boolean indexing, and iteration may vary based on factors such as dataset size, complexity of conditions, and the number of conditions. 

---
## 20: Basic parallelisation (Brett)
- Use concurrent futures for basic parallelisation.
- Useful for monotonous tasks (i.e.: reading through >100,000 fasta files) to speed it up.

#### Without
```
def get_species_names_list(all_cluster_files, threads):
    species_names = set()
    for file in all_cluster_files:
        species_names.update(job.result())
    return sorted(list(species_names))
```

#### With
```
def get_species_names_list(all_cluster_files, threads):
    species_names = set()
    with ThreadPoolExecutor(max_workers = threads) as executor:
        jobs = [executor.submit(get_species_name_from_file, file) for file in all_cluster_files]
        for job in jobs:
            species_names.update(job.result())
        return sorted(list(species_names))
```
---

---
## 21: Polars lazy/streaming
- Use lazy polars to optimise CPU/memory usage (e.g. skip reading column from filter if it is not used)
- Use polars streaming to analyse data that wont fit in memory
- Can sink to file when the result wont fit in memory either

```
# Create lazy query plan
q1 = (
    pl.scan_csv("docs/data/iris.csv")
    .filter(pl.col("sepal_length") > 5)
    .group_by("species")
    .agg(pl.col("sepal_width").mean())
)
# Collect output df in memory
df = q1.collect(streaming=True)
# Or sink output to file
q1.sink_csv("docs/iris_analysed.csv")
```
---

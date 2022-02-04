import pandas as pd
from glob import glob

from jug import TaskGenerator

@TaskGenerator
def process(fna, ofile):
    from collections import Counter
    from fasta import fasta_iter
    from time import sleep
    from os import makedirs
    sleep(4) # Pretend we are doing a lot of work

    makedirs('outputs', exist_ok=True)

    r = []
    hs = []
    for h,seq in fasta_iter(fna):
        r.append(pd.Series(Counter(seq)))
        hs.append(h)
    r = pd.DataFrame(r, index=hs)
    r.to_csv(ofile, sep='\t')
    return ofile

@TaskGenerator
def gc_fraction(c):
    c = pd.read_table(c, index_col=0)
    gcf = c.eval('GCf = (G+C)/(G+C+T+A)')['GCf']
    return gcf.describe()

@TaskGenerator
def build_table(partials):
    return pd.DataFrame(partials).T

@TaskGenerator
def plot_results(data):
    from matplotlib import pyplot as plt
    from matplotlib import style
    style.use('default')

    data = data.rename(index=lambda ix: ix.split('/')[1].split('.')[1], )
    data = data.sort_values(by='mean')

    fig,ax = plt.subplots()
    ax.plot(data['25%'], label='Q1')
    ax.plot(data['50%'], label='Q2')
    ax.plot(data['75%'], label='Q3')
    ax.legend(loc='best')

    ax.set_ylabel('GC fraction')
    ax.set_xlabel('Habitat')

    for x in ax.get_xticklabels():
        x.set_rotation(90)

    fig.tight_layout()
    fig.savefig('gc-per-habitat-jugfile.pdf')

ifiles = glob('demo-data/*.fna.gz')
ifiles.sort()

partials = {}
for i, fna in enumerate(ifiles):
    ofile = f'outputs/chunk{i:02}.tsv'
    ofile = process(fna, ofile)
    partials[fna] = gc_fraction(ofile)

final = build_table(partials)
plot_results(final)


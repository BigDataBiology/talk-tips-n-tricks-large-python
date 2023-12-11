import gzip

with gzip.open('output.txt.gz', 'wt', compresslevel=0) as out:
    for i in range(10000):
        out.write(f'{i}^2 = {i**2}\n')

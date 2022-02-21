from igraph import Graph
from sklearn.neighbors import kneighbors_graph
from scipy import sparse
import numpy as np

embedding = np.random.rand(1000, 100)
kmer = np.random.rand(1000, 100)

embedding_matrix = kneighbors_graph(embedding,
    n_neighbors=50,
    mode='distance',
    p=2)

kmer_matrix = kneighbors_graph(
    kmer,
    n_neighbors=50,
    mode='distance',
    p=2)

# We want to intersect the matrices, so we make kmer_matrix into a
# matrix of 1.s and multiply
kmer_matrix.eliminate_zeros()
kmer_matrix.data.fill(1.)
embedding_matrix = embedding_matrix.multiply(kmer_matrix)
embedding_matrix.eliminate_zeros()

embedding_matrix.data[embedding_matrix.data <= 1e-6] = 0
X, Y, V = sparse.find(embedding_matrix)
# Find above diagonal positions:
above_diag = Y > X
X = X[above_diag]
Y = Y[above_diag]

edges = [(x, y) for x, y in zip(X, Y)]
edge_weights = V[above_diag]

g = Graph()
g.add_vertices(1000)
g.add_edges(edges)

result = g.community_infomap(edge_weights=edge_weights, trials=1)
print(result)
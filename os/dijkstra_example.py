# define function 'xrange'

def xrange(x):

    return iter(range(x))

genes = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7']
Matrix = [[0 for x in xrange(7)] for x in xrange(7)]
Matrix[1-1][4-1] = 12
Matrix[4-1][1-1] = 12
Matrix[1-1][3-1] = 23
Matrix[3-1][1-1] = 23
Matrix[2-1][3-1] = 5
Matrix[3-1][2-1] = 5
Matrix[2-1][7-1] = 16
Matrix[7-1][2-1] = 16
Matrix[3-1][5-1] = 17
Matrix[5-1][3-1] = 17
Matrix[3-1][4-1] = 9
Matrix[4-1][3-1] = 9
Matrix[4-1][5-1] = 18
Matrix[5-1][4-1] = 18
Matrix[4-1][6-1] = 25
Matrix[6-1][4-1] = 25
Matrix[5-1][6-1] = 7
Matrix[6-1][5-1] = 7
Matrix[5-1][7-1] = 22
Matrix[7-1][5-1] = 22

# make a sparse matrix
from scipy.sparse import csr_matrix
Matrix_sparse = csr_matrix(Matrix)

# run Dijkstra's algorithm, starting at index 0
from scipy.sparse.csgraph import dijkstra
distances, predecessors = dijkstra(Matrix_sparse, indices=0, return_predecessors=True)

# print out the distance to g7
print("distance to g7=",distances[6])

# print out the path
path = []
i = 6
while i != 0:
    path.append(genes[i])
    i = predecessors[i]
path.append(genes[0])
print("path=",path[::-1])

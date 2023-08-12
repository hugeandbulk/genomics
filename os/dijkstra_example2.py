# define function 'xrange'

def xrange(x):
   
    return iter(range(x))

genes = ['g0', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11']

# create a list containing 12 lists initialised to 0
Matrix = [[0 for x in xrange(12)] for x in xrange(12)]
Matrix[0][1] = -105
Matrix[0][2] = -110
Matrix[1][3] = -132
Matrix[1][4] = -126
Matrix[3][5] = -150
Matrix[4][6] = -128
Matrix[4][7] = -166
Matrix[4][8] = -132
Matrix[4][9] = -118
Matrix[5][6] = -128
Matrix[5][7] = -166
Matrix[5][8] = -132
Matrix[5][9] = -118
Matrix[6][10] = -196
Matrix[7][8] = -132
Matrix[8][9] = -118
Matrix[9][10] = -196
Matrix[2][4] = -126
Matrix[2][5] = -150
Matrix[10][11] = -100

# make a sparse matrix
from scipy.sparse import csr_matrix
Matrix_sparse = csr_matrix(Matrix)

# run Dijkstra's algorithm, starting at index 0
from scipy.sparse.csgraph import dijkstra
distances, predecessors = dijkstra(Matrix_sparse, indices=0, return_predecessors=True, directed=True)

# print out the distance to g11
print("distance to g11=",distances[11])

# print out the path
path = []
i = 11
while i != 0:
    path.append(genes[i])
    i = predecessors[i]
path.append(genes[0])
print("path=",path[::-1])


import networkx as nx
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def generate_random_3Dgraph(n_nodes, radius, seed=None):

    if seed is not None:
        random.seed(seed)

    # Generate a dict of positions
    pos = {i: (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(n_nodes)}

    # Create random 3D network
    G = nx.random_geometric_graph(n_nodes, radius)

    return G



g = generate_random_3Dgraph(120, 3)
print(g["pos"])

# d = np.asarray(list(p.values()))


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# ax.scatter(d[:,0], d[:,1], d[:,2])

# plt.show()


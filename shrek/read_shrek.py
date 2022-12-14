import numpy as np
from matplotlib import pyplot as plt
from random import sample

from stl import mesh

def read_stl(fname):
    shrek = mesh.Mesh.from_file(fname)
    shrek = list(shrek.points.reshape(-1, 3))
    shrek = sample(shrek, 800)
    return np.array(shrek)

def read_pcd(fname):
    ret = []
    with open(fname, 'r') as rf:
        for l in rf:
            if not l[0].isnumeric():
                continue
            ret.append([float(x) for x in l.split()])
    ret = sample(ret, int(1200))
    return np.array(ret)

if __name__ == '__main__':
    fname = 'shrek.pcd'
    # shrek = read_stl(fname)
    shrek = read_pcd(fname)

    print(shrek.shape)
    np.save(fname, shrek)

    fig = plt.figure()
    prot_ax = fig.add_subplot(1, 1, 1, projection='3d')
    prot_ax.scatter(*shrek.T)
    prot_ax.set_axis_off()

    plt.show()


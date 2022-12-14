import numpy as np
from matplotlib import pyplot as plt
from random import sample

from stl import mesh

def read_stl(fname):
    shrek = mesh.Mesh.from_file(fname)
    shrek = shrek.points.reshape(-1, 3)
    shrek = sample(shrek, 800)
    return shrek

def read_pcd(fname):
    ret = []
    with open(fname, 'r') as rf:
        for l in rf:
            if not l[0].isnumeric():
                continue
            ret.append([float(x) for x in l.split()])
    ret = sample(ret, int(800))
    np.save(fname, ret)
    return np.array(ret)

if __name__ == '__main__':
    # shrek = read_pcd('./shrek.pcd')
    shrek = read_stl('shrek_wazowski.stl')

    print(shrek.shape)

    fig = plt.figure()
    prot_ax = fig.add_subplot(1, 1, 1, projection='3d')
    prot_ax.scatter(*shrek.T)

    plt.show()


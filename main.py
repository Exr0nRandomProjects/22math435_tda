# mash of examples from github/scikit-tda/ripser.py and github/scikit-tda/persim
import numpy as np
from matplotlib import pyplot as plt
from ripser import ripser
from persim import plot_diagrams
from glob import glob

# import tadasets
#
# torus = tadasets.torus(n=1200, c=2, a=1, noise=0.05)
# np.save('torus.npy', torus)
# swiss_roll = tadasets.swiss_roll(n=2000, r=4, ambient=10, noise=1.2)
# dsphere = tadasets.dsphere(n=1000, d=12, r=3.14, ambient=14, noise=0.14)
# infty_sign = tadasets.infty_sign(n=3000, noise=0.1)


def make_plot_from_fname(data):

    fig = plt.figure(figsize=plt.figaspect(1/2.))
    barcode_ax = fig.add_subplot(1, 2, 1)
    prot_ax = fig.add_subplot(1, 2, 2, projection='3d')

    prot_ax.scatter(*data.T)
    prot_ax.set_axis_off()


# data = np.random.random((100,2))
    ripped = ripser(data)
    lifespans = ripped['dgms']
# print(ripped['idx_perm'][10])
    print([x for i, x in enumerate(ripped['idx_perm']) if x != i])


    barcode_scatters = []

    def hover(event):
        if event.inaxes == barcode_ax:
            for betti_dim, barcode_scatter in enumerate(barcode_scatters):
                cont, ind = barcode_scatter.contains(event)
                if cont:
                    affected_simplicies = ind['ind']
                    for affected_id in affected_simplicies:
                        print(lifespans[betti_dim][affected_id])
                    # print(betti_dim, ind['ind'])

    fig.canvas.mpl_connect("motion_notify_event", hover)


    epsilon_max = 0
    for dim, dim_pts in enumerate(lifespans):
        # print(dim_pts.shape)
        sc = barcode_ax.scatter(*dim_pts.T, label=f"$H_{dim}$")
        barcode_scatters.append(sc)
        epsilon_max = max(epsilon_max, np.max(dim_pts[dim_pts != np.inf]))

    bounds = [-0.1 * epsilon_max, 1.1 * epsilon_max]
    barcode_ax.set_xlim(*bounds)
    barcode_ax.set_ylim(*bounds)
    barcode_ax.plot(bounds, bounds, "--")
    barcode_ax.legend()
    plt.show()




pdb_files = glob("./graph_theory_final/out/*.cif.npy")

for pdb_file in pdb_files:
    data = np.load(pdb_file)
    make_plot_from_fname(data)



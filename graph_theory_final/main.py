# mash of examples from github/scikit-tda/ripser.py and github/scikit-tda/persim
import numpy as np
from matplotlib import pyplot as plt
from ripser import ripser
from persim import plot_diagrams
from glob import glob
from tqdm import tqdm
from os.path import basename
from multiprocessing import get_context, Pool
from os import sys

# import tadasets
#
# torus = tadasets.torus(n=1200, c=2, a=1, noise=0.05)
# np.save('torus.npy', torus)
# swiss_roll = tadasets.swiss_roll(n=2000, r=4, ambient=10, noise=1.2)
# dsphere = tadasets.dsphere(n=1000, d=12, r=3.14, ambient=14, noise=0.14)
# infty_sign = tadasets.infty_sign(n=3000, noise=0.1)


def pickle_memoize(fname, creation_callback, verbose=False):
    """
    Try to read data from the pickle at `fname`, and save the output of
    `creation_callback` to `fname` as a pickle if `fname` doesn't exist.
    """
    from pickle import load, dump, UnpicklingError
    if verbose: print(f"looking for pickle file '{fname}'...")
    try:
        with open(fname, 'rb') as rf:
            if verbose: print(f"    found pickle file '{fname}'! :)) loading it...")
            return load(rf)
    except (FileNotFoundError, UnpicklingError):
        if verbose: print(f"    did not find pickle file '{fname}' or it was corrupted :( making it...")
    # except UnpicklingError:
    #     if verbose: print(f"    pickle file was corrupted! remaking...")
        got = creation_callback()
        try:
            with open(fname, 'wb') as wf:
                dump(got, wf)
            if verbose: print(f"    successfully made pickle file '{fname}'! :)")
        except TypeError as err:
            from sys import stderr
            if verbose: print("couldn't pickle the object! :(", err, file=stderr)
        return got


def make_lifespans(spatial_points):
# data = np.random.random((100,2))
    ripped = ripser(spatial_points, maxdim=1)
    print("made lifespans!!")
    lifespans = ripped['dgms']
    print("made lifespans!!")
# print(ripped['idx_perm'][10])
    # print([x for i, x in enumerate(ripped['idx_perm']) if x != i])

    return lifespans


def make_plot_from_fname(lifespans, data, pdb_file):

    fig = plt.figure(figsize=plt.figaspect(1/2.))
    barcode_ax = fig.add_subplot(1, 2, 1)
    prot_ax = fig.add_subplot(1, 2, 2, projection='3d')

    prot_ax.scatter(*data.T)
    prot_ax.set_axis_off()


    barcode_scatters = []

    # def hover(event):
    #     if event.inaxes == barcode_ax:
    #         for betti_dim, barcode_scatter in enumerate(barcode_scatters):
    #             cont, ind = barcode_scatter.contains(event)
    #             if cont:
    #                 affected_simplicies = ind['ind']
    #                 for affected_id in affected_simplicies:
    #                     print(lifespans[betti_dim][affected_id])
    #                 # print(betti_dim, ind['ind'])
    #
    # fig.canvas.mpl_connect("motion_notify_event", hover)


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
    barcode_ax.legend(loc="lower right")
    plt.title(pdb_file)
    # plt.savefig(f'imgs/junk/{pdb_file.replace("/", "-")}.png')
    plt.savefig(f'imgs/NEUROTRANSMITTER/NEW_{pdb_file.replace("/", "-")}.png')


    # plt.show()



def pool_worker(fname):
    print("entered pool worker!")
    spatial_points = np.load(fname)
    # spatial_points = np.random.rand(500, 3)
    print("loaded points")
    # lifespans = pickle_memoize(f"temp/{basename(fname)}.lifespans_pkl", lambda: make_lifespans(spatial_points))
    lifespans = make_lifespans(spatial_points)
    print("computed barcodes")
    make_plot_from_fname(lifespans, spatial_points, fname)
    print("plotted!")

if __name__ == '__main__':
    pdb_files = glob("./out/NEUROTRANSMITTER/*.cif.npy")
    # spatial_points_all = (np.load(pdb_file) for pdb_file in pdb_files)
    # with get_context("spawn").Pool(1) as p:
    for i, fname in enumerate(tqdm(pdb_files[:100])):
        # if i % 4 == int(sys.argv[1]) or True:
        if True:
            print(f"starting on {i} {fname}")
            pool_worker(fname)
    # list(tqdm(map(pool_worker, pdb_files), total=len(pdb_files)))
        # print("pool initialized!")
        # for lifespans, points, pdb_file in zip(
        #         tqdm(p.imap(make_lifespans, spatial_points_all), total=len(pdb_files)),
        #         spatial_points_all,
        #         pdb_files
        #         ):
        #     make_plot_from_fname(lifespans, points, pdb_file)
    # print("finished computing homology")



# with tqdm(total=len(pdb_files)) as pbar:
#     print("initializing...")
#     for pdb_file in tqdm(pdb_files[:8]):
#         pbar.set_description(pdb_file)
#         pbar.update(1)
#
#         spatial_points = np.load(pdb_file)
#
#         if spatial_points.size > 1000:
#             continue
#
#         print(spatial_points.size)
#
#         # lifespans = pickle_memoize(f"temp/{basename(pdb_file)}.lifespans_pkl", lambda: make_lifespans(spatial_points))
#         lifespans = make_lifespans(spatial_points)
#
#         make_plot_from_fname(lifespans, spatial_points, pdb_file)
#


# mash of examples from github/scikit-tda/ripser.py and github/scikit-tda/persim
import numpy as np
from matplotlib import pyplot as plt
from ripser import ripser
from persim import plot_diagrams

import tadasets

torus = tadasets.torus(n=200, c=2, a=1, noise=0.1)
# swiss_roll = tadasets.swiss_roll(n=2000, r=4, ambient=10, noise=1.2)
# dsphere = tadasets.dsphere(n=1000, d=12, r=3.14, ambient=14, noise=0.14)
# infty_sign = tadasets.infty_sign(n=3000, noise=0.1)

print(torus.shape)

# ax = plt.figure().add_subplot(projection='3d')
# ax.scatter(*torus.T)
# plt.show()

print("ripping...")

fig, barcode_ax = plt.subplots()
prot_ax = fig.add_subplot(projection='3d')

prot_ax.scatter(*torus.T)
prot_ax.set_axis_off()

# data = np.random.random((100,2))
ripped = ripser(torus)
lifespans = ripped['dgms']
epsilon_max = 0
for dim, dim_pts in enumerate(lifespans):
    print(dim_pts.shape)
    barcode_ax.scatter(*dim_pts.T, label=f"$H_{dim}")
    epsilon_max = max(epsilon_max, np.max(dim_pts))

barcode_ax.set_xlim(0, epsilon_max)
barcode_ax.set_ylim(0, epsilon_max)
barcode_ax.plot([0, epsilon_max], [0, epsilon_max], style="--")
barcode_ax.legend()
plt.show()
# print([x.shape for x in diagrams])
# plot_diagrams(diagrams, show=True)
plt.savefig("out.png")


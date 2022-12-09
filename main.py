# mash of examples from github/scikit-tda/ripser.py and github/scikit-tda/persim
import numpy as np
from matplotlib import pyplot as plt
from ripser import ripser
from persim import plot_diagrams

import tadasets

torus = tadasets.torus(n=4000, c=2, a=1, noise=0.1)
swiss_roll = tadasets.swiss_roll(n=2000, r=4, ambient=10, noise=1.2)
dsphere = tadasets.dsphere(n=1000, d=12, r=3.14, ambient=14, noise=0.14)
infty_sign = tadasets.infty_sign(n=3000, noise=0.1)

print(torus.shape)

# ax = plt.figure().add_subplot(projection='3d')
# ax.scatter(*torus.T)
# plt.show()
#
print("ripping...")

# data = np.random.random((100,2))
diagrams = ripser(torus)['dgms']
print([x.shape for x in diagrams])
plot_diagrams(diagrams, show=True)
plt.savefig("out.png")


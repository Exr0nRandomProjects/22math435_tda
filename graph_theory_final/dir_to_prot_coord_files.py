from protein_to_points import protein_to_points
import os
import numpy as np

files = os.listdir("./inp/ION_CHANNELS")

for f in files:
    d = protein_to_points(f"./inp/ION_CHANNELS/{f}")
    np.save(f"./out/ION_CHANNELS/{f}.npy", d)

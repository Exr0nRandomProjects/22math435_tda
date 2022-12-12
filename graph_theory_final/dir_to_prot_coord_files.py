from protein_to_points import protein_to_points
import os
import numpy as np

files = os.listdir("./inp")

for f in files:
    d = protein_to_points(f"./inp/{f}")
    np.save(f"./out/{f}.npy", d)

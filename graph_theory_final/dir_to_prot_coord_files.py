from protein_to_points import protein_to_points
import os
import numpy as np

files = os.listdir("./inp/IN_THE_DOC")

for f in files:
    d = protein_to_points(f"./inp/IN_THE_DOC/{f}")
    np.save(f"./out/IN_THE_DOC/{f}.npy", d)

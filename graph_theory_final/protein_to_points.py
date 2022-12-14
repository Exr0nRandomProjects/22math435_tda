import os

from Bio.PDB import MMCIFParser, PDBList, PDBParser
from Bio.PDB.DSSP import DSSP
from typing import *
import copy
import math
from parse import parse

import numpy as np
from Bio.PDB.ResidueDepth import get_surface
from tqdm import tqdm
import time
import matplotlib.pyplot as plt
import numpy as np

from main import *

def protein_to_points(filepath):
    def load(protein, filepath=False):
        """
        Arguments:
        - `protein`: A RCSB PDB code (or a filename if `filepath` is True.)
        - `filepath`: Specifies if the protein provided is a PDB code or file path.

        Returns: tuple containing the Bio.PDB.Structure and DSSP result of the protein. Not meant for direct consumption but to be fed into other tesslib functions.
        """

        path = None
        if not filepath:
            if len(protein) != 4:
                raise ValueError("Invalid RCSB protein ID (length != 4).")

            path = f"./pdb/{protein.lower()}.cif"
            if not os.path.isfile(path):
                pdbl = PDBList()
                pdbl.retrieve_pdb_file(protein, pdir="./pdb")
        else:
            if not os.path.isfile(protein):
                raise ValueError("Invalid path provided.")
            path = protein

        parser = MMCIFParser()
        if (protein.endswith(".pdb")):
            parser = PDBParser()

        structure = parser.get_structure(protein if not filepath else os.path.basename(protein), path)

        if not filepath:
            dssp = DSSP(structure[0], f"./pdb/{protein.lower()}.cif")
        else:
            dssp = DSSP(structure[0], protein)
        return (structure, dssp)


    def enc_point(p):
        return (p[0], p[1], p[2])

    class Point:
        """A single point, usually representing an atom.

        Should solely be used inside a `Protein`.
        """

        def __init__(self, coords: List[float]):
            """Initializes a point from coordinates.

            Arguments:
            - `coords`: a list of length 3 containing coordinates.
            """
            self.coords = coords
            self.data = { }

        def __getitem__(self, key: str) -> Any:
            """Returns the value associated with a given attribute.

            Arguments:
            - `key`: the name of the desired attribute.

            Returns: the value of the attribute for the Point.
            """
            return self.data[key]

        def __setitem__(self, key: str, val: Any):
            """Sets the value associated with a given attribute.

            Arguments:
            - `key`: the name of the desired attribute.
            - `val`: the value of the attribute for the Point.).reshape(4,3)
            """
            self.data[key] = val

        def loc(self):
            return (self.coords[0], self.coords[1], self.coords[2])


    class Protein:
        """Class representing a protein and its associated attributes/tessellations."""

        def __init__(self, structure, dssp, multimer = False):
            """Initializes a protein from its atoms.
            Arguments:
            - `structure`: the Structure object returned by BioPython's PDB parser.
            - `dssp`: the DSSP object returned by BioPython's PDB.DSSP module
            - `multimer`: whether to use all chains or just the first. Defaults to False.
            """
            self.points = []
            orig = None
            ca_count = 0
            model = next(structure.get_models())
            chains = [next(model.get_chains())] if not multimer else model.get_chains()
            for chain_num, chain in enumerate(chains):
                for atom in chain.get_atoms():
                    self.points.append(Point(atom.get_coord()))
                    self.points[-1]["type"] = atom.get_name()
                    self.points[-1]["amino"] = atom.get_parent()
                    code = None
                    if atom.get_name() == "CA" and atom.get_parent().get_resname() != "CA":
                        idx = (ca_count := ca_count + 1)-1
                        code = dssp[list(dssp.keys())[idx]][2]
                        self.points[-1]["surface"] = dssp[list(dssp.keys())[idx]][3]
                    else:
                        self.points[-1]["surface"] = None
                    for i in ("-", "G", "H", "I", "T", "E", "B", "S"):
                        self.points[-1][f"sec_{i}"] = 1 if code == i else 0
                    self.points[-1]["cells"] = []
                    self.points[-1]["tet_classes"] = []
                    if multimer:
                        self.points[-1]["chain"] = chain_num
                self.tessellations = []
                self.filtered = []

        def __iter__(self):
            self.__n = 0
            return self

        def __next__(self) -> Point:
            """Gets the next `Point` composing the `Protein`.

            Returns: the next `Point`.
            """
            if self.__n < len(self.points):
                self.__n += 1
                return self.points[self.__n - 1]
            else:
                raise StopIteration

        def attributes(self) -> List[str]:
            """Returns a list of all known attributes."""
            return list(self.points[0].data.keys())

        def add_attribute(self, name: str, rule: Callable[[Point], Any]):
            """Adds an attribute to each point.

            Arguments:
            - `name`: the name of the attribute
            - `rule`: a function that returns the attribute of a given `Point`. Can theoretically return anything; really should return a string or a float.
            """
            for point in self.points:
                point[name] = rule(point)


    def isresidue(point: Point) -> bool:
        """Checks whether or not the point is associated with a residue.

        Helpful for filtering out ions that still have carbon alphas!

        Arguments:
        - `point`: the Point in question.

        Returns: True if residue, False if not.
        """
        return point["amino"].get_resname() in (
            "ILE",
            "VAL",
            "LEU",
            "PHE",
            "CYS",
            "MET",
            "ALA",
            "GLY",
            "THR",
            "SER",
            "TRP",
            "TYR",
            "PRO",
            "HIS",
            "GLU",
            "GLN",
            "ASP",
            "ASN",
            "LYS",
            "ARG",
        )

    s, d = load(filepath, filepath=True)
    p = Protein(s, d)


    def rule(x):
        # print(vars(x))
        # print(x)
        # if x["type"] != "CA" and x["type"] != "CB":
        #     return False

        # if not isresidue(x):
        #     return False

        # return True
        # return isresidue(x)
        # return x["surface"] != None
        cond = x["type"] == "CA" and isresidue(x)

        # if cond:
        #     print(x["surface"])

        return cond # and x["surface"] >= 0.5

    # rule = lambda x: x["type"] == "CA" and isresidue(x)
    # rule = lambda x: x["type"] == "CA" or isresidue(x)
    # rule = lambda x: isresidue(x)

    filtered = list(filter(rule, p.points))

    d = np.asarray([point.coords for point in filtered])
    print(len(d))

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')

    # ax.scatter(d[:,0], d[:,1], d[:,2])

    # plt.show()

    return d

if __name__ == '__main__':
# Q4KMQ2
    # neg_points = protein_to_points("./inp/IN_THE_DOC/AF-Q9Y4I1.cif")
    # pos_points = protein_to_points("./inp/NEUROTRANSMITTER/AF-P23978.cif")
    pos_points = protein_to_points("./inp/NEUROTRANSMITTER/AF-Q08469.cif")

    lifespans = make_lifespans(pos_points)
    make_plot_from_fname(lifespans, pos_points, "yuh.")

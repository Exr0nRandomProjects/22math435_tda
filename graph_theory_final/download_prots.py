import wget


cat_name = "IN_THE_DOC"
# ids = "P05177,P10635,P10632,Q96C23,Q9HCG7,P04062,P47712,P08684,P05108,O15269,O95237,O95470,P07902,P22570,P49619,Q5VZY2,Q6ZNC8,Q9HAY6,Q9HBH5,Q9NUN7,Q02928,P11712,Q6GTS8,O14495,P33260,P51589,P52429,Q16760,Q53H12,Q5KSL6,Q8NEB5,Q8NFR3,Q8NFU5,Q8TDN7,Q92781,Q969W0,Q9HA82,Q9NUV7,Q9Y6T7,Q9NZ01,Q16678,A5PLL7,Q02083,Q8N9I5,Q96SQ9,Q9UJ83,O60218,O75912,P52824,Q13574,Q14376,Q16739,Q643R3,Q6P531,Q7Z449,Q8N5B7,Q8NBN7,Q9HCS2,Q9UHH9,Q9UJ14,O14494,O75907,P00352,P47895,P50053,P51687,Q06136,Q6P1A2,Q6ZMG9,Q6ZWT7,Q8NF37,Q96G23,Q96N66,Q96NR8,Q9NY59,Q9UHE5,Q9UHK6,P04798,A6NGU5,O15270,O43688,O60906,O75452,O94788,P23743,P27544,P49675,P51570,Q13510,Q5QJU3,Q7L5N7,Q86XP1,Q8IU89,Q8IZV5,Q8N3Y7,Q8TC12,Q99999,Q9NR71,O00154,P98187".split(",")
ids = "Q9NY46,Q9Y5Y9,Q9UI33,Q15858,O00180,Q7Z418,Q9Z2T2,Q91WD2,Q4KMQ2,Q7Z3S7,P61981,Q6UXB4,P50148,Q2PKF4,Q12965,O00160,Q8WXR4,Q9Y2K3,Q9Y623,Q9Y4I1,Q13459,Q9HD67".split(",")


for ID in ids:
    try:
        wget.download(
            f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-model_v2.cif",
            f"./inp/{cat_name}/AF-{ID}.cif")
    except:
        print("no good!", ID)







































exit(0)

from parse import parse
import rich
import re
import io
import rich.console
from rich import print
import click
import pandas as pd
from bioservices import UniProt
from pypdb.clients.pdb.pdb_client import get_pdb_file, PDBFileType
import os
import wget
import requests
import gzip
import tqdm
import traceback
from requests.adapters import HTTPAdapter, Retry
import shutil

def main():
    # if not os.path.exists("stream.tar.gz"):
    #     url = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Clength%2Cxref_prints%2Cft_transmem%2Cxref_alphafolddb%2Cxref_pdb&format=tsv&query=%28go%3A0008305%29"
    #     wget.download(url, "stream.tar.gz")
    # df = pd.read_table(gzip.open("stream.tar.gz"))
    # pdb_out_dir = "/mnt/db/data/input/interns/int_pdbs/pos"
    # fasta_out_dir = "/mnt/db/data/input/interns/int_fastas/pos"
    # fetched = 0
    # found = 0
    # # df.to_csv("./whee.csv")
    # # exit(1)
    # for i, row in tqdm.tqdm(df.iterrows()):
    #     if "disintegrin" in row["Protein names"].lower(): continue
    #     # if type(row["PDB"]) == str:
    #     #     found+=1
    #     #     for ident in row["PDB"][:-1].split(';'):
    #     #         if not os.path.exists(f"{pdb_out_dir}/{ident}.cif"):
    #     #             pdb_str = get_pdb_file(ident.lower(), filetype=PDBFileType.CIF, compression=True)
    #     #             pdb_file = open(f"{pdb_out_dir}/{ident}.cif", "wb")
    #     #             pdb_file.write(pdb_str)
    #     #             pdb_file.close()
    #     #             fetched += 1
    #     # # elif type(row["AlphaFoldDB"]) == str:
    #     # #     found+=1
    #     # #     for ident in row["AlphaFoldDB"][:-1].split(';'):
    #     # #         if not os.path.exists(f"{pdb_out_dir}/AF_{ident}.cif"):
    #     # #             wget.download(f"https://alphafold.ebi.ac.uk/files/AF-{ident}-F1-model_v2.cif", f"{pdb_out_dir}/AF_{ident}.cif")
    #     # #             fetched += 1
    #     # else:
    #     try:
    #         if not os.path.exists(f"{pdb_out_dir}/AF_{row['Entry']}.cif"):
    #             wget.download(f"https://alphafold.ebi.ac.uk/files/AF-{row['Entry']}-F1-model_v3.cif", f"{pdb_out_dir}/AF_{row['Entry']}.cif")
    #             print(f"Downloading {row['Entry']}")
    #             # wget.download(f"https://swissmodel.expasy.org/repository/uniprot/{row['Entry']}.pdb", f"{pdb_out_dir}/SM_{ident}.pdb")
    #     except:
    #         traceback.print_exc()
    #         pass


    # print("Exhausted results!")
    # print(found)
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    ranges = [(0, 500, 200), # 200
              (700, 840, 5000),
              (840, 1020, 800),
              (1020, 1100, 2500),
              (1100, 1230, 2500)]

    for low,up,num_integrins in ranges:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=NOT%20(go:0008305)%20AND%20(length:[{low}%20TO%20{up}])&format=tsv&fields=accession%2Cid%2Cprotein_name%2Clength%2Cft_transmem%2Cxref_alphafolddb%2Cxref_pdb"

        pdb_out_dir = "/mnt/db/data/input/interns/int_pdbs/raw/nneg"
        fasta_out_dir = "/mnt/db/data/input/interns/int_fastas/neg"
        fetched = 0
        print("ok srsly")

        for batch, total in get_batch(url):
            print('batch??')
            df2 = pd.read_table(io.StringIO(batch.text))

            for i, row in tqdm.tqdm(df2.iterrows()):
                if "disintegrin" in row["Protein names"].lower(): continue
                if type(row["PDB"]) == str:
                    for ident in row["PDB"][:-1].split(';'):
                        if not os.path.exists(f"{pdb_out_dir}/{ident}.cif"):
                            pdb_str = get_pdb_file(ident.lower(), filetype=PDBFileType.CIF, compression=True)
                            pdb_file = open(f"{pdb_out_dir}/{ident}.cif", "wb")
                            pdb_file.write(pdb_str)
                            pdb_file.close()
                            fetched += 1
                    # else:
                    #     print("found")
                    if not os.path.exists(f"{pdb_out_dir}/AF_{row['Entry']}.cif"):
                        try:
                            wget.download(f"https://alphafold.ebi.ac.uk/files/AF-{row['Entry']}-F1-model_v2.cif", f"{pdb_out_dir}/AF_{ident}.cif")
                            fetched += 1
                        except:
                            pass

            if fetched > num_integrins:
                break

        # elif type(row["Transmembrane"]) == str:
        #     fasta_req = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=5&sort=score&accession={row['Entry']}"
        #     r = requests.get(fasta_req, headers={"Accept": "text/x-fasta"})
        #     region = [int(i) for i in parse("TRANSMEM {}..{} {}", row["Transmembrane"])[:2]]
        #     raw_seq = r.text
        #     if region[0] > len(raw_seq) - region[1]:
        #         raw_seq = raw_seq[:region[0]-1]
        #     else:
        #         raw_seq = raw_seq[region[1]:]
        #     fasta = open(f"{fasta_out_dir}/{row['Entry']}.fasta", "w")
        #     fasta.write(raw_seq)
        #     fasta.close()


import pickle
if __name__ == "__main__":
    # pdb_orig_out_dir = "/mnt/db/data/input/interns/int_pdbs/pos"
    # pdb_out_dir = "/mnt/db/data/input/interns/int_pdbs/pos2"
    # l = pickle.load(open("/mnt/db/data/out.pkl", "rb"))
    # for acc in l:
    #     if not os.path.exists(f"{pdb_orig_out_dir}/AF_{acc}.cif"):
    #         try:
    #             wget.download(f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v3.cif", f"{pdb_out_dir}/AF_{acc}.cif")
    #             print(f"Downloading {acc}")
    #         except:
    #             pass
    #     else:
    #         shutil.copy(f"{pdb_orig_out_dir}/AF_{acc}.cif", pdb_out_dir)
    main()




import requests
import re
import os
import Compound

def is_valid(compound: Compound):
    central_atom = compound.find_central_atom()
    if not central_atom:
        return False
    # Check the compound has more than one unique element (central + ligands)
    if compound.numUniqueElements() < 2:
        return False

    return compound.isCovalent()

def get_cid(query):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/cids/TXT"
    r = requests.get(url)
    if r.ok and r.text.strip().isdigit():
        return r.text.strip()
    return None

def download_sdf(cid, name, out_dir="sdf_single_center"):
    os.makedirs(out_dir, exist_ok=True)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    r = requests.get(url)
    if r.ok:
        filename = os.path.join(out_dir, f"{name}.sdf")
        with open(filename, "wb") as f:
            f.write(r.content)
        print(f"[✓] Downloaded: {filename}")
    else:
        print(f"[!] Failed to download SDF for CID {cid}")

def process_file(filename):
    with open(filename, "r") as f:
        formulas = [line.strip() for line in f if line.strip()]
    for formula in formulas:
        cmpd = Compound(formula)
        if is_valid(cmpd):
            cid = get_cid(formula)
            if cid:
                safe_name = re.sub(r'[^\w\-]', '_', formula)
                download_sdf(cid, safe_name)
            else:
                print(f"[!] CID not found for: {formula}")
        else:
            print(f"[x] Skipped (not single-center covalent): {formula}")

# Example call:
# process_file("your_compounds.txt")
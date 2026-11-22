"""Prepare EGFR wild-type and mutant structures from the same 1M17 scaffold."""
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile

SRC = "1M17.pdb"
OUT = "prepared"
os.makedirs(OUT, exist_ok=True)

# clinical name -> (wild-type residue, residue number in 1M17, mutant residue)
# 1M17 numbering = clinical numbering - 24
MUTS = {
    "L858R": ("LEU", 834, "ARG"),
    "T790M": ("THR", 766, "MET"),
    "C797S": ("CYS", 773, "SER"),
    "G719S": ("GLY", 695, "SER"),
}
VARIANTS = {
    "WT": [],
    "L858R": ["L858R"],
    "T790M": ["T790M"],
    "C797S": ["C797S"],
    "G719S": ["G719S"],
    "L858R_T790M": ["L858R", "T790M"],
}

def prepare(name, mut_names):
    fixer = PDBFixer(filename=SRC)
    fixer.removeHeterogens(keepWater=False)      # removes erlotinib and waters
    if mut_names:
        specs = [f"{MUTS[m][0]}-{MUTS[m][1]}-{MUTS[m][2]}" for m in mut_names]
        fixer.applyMutations(specs, "A")
    fixer.findMissingResidues()
    fixer.missingResidues = {}                    # do not invent missing loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    path = os.path.join(OUT, f"EGFR_{name}.pdb")
    with open(path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    print("wrote", path)

for name, muts in VARIANTS.items():
    prepare(name, muts)

# co-crystal erlotinib, used later as the reference pose for re-docking
with open(SRC) as src, open(os.path.join(OUT, "erlotinib_crystal.pdb"), "w") as dst:
    for line in src:
        if line.startswith("HETATM") and line[17:20].strip() == "AQ4":
            dst.write(line)
    dst.write("END\n")
print("wrote erlotinib_crystal.pdb")

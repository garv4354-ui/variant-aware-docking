import os, glob, time, subprocess, urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem

DRUGS = ["erlotinib", "gefitinib", "afatinib", "dacomitinib", "osimertinib",
         "lapatinib", "brigatinib", "aspirin", "caffeine", "ibuprofen"]

os.makedirs("ligands", exist_ok=True)
os.makedirs("receptors", exist_ok=True)

for name in DRUGS:
    sdf2d = f"ligands/{name}_2d.sdf"
    if not os.path.exists(sdf2d):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF"
        urllib.request.urlretrieve(url, sdf2d)
        time.sleep(0.5)
    mol = Chem.SDMolSupplier(sdf2d, removeHs=False)[0]
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        print("FAILED to embed", name); continue
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    sdf3d = f"ligands/{name}.sdf"
    w = Chem.SDWriter(sdf3d); w.write(mol); w.close()
    subprocess.run(["obabel", sdf3d, "-O", f"ligands/{name}.pdbqt"],
                   check=True, capture_output=True)
    print("ligand ready:", name)

for pdb in sorted(glob.glob("../03_structures/prepared/EGFR_*.pdb")):
    tag = os.path.basename(pdb).replace(".pdb", "")
    subprocess.run(["obabel", "-ipdb", pdb, "-xr", "-O", f"receptors/{tag}.pdbqt"],
                   check=True, capture_output=True)
    print("receptor ready:", tag)

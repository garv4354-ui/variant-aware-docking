load 03_structures/prepared/EGFR_WT.pdb, wt
load 03_structures/AF_aligned_trimmed.pdb, af
hide everything
bg_color white
show cartoon, wt or af
color grey70, wt
spectrum b, red_yellow_green_cyan_blue, af, minimum=50, maximum=90
set ray_shadow, off
png results/figures/alphafold_vs_crystal.png, width=1600, height=1200, dpi=300, ray=1
quit

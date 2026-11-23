load 03_structures/prepared/EGFR_WT.pdb, wt
load 03_structures/prepared/EGFR_T790M.pdb, mut
load 03_structures/prepared/erlotinib_crystal.pdb, xtal
align mut, wt
hide everything
bg_color white
show cartoon, wt or mut
color grey80, wt
color palecyan, mut
show sticks, resi 766
color magenta, wt and resi 766 and elem C
color orange, mut and resi 766 and elem C
show sticks, xtal
color yellow, xtal and elem C
hide (hydro)
zoom resi 766, 10
set ray_shadow, off
png results/figures/t790m_site.png, width=1600, height=1200, dpi=300, ray=1
quit

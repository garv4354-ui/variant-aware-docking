load 03_structures/prepared/EGFR_WT.pdb, wt
load results/figures/wt_erlotinib_pose.pdb, pose
load 03_structures/prepared/erlotinib_crystal.pdb, xtal
hide everything
bg_color white
show cartoon, wt
color grey80, wt
select pocket, wt within 5 of xtal
show sticks, pocket
color grey50, pocket and elem C
show sticks, wt and resi 766
color magenta, wt and resi 766 and elem C
show sticks, pose
show sticks, xtal
color cyan, pose and elem C
color yellow, xtal and elem C
hide (hydro)
zoom xtal, 8
set ray_shadow, off
png results/figures/erlotinib_pose.png, width=1600, height=1200, dpi=300, ray=1
quit

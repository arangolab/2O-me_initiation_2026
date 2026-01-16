# Load your structure (adjust if needed)
fetch 8pj3, async=0
bg_color white
hide everything

# Selections
select mRNA, chain 7 and resi -3+-2+-1
select 18S, chain A and resi 961
select eIF2a, chain r and resi 54
select S5, chain V and resi 130

# Set ribbon width to make the ribbons thicker
set ribbon_width, 20

# Show ribbons for the selected chains
show ribbon, mRNA
show ribbon, 18S
show ribbon, eIF2a
show ribbon, S5

# Show sticks for the selected chains
show sticks, mRNA
show sticks, 18S
show sticks, eIF2a
show sticks, S5

# Colors for ribbons (make sure to color ribbons based on chain, not atoms)
color bluewhite, mRNA
color yelloworange, 18S
color black, eIF2a
color wheat, S5

# Atom-specific colors (not affecting ribbons)
color red, elem O
color blue, elem N
color orange, elem P

# Color ribbons based on the chain selections
set ribbon_color, bluewhite, mRNA  
set ribbon_color, yelloworange, 18S  
set ribbon_color, black, eIF2a  
set ribbon_color, wheat, S5 

# Visual settings
set specular, 0.0         
set shininess, 10         
set ambient, 0.3          
set reflect, 0.0          
set direct, 0.5          
set ray_trace_mode, 1

# Hydrogen bonds (between selections)
dist hbonds_mRNA_tRNA, (donor and mRNA), (acceptor and tRNA), mode=2, cutoff=3.2
dist hbonds_tRNA_mRNA, (donor and tRNA), (acceptor and mRNA), mode=2, cutoff=3.2
dist hbonds_mRNA_eIF1A, (donor and eIF2a), (acceptor and mRNA), mode=2, cutoff=3.2
dist hbonds_eIF2a_mRNA, (donor and mRNA), (acceptor and eIF2a), mode=2, cutoff=3.2

# Dash styling (optional)
set dash_color, yellow
set dash_width, 2.0
set dash_gap, 0.25

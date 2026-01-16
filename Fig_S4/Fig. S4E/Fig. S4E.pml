# Load your structure (adjust if needed)
fetch 8pj1, async=0
bg_color white
hide everything

# Selections
select mRNA, chain 7 and resi 1+2+3
select tRNA, chain w and resi 34+35+36
select 18S, chain A and resi 1830

# Set ribbon width to make the ribbons thicker
set ribbon_width, 20

# Show ribbons for the selected chains
show ribbon, mRNA
show ribbon, tRNA
show ribbon, 18S

# Show sticks for the selected chains
show sticks, mRNA
show sticks, tRNA
show sticks, 18S

# Colors for ribbons (make sure to color ribbons based on chain, not atoms)
color bluewhite, mRNA
color green, tRNA
color yelloworange, 18S

# Atom-specific colors (not affecting ribbons)
color red, elem O
color blue, elem N
color orange, elem P

# Color ribbons based on the chain selections
set ribbon_color, bluewhite, mRNA  
set ribbon_color, green, tRNA  
set ribbon_color, yelloworange, 18S  

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
dist hbonds_mRNA_18S, (donor and 18S), (acceptor and mRNA), mode=2, cutoff=3.2

# Dash styling (optional)
set dash_color, yellow
set dash_width, 2.0
set dash_gap, 0.25

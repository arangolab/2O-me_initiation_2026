# Load your structure (adjust if needed)fetch 8pj3, async=0bg_color whitehide everything# Selectionsselect mRNA, chain 7 and resi 1+2+3select tRNA, chain w and resi 34+35+36select eIF1A, chain q and resi 6+7+8
# Set ribbon width to make the ribbons thicker
set ribbon_width, 20
# Show ribbons for the selected chains
show ribbon, mRNAshow ribbon, tRNAshow ribbon, eIF1A# Show sticks for the selected chainsshow sticks, mRNAshow sticks, tRNAshow sticks, eIF1A# Colors for ribbons (make sure to color ribbons based on chain, not atoms)
color bluewhite, mRNA
color green, tRNA
color black, eIF1A

# Atom-specific colors (not affecting ribbons)
color red, elem O
color blue, elem N
color orange, elem P

# Color ribbons based on the chain selections
set ribbon_color, bluewhite, mRNA  
set ribbon_color, green, tRNA  
set ribbon_color, black, eIF1A  # Visual settingsset specular, 0.0         
set shininess, 10         
set ambient, 0.3          
set reflect, 0.0          
set direct, 0.5          set ray_trace_mode, 1# Hydrogen bonds (between selections)#dist hbonds_mRNA_tRNA, (donor and mRNA), (acceptor and tRNA), mode=2, cutoff=3.2
#dist hbonds_tRNA_mRNA, (donor and tRNA), (acceptor and mRNA), mode=2, cutoff=3.2#dist hbonds_mRNA_eIF1A, (donor and eIF1A), (acceptor and mRNA), mode=2, cutoff=3.2# Dash styling (optional)#set dash_color, yellow#set dash_width, 2.0#set dash_gap, 0.25
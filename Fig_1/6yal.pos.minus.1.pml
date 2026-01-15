# Load your structure (adjust if needed)fetch 6yal, async=0bg_color whitehide everything# Selectionsselect mRNA, chain 3 and resi 51+52+53+54select tRNA, chain 1 and resi 34+35+36+37# Set ribbon width to make the ribbons thickerset ribbon_width, 20# Show ribbons for the selected chainsshow ribbon, mRNAshow ribbon, tRNA# Show sticks for the selected chainsshow sticks, mRNAshow sticks, tRNA# Colors for ribbons (make sure to color ribbons based on chain, not atoms)color bluewhite, mRNAcolor green, tRNA# Atom-specific colors (not affecting ribbons)color red, elem Ocolor blue, elem Ncolor orange, elem P# Color ribbons based on the chain selectionsset ribbon_color, bluewhite, mRNA  set ribbon_color, green, tRNA  # Visual settingsset specular, 0.0         set shininess, 10         set ambient, 0.3          set reflect, 0.0          set direct, 0.5          set ray_trace_mode, 1# Hydrogen bonds (between selections)#dist hbonds_mRNA_tRNA, (donor and mRNA), (acceptor and tRNA), mode=2, cutoff=3.2#dist hbonds_tRNA_mRNA, (donor and tRNA), (acceptor and mRNA), mode=2, cutoff=3.2

# Measure Angle for Acceptor O14 (chain 1, resi 37)
#select acceptor_atom, chain 1 and resi 37 and name O14
#select donor_atom, chain 3 and resi 51 and name O2'
#select base_neighbor, chain 1 and resi 37 and name C14
#distance hbond_dist, donor_atom, acceptor_atom
#angle hbond_base, donor_atom, acceptor_atom, base_neighbor

# Measure Angle for Acceptor ODB (chain 1, resi 37)
#select acceptor_atom, chain 1 and resi 37 and name ODB
#select donor_atom, chain 3 and resi 51 and name O2'
#select base_neighbor, chain 1 and resi 37 and name C13
#distance hbond_dist, donor_atom, acceptor_atom
#angle hbond_base, donor_atom, acceptor_atom, base_neighbor
# Dash styling (optional)#set dash_color, yellow#set dash_width, 2.0#set dash_gap, 0.25
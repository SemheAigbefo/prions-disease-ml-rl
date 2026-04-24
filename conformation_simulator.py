"""
Conformation Simulator for Prion Disease Research
This module based on the HP Lattice Model simulatea folding
Proteins are chains of H and P monomers placed on a 2D grid. 
The energy function counts H-H contacts between non-bonded neighbors
 — each contact contributes -1 to energy, so lower energy means more stable. 
 You implement move operations — pivot moves, end moves, 
 crankshaft moves — that define how the chain transitions between conformations one step at a time. 
 This simulator is your gym environment, the thing everything else runs on.

"""
#When two H residues find each other, the system releases energy — it becomes more stable. Lower energy = more stable = more folded.
def is_valid (chain):
    chain_set = set(chain)
    return len(chain_set) == len(chain)


def energy_func(sequence, chain):
    energy = 0
    
    for i in range(len(chain)):
        for j in range(i+1, len(chain)):
            if sequence[i] == 'H' and sequence[j] == 'H':
                grid_dist = abs(chain[i][0] - chain[j][0]) + abs(chain[i][1] - chain[j][1])
                if grid_dist == 1:
                    if j > i + 1:
                        energy -= 1
    
    return energy
                
# # Case 1: straight chain, no contacts possible
# sequence = "HPH"
# chain = [(0,0), (1,0), (2,0)]
# print(energy_func(sequence, chain))  # should print 0

# # Case 2: bent chain, two H beads touch
# sequence = "HPH"
# chain = [(0,0), (0,1), (1,0)]
# print(energy_func(sequence, chain))  # should print -1

# sequence = "HPHH"
# chain = [(0,0), (0,1), (1,1), (1,0)]
# print(energy_func(sequence, chain))

# Valid - no overlaps
print(is_valid([(0,0), (0,1), (1,1)]))  # should print True

# Invalid - bead 0 and bead 2 share position
print(is_valid([(0,0), (0,1), (0,0)]))  # should print False
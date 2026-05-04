"""
Conformation Simulator for Prion Disease Research
This module based on the HP Lattice Model simulatea folding
Proteins are chains of H and P monomers placed on a 2D grid. 
The energy function counts H-H contacts between non-bonded neighbors
 — each contact contributes -1 to energy, so lower energy means more stable. 
 You implement move operations — pivot moves, end moves, 
 crankshaft moves — that define how the chain transitions between conformations one step at a time. 
 This simulator is the thing everything else runs on.

"""
#When two H residues find each other, the system releases energy, it becomes more stable. lower energy = more stable = more folded.
def is_valid (chain): 
    """
    a chain is valid if it has no duplicates
    """
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

import random

def pivot_move(chain):
    n = len(chain)
    
    # need at least 3 beads for a pivot to make sense
    if n < 3:
        return None
    
    # pick a random bead that isn't the first or last
    pivot_idx = random.randint(1, n - 2)
    px, py = chain[pivot_idx] #index x and y of the random bead
    
    # three possible rotations around the pivot point
    rotations = [
        lambda x, y: (px + (y - py), py - (x - px)),   # 90 clockwise
        lambda x, y: (px - (y - py), py + (x - px)),   # 90 counterclockwise
        lambda x, y: (2*px - x, 2*py - y)              # 180 flip
    ]
    
    # pick a random rotation
    rotate = random.choice(rotations)
    
    # randomly pick which side of the pivot gets rotated
    # either everything before the pivot, or everything after
    new_chain = []
    for i in range(n):
        if i < pivot_idx: 
            # rotate the left side, keep right side
            new_chain.append(rotate(chain[i][0], chain[i][1]))
        else:
            new_chain.append(chain[i])
    
    # only return the new chain if it doesn't self intersect
    if is_valid(new_chain):
        return new_chain
    return None


def end_move(chain):
    n = len(chain)
    
    # randomly decide whether to move the first or last bead
    move_first = random.random() < 0.5
    
    if move_first:
        # the bead that moves and its reference neighbor
        moving_bead = chain[0]
        reference = chain[1]
    else:
        moving_bead = chain[-1]
        reference = chain[-2]
    
    # try all four directions adjacent to the reference bead
    directions = [(1,0), (-1,0), (0,1), (0,-1)]
    
    valid_moves = []
    for dx, dy in directions:
        candidate = (reference[0] + dx, reference[1] + dy)
        
        # candidate can't be where the moving bead already is
        if candidate == moving_bead:
            continue
            
        # build the new chain with candidate replacing moving bead
        if move_first:
            new_chain = [candidate] + list(chain[1:])
        else:
            new_chain = list(chain[:-1]) + [candidate]
            
        if is_valid(new_chain):
            valid_moves.append(new_chain)
    
    # return a random valid move or None if none exist
    if valid_moves:
        return random.choice(valid_moves)
    return None

import math
import random

def monte_carlo(sequence, chain, steps=10000, T_start=2.0, T_end=0.01):
    
    n = len(sequence)
    current_chain = [pos for pos in chain]
    current_energy = energy_func(sequence, current_chain)
    
    # best conformation found so far
    best_chain = current_chain
    best_energy = current_energy
    
    # this is your scientific data - energy at every single step
    trajectory = [current_energy]
    
    for step in range(steps):
        
        # calculate current temperature using cooling schedule
        T = T_start * (T_end / T_start) ** (step / steps)
        
        # randomly propose either a pivot or end move
        if random.random() < 0.7:
            new_chain = pivot_move(current_chain)
        else:
            new_chain = end_move(current_chain)
            
        # if move returned None, skip this step
        if new_chain is None:
            trajectory.append(current_energy)
            continue
        
        # calculate energy of proposed conformation
        new_energy = energy_func(sequence, new_chain)
        delta_e = new_energy - current_energy
        
        # acceptance decision
        if delta_e <= 0:
            # always accept better or equal moves
            current_chain = new_chain
            current_energy = new_energy
        else:
            # accept worse move with probability e^(-delta_e/T)
            if random.random() < math.exp(-delta_e / T):
                current_chain = new_chain
                current_energy = new_energy
        
        # record energy regardless of accept or reject
        trajectory.append(current_energy)
        
        # track the best conformation found across entire run
        if current_energy < best_energy:
            best_energy = current_energy
            best_chain = [pos for pos in current_chain]
    
    return {
        'best_chain': best_chain,
        'best_energy': best_energy,
        'final_chain': current_chain,
        'final_energy': current_energy,
        'trajectory': trajectory
    }
import matplotlib.pyplot as plt

# starting conformation - straight chain
sequence = "HPHPPHHPHPPHPHHPPHPH"
chain = [(i, 0) for i in range(len(sequence))]

# run one trial
result = monte_carlo(sequence, chain, steps=10000)
for trial in range(5):
    result = monte_carlo(sequence, chain, steps=10000)
    print(f"Trial {trial+1}: best energy = {result['best_energy']}")

# plot the energy trajectory
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def plot_trajectories(sequence, chain, n_trials=20, steps=10000):
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    
    all_results = []
    for _ in range(n_trials):
        result = monte_carlo(sequence, chain, steps=steps)
        all_results.append(result)
    
    final_energies = [r['final_energy'] for r in all_results]
    best_overall = min(r['best_energy'] for r in all_results)
    
    # --- LEFT PLOT: trajectories with color encoding ---
    ax1 = axes[0]
    
    # color each line by its final energy
    # darker = deeper trap, green = reached global minimum
    for r in all_results:
        traj = r['trajectory']
        x = range(0, len(traj), 20)  # sample every 20 steps
        y = traj[::20]
        
        if r['final_energy'] == best_overall:
            ax1.plot(x, y, color='#1D9E75', linewidth=1.5, 
                    alpha=0.9, zorder=3)
        else:
            # purple with alpha so overlapping lines show density
            ax1.plot(x, y, color='#7F77DD', linewidth=0.8, 
                    alpha=0.3, zorder=2)
    
    # add temperature curve on second y axis
    ax1_twin = ax1.twinx()
    steps_range = range(steps)
    temps = [2.0 * (0.01/2.0)**(s/steps) for s in steps_range]
    ax1_twin.plot(steps_range, temps, color='#EF9F27', 
                 linewidth=1, alpha=0.5, linestyle='--')
    ax1_twin.set_ylabel('temperature', color='#EF9F27', fontsize=10)
    ax1_twin.tick_params(colors='#EF9F27')
    ax1_twin.set_ylim(0, 2.5)
    
    ax1.set_xlabel('step')
    ax1.set_ylabel('energy')
    ax1.set_title(f'trajectories — {n_trials} trials\n'
                  f'green = reached global min, purple = trapped')
    ax1.axhline(y=best_overall, color='#1D9E75', 
               linestyle=':', linewidth=1, alpha=0.7)
    ax1.text(steps*0.98, best_overall+0.1, 
            f'global min: {best_overall}', 
            ha='right', fontsize=9, color='#1D9E75')
    
    # --- RIGHT PLOT: histogram of final energies ---
    ax2 = axes[1]
    
    energy_counts = {}
    for e in final_energies:
        energy_counts[e] = energy_counts.get(e, 0) + 1
    
    energies = sorted(energy_counts.keys())
    counts = [energy_counts[e] for e in energies]
    colors = ['#1D9E75' if e == best_overall else '#7F77DD' 
             for e in energies]
    
    bars = ax2.bar(energies, counts, color=colors, 
                  edgecolor='white', linewidth=0.5, width=0.6)
    
    # label each bar with count
    for bar, count in zip(bars, counts):
        ax2.text(bar.get_x() + bar.get_width()/2, 
                bar.get_height() + 0.2,
                str(count), ha='center', va='bottom', 
                fontsize=10, fontweight='bold')
    
    ax2.set_xlabel('final energy')
    ax2.set_ylabel('number of runs')
    ax2.set_title(f'where chains get trapped\n'
                 f'green = global minimum, purple = metastable trap')
    ax2.set_xticks(energies)
    
    # annotate the most common trap
    most_common = max(energy_counts, key=energy_counts.get)
    if most_common != best_overall:
        ax2.annotate('most common trap', 
                    xy=(most_common, energy_counts[most_common]),
                    xytext=(most_common+0.5, energy_counts[most_common]+1),
                    arrowprops=dict(arrowstyle='->', color='#D85A30'),
                    color='#D85A30', fontsize=9)
    
    trapped = sum(1 for e in final_energies if e > best_overall)
    fig.suptitle(
        f'HP lattice metastable state analysis — '
        f'{trapped}/{n_trials} runs trapped '
        f'({round(trapped/n_trials*100)}%)',
        fontsize=13, fontweight='bold', y=1.02
    )
    
    plt.tight_layout()
    plt.savefig('metastable_analysis.png', dpi=150, bbox_inches='tight')
    print(f"\nResults across {n_trials} trials:")
    print(f"Global minimum found: {best_overall}")
    print(f"Runs trapped in local minimum: {trapped}/{n_trials}")
    print(f"Trapping rate: {round(trapped/n_trials*100)}%")
    print(f"Final energy distribution: {energy_counts}")
def plot_small_multiples(sequence, chain, n_trials=20, steps=50000):
    
    results = []
    for i in range(n_trials):
        result = monte_carlo(sequence, chain, steps=steps)
        results.append(result)
        print(f"Trial {i+1}/{n_trials}: final energy = {result['final_energy']}")
    
    best_overall = min(r['best_energy'] for r in results)
    
    # 4 columns, enough rows to fit all trials
    n_cols = 4
    n_rows = -(-n_trials // n_cols)  # ceiling division
    
    fig, axes = plt.subplots(n_rows, n_cols, 
                             figsize=(16, n_rows * 3),
                             sharex=False, sharey=True)
    axes = axes.flatten()
    
    for i, (result, ax) in enumerate(zip(results, axes)):
        
        escaped = result['final_energy'] == best_overall
        color   = '#1D9E75' if escaped else '#7F77DD'
        border  = '#1D9E75' if escaped else '#7F77DD'
        
        # sample every 100 steps so lines are clean
        traj = result['trajectory'][::100]
        steps_x = list(range(0, len(result['trajectory']), 100))
        
        ax.plot(steps_x, traj, color=color, linewidth=1.2)
        ax.axhline(y=best_overall, color='#EF9F27', 
                  linestyle=':', linewidth=0.8, alpha=0.6)
        
        # color the border to immediately show escaped vs trapped
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(2)
        
        status = 'ESCAPED' if escaped else 'TRAPPED'
        ax.set_title(f'run {i+1}  |  {status}  |  final: {result["final_energy"]}', 
                    fontsize=8, color=color, fontweight='bold')
        ax.set_xlabel('step', fontsize=7)
        ax.set_ylabel('energy', fontsize=7)
        ax.tick_params(labelsize=7)
        
        # mark exactly where it got trapped
        if not escaped:
            final_e = result['final_energy']
            # find the step where it first reached final energy and stayed
            traj_full = result['trajectory']
            for step_idx in range(len(traj_full)-1, -1, -1):
                if traj_full[step_idx] != final_e:
                    trap_step = step_idx + 1
                    break
            else:
                trap_step = 0
            ax.axvline(x=trap_step, color='#D85A30', 
                      linestyle='--', linewidth=0.8, alpha=0.7)
            ax.text(trap_step, final_e + 0.2, 'trapped here', 
                   fontsize=6, color='#D85A30')
    
    # hide any unused panels
    for j in range(n_trials, len(axes)):
        axes[j].set_visible(False)
    
    escaped_count = sum(1 for r in results if r['final_energy'] == best_overall)
    fig.suptitle(
        f'Individual run trajectories — {escaped_count}/{n_trials} escaped to global min ({best_overall})\n'
        f'Orange dotted line = global minimum   |   Red dashed line = moment of trapping',
        fontsize=11, fontweight='bold', y=1.01
    )
    
    plt.tight_layout()
    plt.savefig('small_multiples.png', dpi=150, bbox_inches='tight')
    print(f"\nSaved to small_multiples.png")
    print(f"Escaped: {escaped_count}/{n_trials}")

sequence = "HPHPPHHPHPPHPHHPPHPH"
chain = [(i, 0) for i in range(len(sequence))]
plot_small_multiples(sequence, chain, n_trials=20, steps=50000)

sequence = "HPHPPHHPHPPHPHHPPHPH"
chain = [(i, 0) for i in range(len(sequence))]
plot_trajectories(sequence, chain, n_trials=20, steps=50000)
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
# print(is_valid([(0,0), (0,1), (1,1)]))  # should print True

# # Invalid - bead 0 and bead 2 share position
# print(is_valid([(0,0), (0,1), (0,0)]))  # should print False
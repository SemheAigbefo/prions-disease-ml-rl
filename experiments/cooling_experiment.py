import sys
sys.path.append('..')
from simulator.hp_lattice import monte_carlo
import matplotlib.pyplot as plt

sequence = "HPHPPHHPHPPHPHHPPHPH"
chain = [(i, 0) for i in range(len(sequence))]

# test four different cooling schedules
# same number of steps, different T_end
# slower cooling = T_end stays higher longer
cooling_schedules = [
    {"T_start": 2.0, "T_end": 0.5,  "label": "fast cool"},
    {"T_start": 2.0, "T_end": 0.1,  "label": "medium cool"},
    {"T_start": 2.0, "T_end": 0.01, "label": "slow cool"},
    {"T_start": 2.0, "T_end": 0.001,"label": "very slow cool"},
]

n_trials = 30
steps = 20000
results_summary = []

for schedule in cooling_schedules:
    print(f"\nRunning {schedule['label']}...")
    energies = []
    for t in range(n_trials):
        result = monte_carlo(
            sequence, chain,
            steps=steps,
            T_start=schedule['T_start'],
            T_end=schedule['T_end']
        )
        energies.append(result['final_energy'])
        print(f"  trial {t+1}/{n_trials}: {result['final_energy']}")

    best = min(energies)
    escaped = sum(1 for e in energies if e == best)
    avg = sum(energies) / len(energies)

    results_summary.append({
        'label': schedule['label'],
        'T_end': schedule['T_end'],
        'escape_rate': escaped/n_trials,
        'avg_energy': avg,
        'best': best,
        'distribution': energies
    })

    print(f"  escape rate: {escaped}/{n_trials} ({round(escaped/n_trials*100)}%)")
    print(f"  average final energy: {round(avg, 2)}")

# plot results
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# left: escape rate by cooling schedule
labels = [r['label'] for r in results_summary]
escape_rates = [r['escape_rate']*100 for r in results_summary]
colors = ['#D85A30', '#EF9F27', '#7F77DD', '#1D9E75']

axes[0].bar(labels, escape_rates, color=colors, 
            edgecolor='white', linewidth=0.5)
axes[0].set_ylabel('escape rate (%)')
axes[0].set_title('barrier crossing rate vs cooling schedule\n'
                  'does slower cooling help escape metastable traps?')
axes[0].set_ylim(0, 100)
for i, v in enumerate(escape_rates):
    axes[0].text(i, v+1, f'{round(v)}%', 
                ha='center', fontsize=10, fontweight='bold')

# right: energy distribution per schedule as box plot
dist_data = [r['distribution'] for r in results_summary]
bp = axes[1].boxplot(dist_data, labels=labels, patch_artist=True)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[1].set_ylabel('final energy')
axes[1].set_title('energy distribution per cooling schedule\n'
                  'lower = better, less spread = more consistent')
axes[1].axhline(y=min(r['best'] for r in results_summary),
               color='#1D9E75', linestyle=':', linewidth=1)

fig.suptitle(
    f'Cooling schedule experiment — {n_trials} trials each, {steps} steps\n'
    f'Sequence: {sequence}',
    fontsize=11, fontweight='bold'
)
plt.tight_layout()
plt.savefig('../results/figures/cooling_experiment.png', 
            dpi=150, bbox_inches='tight')
plt.show()
print("\nSaved to results/figures/cooling_experiment.png")
from programming1 import generate_system, timed_report
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Configuration
dimensions = range(1, 51)  # 1 to 50
num_samples_random = 5  # 5 samples for well-conditioned and ill-conditioned
num_samples_exact = 1   # 1 sample for exact
m_timing = 10000

# Fields to plot (from postprocessing + runtime_sec)
fields = [
    "pivot_growth",
    "cond_estimate",
    "cond_to_tcond",
    "ferr",
    "ferr_to_terr",
    "terr_to_cond_macheps",
    "scaled_berr",
    "berr_to_macheps",
    "runtime_sec"
]

# Field display names
field_names = {
    "pivot_growth": "Pivot Growth Factor",
    "cond_estimate": "Condition Number Estimate (1/RCOND)",
    "cond_to_tcond": "Cond Estimate / True Cond",
    "ferr": "Forward Error Bound (FERR)",
    "ferr_to_terr": "FERR / True Error",
    "terr_to_cond_macheps": "True Error / (macheps * cond)",
    "scaled_berr": "Scaled Backward Error",
    "berr_to_macheps": "BERR / macheps",
    "runtime_sec": "Runtime (seconds)"
}

print("="*80)
print("Generating data for visualization...")
print("="*80)

# Data structure: data[case_type][field][dimension] = {'gepp': [...], 'gecp': [...]}
data = {
    'well_conditioned': {field: {n: {'gepp': [], 'gecp': []} for n in dimensions} for field in fields},
    'ill_conditioned': {field: {n: {'gepp': [], 'gecp': []} for n in dimensions} for field in fields},
    'exact': {field: {n: {'gepp': [], 'gecp': []} for n in dimensions} for field in fields}
}

# Collect data
for n in dimensions:
    print(f"\nProcessing dimension n={n}...")

    # Well-conditioned systems (5 samples)
    for seed in range(num_samples_random):
        try:
            _, A, x, b, c = generate_system(n, well_conditioned=True, exact=False, seed=seed)
            results = timed_report(n, A, x, b, c, m=m_timing)

            for field in fields:
                data['well_conditioned'][field][n]['gepp'].append(results['gepp'][field])
                data['well_conditioned'][field][n]['gecp'].append(results['gecp'][field])
        except Exception as e:
            print(f"  Warning: Well-conditioned seed={seed} failed: {e}")

    # Ill-conditioned systems (5 samples)
    for seed in range(num_samples_random):
        try:
            _, A, x, b, c = generate_system(n, well_conditioned=False, exact=False, seed=seed)
            results = timed_report(n, A, x, b, c, m=m_timing)

            for field in fields:
                data['ill_conditioned'][field][n]['gepp'].append(results['gepp'][field])
                data['ill_conditioned'][field][n]['gecp'].append(results['gecp'][field])
        except Exception as e:
            print(f"  Warning: Ill-conditioned seed={seed} failed: {e}")

    # Exact system (1 sample)
    try:
        _, A, x, b, c = generate_system(n, well_conditioned=True, exact=True, seed=0)
        results = timed_report(n, A, x, b, c, m=m_timing)

        for field in fields:
            data['exact'][field][n]['gepp'].append(results['gepp'][field])
            data['exact'][field][n]['gecp'].append(results['gecp'][field])
    except Exception as e:
        print(f"  Warning: Exact case failed: {e}")

print("\n" + "="*80)
print("Creating plots...")
print("="*80)

# Create plots
for field in fields:
    print(f"\nPlotting {field}...")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(field_names[field], fontsize=16)

    cases = [
        ('well_conditioned', 'Well-Conditioned Systems', axes[0]),
        ('ill_conditioned', 'Ill-Conditioned Systems', axes[1]),
        ('exact', 'Exact Matrix', axes[2])
    ]

    for case_type, case_name, ax in cases:
        # Extract data for this case
        dims = []
        gepp_vals = []
        gecp_vals = []

        for n in dimensions:
            gepp_data = data[case_type][field][n]['gepp']
            gecp_data = data[case_type][field][n]['gecp']

            if gepp_data and gecp_data:  # Only plot if data exists
                for val in gepp_data:
                    dims.append(n)
                    gepp_vals.append(val)

                for val in gecp_data:
                    gecp_vals.append(val)

        # Create separate x coordinates for gepp and gecp to avoid overlap
        if dims:
            dims_array = np.array(dims)
            gepp_vals_array = np.array(gepp_vals)
            gecp_vals_array = np.array(gecp_vals)

            # Filter out infinite values for better visualization
            finite_gepp = np.isfinite(gepp_vals_array)
            finite_gecp = np.isfinite(gecp_vals_array)

            # Plot GEPP with X markers
            ax.scatter(dims_array[finite_gepp], gepp_vals_array[finite_gepp],
                      marker='x', s=50, c='blue', label='GEPP', alpha=0.6)

            # Plot GECP with O markers
            ax.scatter(dims_array[finite_gecp], gecp_vals_array[finite_gecp],
                      marker='o', s=50, c='red', facecolors='none', label='GECP', alpha=0.6)

        ax.set_xlabel('Dimension (n)')
        ax.set_ylabel(field_names[field])
        ax.set_title(case_name)
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Use log scale for y-axis if values span multiple orders of magnitude
        if dims:
            vals_all = np.concatenate([gepp_vals_array[finite_gepp], gecp_vals_array[finite_gecp]])
            if len(vals_all) > 0 and np.max(vals_all) / np.min(vals_all) > 100:
                ax.set_yscale('log')

    plt.tight_layout()
    filename = f'data_collection/plot_{field}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")

print("\n" + "="*80)
print("Visualization complete! Generated 9 plot files.")
print("="*80)

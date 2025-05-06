# -*- coding: utf-8 -*-
"""
Simplified SiZer Map Example (Chaudhuri and Marron, 1999):
- Blue = Increasing (1)
- Red = Decreasing (2)
- Gray = Inconclusive (3)
- Purple = Insufficient Data (4)
This code provides an illustrative example of a SiZer map for educational purposes.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter

# --- Step 1: Generate synthetic data ---
np.random.seed(42)
t = np.linspace(0, 100, 500)
data = np.sin(t / 10) + 0.5 * np.sin(t / 5) + np.random.normal(0, 0.1, len(t))

# --- Step 2: Define smoothing scales ---
scales = np.logspace(np.log10(2.5), np.log10(20), 10)

# --- Step 3: Compute SiZer map ---
n_boot = 200
sizer_map = np.zeros((len(scales), len(t)))

for i, h in enumerate(scales):
    smoothed = gaussian_filter1d(data, sigma=h, mode='reflect')
    deriv = np.gradient(smoothed, t)

    residuals = data - smoothed
    boot_slopes = np.array([
        np.gradient(
            gaussian_filter1d(smoothed + np.random.choice(residuals, len(data), replace=True), sigma=h, mode='reflect'),
            t
        )
        for _ in range(n_boot)
    ])

    ci_lower = np.percentile(boot_slopes, 2.5, axis=0)
    ci_upper = np.percentile(boot_slopes, 97.5, axis=0)

    # Adjusted ESS threshold for visible purple areas
    ess = h * len(t) / (t[-1] - t[0])
    sizer_map[i, :] = np.select(
        [ess < 15, ci_lower > 0, ci_upper < 0, True],
        [4, 1, 2, 3]
    )

# --- Step 4: Plotting ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [1, 2]})

# Top: Noisy signal with smoothed trend
ax1.plot(t, data, 'k.', alpha=0.3, label='Noisy data')
ax1.plot(t, gaussian_filter1d(data, sigma=5, mode='reflect'), 'r-', lw=2, label='Smoothed trend')
ax1.set_title("Time Series with Noise", fontsize=14)
ax1.legend()

# Bottom: SiZer map
colors = ['blue', 'red', 'lightgray', 'purple']
cmap = ListedColormap(colors)
im = ax2.imshow(
    sizer_map,
    aspect='auto',
    extent=[t[0], t[-1], scales[-1], scales[0]],
    cmap=cmap,
    vmin=1,
    vmax=4,
    interpolation='nearest'
)
ax2.set_yscale('log')
ax2.set_ylabel('Scale (Bandwidth)', fontsize=12)
ax2.set_xlabel('Time (t)', fontsize=12)
ax2.set_title('SiZer Map: Trends Across Scales', fontsize=14)

# Y-axis ticks: fewer, plain formatting (no scientific notation)
tick_locations = scales[::-1][::2]
tick_labels = [f'{s:.1f}' for s in tick_locations]
ax2.set_yticks(tick_locations)

# Use a custom formatter to bypass Matplotlib's logarithmic formatting
label_dict = dict(zip(tick_locations, tick_labels))
def fixed_label_formatter(x, pos):
    return label_dict.get(x, f'{x:.1f}')
ax2.yaxis.set_major_formatter(FuncFormatter(fixed_label_formatter))

# Disable minor ticks to prevent any additional formatting
ax2.yaxis.set_minor_locator(plt.NullLocator())

ax2.tick_params(axis='y', which='major', labelsize=10)

# Colorbar (legend on the right)
cbar = fig.colorbar(im, ax=ax2, ticks=[1.5, 2.5, 3.5, 4.5])
cbar.ax.set_yticklabels([
    'Increasing (Blue)',
    'Decreasing (Red)',
    'Inconclusive (Gray)',
    'Insufficient Data (Purple)'
])
cbar.ax.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig('sizer_map_final.png', dpi=300, bbox_inches='tight')
plt.show()
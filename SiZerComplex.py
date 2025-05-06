import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from matplotlib.colors import ListedColormap

# --- Step 1: Generate synthetic data ---
np.random.seed(88)
t = np.linspace(0, 100, 500)
data = np.sin(t / 10) + 0.5 * np.sin(t / 5) + np.random.normal(0, 0.1, len(t))

# --- Step 2: Define scales (bandwidths) ---
scales = np.logspace(np.log10(1), np.log10(30), 12)  # Finer range for better visualization

# --- Step 3: Compute SiZer map ---
n_boot = 200
sizer_map = np.zeros((len(scales), len(t)))

for i, h in enumerate(scales):
    smoothed = gaussian_filter1d(data, sigma=h, mode='reflect')
    window_length = max(7, min(int(2 * h / (t[1] - t[0])) // 2 * 2 + 1, 51))
    slope = savgol_filter(smoothed, window_length=window_length, polyorder=3, deriv=1, delta=t[1]-t[0])
    residuals = data - smoothed
    boot_slopes = np.array([
        savgol_filter(
            gaussian_filter1d(smoothed + np.random.choice(residuals, len(data), replace=True), sigma=h, mode='reflect'),
            window_length=window_length, polyorder=3, deriv=1, delta=t[1]-t[0]
        )
        for _ in range(n_boot)
    ])
    ci_lower = np.percentile(boot_slopes, 2.5, axis=0)
    ci_upper = np.percentile(boot_slopes, 97.5, axis=0)
    ess = h * len(t) / (t[-1] - t[0])
    sizer_map[i, :] = np.select(
        [(ess < 5), (ci_lower > 0), (ci_upper < 0), True],
        [4, 1, 2, 3]
    )

# --- Step 4: Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [1, 2]})

# Top: Smoothed data
ax1.plot(t, data, 'k.', alpha=0.3, label='Noisy data')
ax1.plot(t, gaussian_filter1d(data, sigma=scales[len(scales)//2], mode='reflect'), 'r-', lw=2, label='Smoothed trend')
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

# Y-tick labels
tick_indices = np.linspace(0, len(scales) - 1, 6, dtype=int)
tick_locations = scales[tick_indices][::-1]
tick_labels = [f'{s:.1f}' for s in tick_locations]
ax2.set_yticks(tick_locations)
ax2.set_yticklabels(tick_labels)

# Colorbar
cbar = fig.colorbar(im, ax=ax2, ticks=[1, 2, 3, 4])
cbar.set_ticklabels(['Increasing (Blue)', 'Decreasing (Red)', 'Inconclusive (Gray)', 'Insufficient Data (Purple)'])

plt.tight_layout()

# Save before showing
plt.savefig('sizer_map_example.png', dpi=300, bbox_inches='tight')
plt.show()

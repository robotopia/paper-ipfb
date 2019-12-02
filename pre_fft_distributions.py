import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Generate random "5-bit" values
N = 1000
X = np.random.randn(N,1536)*2
X5 = np.round(X)

# Get the filter coefficients
coeffs = np.loadtxt("r.txt")

# Perform the window-overlap-add
filtered = np.reshape(X5*coeffs, (N*1536//12, 12)).transpose()
polyphased = np.sum(filtered, axis=0)

# Demote the result to 8 bits using the weird rounding scheme
reduced = polyphased/2**14
pos_mask = (reduced > 0)
demoted_symmetric = np.round(reduced)
demoted = demoted_symmetric*pos_mask + np.floor(reduced)*np.logical_not(pos_mask)
#diff = demoted_symmetric - demoted
demoted_all = np.array([demoted_symmetric, demoted]).transpose()

# Calculate (and print out) some means
orig_mean = np.mean(X5)
sym_mean  = np.mean(demoted_symmetric)
asym_mean = np.mean(demoted)
print("Mean of original 5-bit samples: {0}".format(orig_mean))
print("Mean of  symmetrically-rounded 8-bit samples: {0}".format(sym_mean))
print("Mean of asymmetrically-rounded 8-bit samples: {0}".format(asym_mean))

# Make plots
fig, ax = plt.subplots(2,1)
fig.set_size_inches([5.0, 7.0])

# Create a small example to demostrate the rounding scheme used in the FPGAs
x0 = np.arange(-5.5, 5.5)
x0 = np.array([x0, x0+1])
y  = np.arange(-5, 6)
y  = np.array([y, y])

y_offset = 0.2
y0 = y + 0.2
y1 = y - 0.2

x1 = np.array([[-5, -4, -3, -2, -1, 0,   0.5, 1.5, 2.5, 3.5, 4.5],
               [-4, -3, -2, -1,  0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]])

x0_open   = np.arange(-5.5, 6.5)
y0_open   = np.array([-5, -4, -3, -2, -1, 0, 0, 1, 2, 3, 4, 5]) + y_offset
x0_closed = np.arange(-4.5, 5.5)
y0_closed = np.array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) + y_offset

x1_open   = np.array([-4, -3, -2, -1,  0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
y1_open   = np.arange(-5, 6) - y_offset
x1_closed = np.array([-5, -4, -3, -2, -1, 0, 0.5, 1.5, 2.5, 3.5, 4.5])
y1_closed = y1_open

ax[0].plot(x0, y0, color='C0', label="Symmetric rounding")
ax[0].plot(x0_open,   y0_open,   color='white', linestyle='None', marker="o", fillstyle='full')
ax[0].plot(x0_open,   y0_open,   color='C0',    linestyle='None', marker="o", fillstyle='none')
ax[0].plot(x0_closed, y0_closed, color='C0',    linestyle='None', marker="o", fillstyle='full')

ax[0].plot(x1, y1, color='C1', label="Asymmetric rounding")
ax[0].plot(x1_open,   y1_open,   color='white', linestyle='None', marker="o", fillstyle='full')
ax[0].plot(x1_open,   y1_open,   color='C1',    linestyle='None', marker="o", fillstyle='none')
ax[0].plot(x1_closed, y1_closed, color='C1',    linestyle='None', marker="o", fillstyle='full')

ax[0].axvline(0,         color='C0', linestyle='--')
ax[0].axvline(asym_mean, color='C1', linestyle='--')

ax[0].set_xticks(np.arange(-5,6))
ax[0].set_yticks(np.arange(-5,6))
ax[0].grid()

ax[0].set_xlim((-5.5,5.5))
ax[0].set_xlabel("$n/2^{14}$")
ax[0].set_ylabel("Rounded value")

# Create custom legend
custom_lines = [Line2D([0], [0], color='C0', lw=4),
                Line2D([0], [0], color='C1', lw=4)]
ax[0].legend(custom_lines, ["Symmetric rounding", "FPGA rounding"], loc='upper left')

# Plot the histogram
h = ax[1].hist(demoted_all, range=(-5.5,5.5), weights=(1/demoted.size)*np.ones(demoted_all.shape), bins=11, rwidth=0.8, log=True)
mistakes = np.logical_and(demoted[:] > 0, demoted[:] != demoted_symmetric[:])

ax[1].set_xticks(np.arange(-5,6))
ax[1].set_xlim((-5.5,5.5))
ax[1].set_ylim((0.01,0.5))
ax[1].set_xlabel("8-bit pre-FFT value")
ax[1].set_ylabel("Normalised frequency")

plt.savefig('pre_fft_distributions.eps', bbox_inches='tight', pad_inches=0.0)

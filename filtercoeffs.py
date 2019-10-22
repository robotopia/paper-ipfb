from scipy import hanning
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, figsize=(5,7))

r = np.loadtxt("r.txt")
rmax = max(r)
print("max coefficient = {0}".format(rmax))

P = 12 # number of taps
N = 128 # number of output channels
M = N*P

m = np.linspace(0, M, M, endpoint=False)
x = (m + 1 - M/2)/N

s = np.sinc(x)
hn = hanning(M+1)[1:]
r1 = s*hn*(rmax-0.125)
rounded = np.round(r1)

'''
# Plots to verify the constructed coefficients vs the actual coefficients.
axs[0].plot(r, label="True coefficients")
axs[0].plot(r1, label="Hanning * sinc")

residuals = r1 - r
axs[1].plot(residuals, label="Residuals")
axs[1].plot(rounded - r, label="Rounded residuals")
print("rounded residuals: min = {0};   max = {1}".format(min(r1-r), max(r1-r)))

axs[0].legend()
axs[1].legend()
plt.show()
'''

# Plots to go in the figure
axs[0].plot(m, hn, 'k:', label="Hanning window")
axs[0].plot(m, s, 'k--', label="Sinc function")
axs[0].plot(m, s*hn, 'k-', label="Analysis filter")

axs[0].set_xlabel("$m$")
axs[0].set_ylabel("$h[m]$")
axs[0].set_xlim([0,1536])

ctr = M//2 - 1
k = (np.arange(M) - ctr)/P * 10 # i.e. 10 kHz
fftr = np.roll(np.abs(np.fft.fft(s*hn)), ctr)
dB = 10*np.log10(fftr/max(fftr))
axs[1].plot(k, dB, 'k-')

axs[1].set_xlabel("$\\nu$ (kHz)")
axs[1].set_ylabel("Normalised power (dB)")
axs[1].set_xlim([-30,30])
axs[1].set_ylim([-60, 5])
plt.grid(True)

axs[0].legend(loc='upper right')
#plt.show()
plt.tight_layout()
plt.savefig("filter.eps")

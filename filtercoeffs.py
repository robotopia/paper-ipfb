import numpy as np
import matplotlib.pyplot as plt

r = np.loadtxt("r.txt")
rmax = max(r)
print("max coefficient = {0}".format(rmax))

def create_fine_PFB_filter(nchans, ntaps):

    P = ntaps # number of taps
    K = nchans # number of output channels
    N = nchans*ntaps

    n = np.linspace(0, N, N, endpoint=False)
    x = (n + 1 - N/2)/nchans

    sinc = np.sinc(x)
    hn   = np.hanning(N+1)[1:]
    h    = sinc*hn
    return h, n, hn, sinc

P = 12 # number of taps
K = 128 # number of output channels
N = K*P # size of filter

h, n, hn, s = create_fine_PFB_filter(K, P)

r1 = h*(rmax-0.125)
rounded = np.round(r1)

fig, axs = plt.subplots(2, figsize=(5,7))

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
axs[0].plot(n, hn, 'k:', label="Hanning window")
axs[0].plot(n, s, 'k--', label="Sinc function")
axs[0].plot(n, s*hn, 'k-', label="Analysis filter")

axs[0].set_xlabel("$n$")
axs[0].set_ylabel("$h[n]$")
axs[0].set_xlim([0,1536])

ctr = N//2 - 1
f = (np.arange(N) - ctr)/P * 10 # i.e. 10 kHz
fft = np.fft.fft(s*hn)
fftr = np.roll(np.abs(fft), ctr)
dB = 10*np.log10(fftr/max(fftr))
axs[1].plot(f, dB, 'k-')
for fshift in [-10, 10]:
    axs[1].plot(f + fshift, dB, '--', color='#c0c0c0')

axs[1].set_xlabel("$\\nu$ (kHz)")
axs[1].set_ylabel("Normalised power (dB)")
axs[1].set_xlim([-30,30])
axs[1].set_ylim([-60, 5])
plt.grid(True)

'''
# Plot FFT of least squares synthesis filter
F = np.loadtxt("rinv.txt")
print(F)
Ffft = np.fft.fft(F)
Ffftr = np.roll(np.abs(Ffft), ctr)
FdB = 10*np.log10(Ffftr/max(Ffftr))
axs[1].plot(f, FdB, 'b-')
'''

axs[0].legend(loc='upper right')
#plt.show()
plt.tight_layout()
plt.savefig("filter.eps", bbox_inches='tight')


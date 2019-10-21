from scipy import hanning
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2)

r = np.loadtxt("r.txt")
rmax = max(r)
print("max coefficient = {0}".format(rmax))

P = 12 # number of taps
N = 128 # number of output channels
M = N*P

x = (np.linspace(0, M, M, endpoint=False) - M/2 + 1)/N

s = np.sinc(x)
hn = hanning(M+1)[1:]
r1 = s*hn*(rmax-0.125)
rounded = np.round(r1)

axs[0].plot(r, label="True coefficients")
axs[0].plot(r1, label="Hanning * sinc")

residuals = r1 - r
axs[1].plot(residuals, label="Residuals")
axs[1].plot(rounded - r, label="Rounded residuals")
print("rounded residuals: min = {0};   max = {1}".format(min(r1-r), max(r1-r)))

axs[0].legend()
axs[1].legend()
plt.show()

from scipy import kaiser, hanning, hamming, bartlett, blackman
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2)

r = np.loadtxt("r.txt")
rnorm = r/max(r)

axs[0].plot(rnorm, label="True coefficients, normalised")

ntaps = 12
nsamps = 128
N = ntaps*nsamps
a = 6
k = kaiser(N, a)
axs[0].plot(kaiser(N,a), label="Kaiser {0}".format(a))

hn = hanning(N)
hm = hamming(N)
br = bartlett(N)
bl = blackman(N)

x = (np.linspace(0, N, N, endpoint=False) - N/2 + 1)/N
s = np.sinc(x*12)
r1 = s*k
axs[0].plot(r1, label="Kaiser * sinc")
axs[0].plot(rnorm/s, label="rnorm/sinc")

axs[1].plot(r1 - rnorm, label="Residuals")

axs[0].legend()
axs[1].legend()
plt.show()

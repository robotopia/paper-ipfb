import numpy as np
import matplotlib.pyplot as plt

# Parameters for PFB
ntaps = 12
K = 128
N = ntaps*K
M = K  # meaning "ciritically sampled"

# Choose the analysis filter
def wH(n):
    return np.sin(np.pi*(n+1)/(N+1))**2

def ws(n):
    return np.sinc((n+1-N/2)/K)

def h(n):
    return wH(n)*ws(n)*np.logical_and(n >= 0, n < N)

# Choose the synthesis filter
def f(n):
    H = np.sum(h(ns))
    return h(-n)

# Calculate λ(s) from h[n] and f[n]
ns = np.arange(N)
ss = np.arange(-12,13)
ms = np.arange(M)
NS, SS, MS = np.meshgrid(ns, ss, ms)
lmda = np.sum(f(NS - MS*M) * h(MS*M - NS + SS*M), axis=2)

# Make plots
plt.figure(1)
ns = np.arange(-N,N)
plt.plot(ns, h(ns), label='Analysis filter')
plt.plot(ns,f(ns), label='Synthesis filter')
plt.legend()

plt.figure(2)
plt.imshow(lmda, aspect='auto', interpolation='None', origin='lower', extent=[-0.5,N+0.5,ss[0]-0.5,ss[-1]+0.5])
cb = plt.colorbar()
plt.xlabel("$n$")
plt.ylabel("$s$")
cb.set_label("$\\lambda(s)$")
#plt.savefig("inverse_condition.eps", bbox_inches='tight')

plt.show()


'''
# Analysis-synthesis response of an impulse
nsamples = K*64
x_ntaps  = nsamples//K - ntaps
x = np.zeros(nsamples)
x[nsamples//3] = 1 # Put an impulse right smack bang in the middle
ns = np.arange(nsamples)
x = np.cos(2*ns)

X = np.array([[np.sum(h(m*M-ns) * x * np.exp(-2j*np.pi*k*ns/K)) for m in range(x_ntaps)] for k in range(K)])
X = np.roll(X, K//2, axis=0)

ks = np.arange(K)
ms = np.arange(x_ntaps)
NS, KS, MS = np.meshgrid(ns, ks, ms)
#xhat = [np.sum([f(n-m*M) / K * np.sum(X[:,m] * np.exp(2j*np.pi*ks*n/K)) for m in range(x_ntaps)]) for n in ns]
xhat = f(NS-MS*M) / K
print(xhat.shape)

fig, ax = plt.subplots(3,1)
im1 = ax[0].imshow(np.abs(X), aspect='auto')
im2 = ax[1].imshow(np.angle(X), cmap='hsv', aspect='auto')
im3 = ax[2].plot(x)
#plt.colorbar(im1)
#plt.colorbar(im2)
plt.show()
'''

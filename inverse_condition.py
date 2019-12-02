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
    return h(n)/H

# Calculate λ(s) from h[n] and f[n]
ns = np.arange(N)
ss = np.arange(-12,13)
lmda = np.array([sum(f(ns) * h(ns + s*M)) for s in ss])

# Analysis-synthesis response of an impulse
nsamples = K*64
x_ntaps  = nsamples//K - ntaps
x = np.zeros(nsamples)
x[nsamples//3] = 1 # Put an impulse right smack bang in the middle
ns = np.arange(nsamples)

X = np.array([[np.sum(h(m*M-ns) * x * np.exp(-2j*np.pi*k*ns/K)) for m in range(x_ntaps)] for k in range(K)])

#ks = np.arange(K)
#xhat = ([f(ns-m*M) / K * np.sum(X[:,m] * np.exp(2j*np.pi*ks*n/K)) for m in range(x_ntaps)])
#print(xhat.shape)

plt.figure(1)
plt.imshow(np.abs(X))
plt.colorbar()
plt.figure(2)
plt.imshow(np.angle(X), cmap='hsv')
plt.colorbar()
plt.show()

# Make plots
plt.plot(ss, lmda, '-xk')
plt.ylim([-0.1, 1.0])
plt.xlabel("$s$")
plt.ylabel("$\\lambda(s)$")
plt.savefig("inverse_condition.eps", bbox_inches='tight')

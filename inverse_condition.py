import numpy as np
import matplotlib.pyplot as plt

ntaps = 12
K = 128
N = ntaps*K
M = K  # meaning "ciritically sampled"

ns = np.arange(N)

def wH(n):
    return np.sin(np.pi*(n+1)/(N+1))**2

def ws(n):
    return np.sinc((n+1-N/2)/K)

def h(n):
    return wH(n)*ws(n)*np.logical_and(n >= 0, n < N)

# Choose the synthesis filter
def f(n):
    return h(n)

res = []

H = np.sum(h(ns))

ss = np.arange(-12,13)
res = np.array([sum((1/H) * f(ns) * h(ns + s*M)) for s in ss])

plt.plot(ss, res, '-x')
plt.xlabel("$s$")
plt.ylabel("$\\lambda(s)$")
plt.savefig("inverse_condition.eps", bbox_inches='tight')

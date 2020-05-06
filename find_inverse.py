import sys
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

# Get the least squared solution
Fhat = np.zeros((ntaps, M))
for n in range(M):
    # Try and find the inverse to the analysis filter, h(n)
    # Step 1: Create H array:
    #    [ h(n+11*M)      0          0         ...    ]
    #    [ h(n+10*M)  h(n+11*M)      0         ...    ]
    #    [ h(n+ 9*M)  h(n+10*M)  h(n+11*M)     ...    ]
    #    [    ...        ...        ...        ...    ]
    #    [    ...     h(n+ 0*M)  h(n+ 1*M)  h(n+ 2*M) ]
    #    [    ...         0      h(n+ 0*M)  h(n+ 1*M) ]
    #    [    ...         0          0      h(n+ 0*M) ]
    ncols = ntaps
    nrows = 2*ntaps - 1
    H = np.zeros((nrows, ncols))
    for i in range(nrows):
        if i < ntaps:
            H[i,0:(i+1)] = [h(n + m*M) for m in range(ntaps - i - 1, ntaps)]
        else:
            H[i,(i-ntaps+1):] = [h(n + m*M) for m in range(nrows - i)]

    # Now, we want to solve the matrix equation HF = D for F,
    # where D is the "Kronecker delta"
    D = np.zeros((nrows, 1))
    D[ntaps-1] = 1

    # The least squares solution is
    HT        = np.transpose(H)
    HTH       = np.matmul(HT, H)
    HTD       = np.matmul(HT, D)
    invHTH    = np.linalg.pinv(HTH)
    Fhat[:,n:n+1] = np.matmul(invHTH, HTD)

plt.figure("Least squares")
F = Fhat.flatten()
ns = np.arange(N)
plt.plot(ns, h(ns), label='Analysis filter')
plt.plot(ns, F, label="Least squares sol'n")
plt.legend()
plt.tight_layout()
plt.savefig("rinv.png")
#plt.show()

# Write out the spectrum to rinv.txt
header =  "Synthesis filter coefficients obtained by least squares solution to inverse analysis filter"
header += "Generated by\n"
header += "python" + " ".join(sys.argv)
np.savetxt("rinv.txt", F, header=header)


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

# Find some other synthesis filters!
filter_types = [24, 18, 12, "Mirror"]
plottable = [0, 1, 2, 3] # The indices for which filters in "filter_types" should be plotted
dashtypes = ['k:', 'k-.', 'k--', 'k-'] # These match up with the "filter_types" array
NSR = np.zeros((len(filter_types), K))

fig, axs = plt.subplots(len(plottable), sharex=True, gridspec_kw={'hspace': 0}) # Plots for the filter coefficients
fig_imp, axs_imp = plt.subplots(len(plottable), sharex=True, gridspec_kw={'hspace': 0})

for a in range(len(filter_types)):

    filter_size = filter_types[a]

    if not isinstance(filter_size, int):
        Ftaps = ntaps
    else:
        Ftaps = filter_size

    Hsize = (Ftaps + ntaps - 1, Ftaps)
    Href = (Hsize[0] // 2, Hsize[1] // 2 - ntaps // 2)

    Fhat = np.zeros((Hsize[1], K))
    Dhat = np.zeros((Hsize[0], K))

    for n in range(K):
        # Try and find the inverse to the analysis filter, h(n)
        # Step 1: Create H array, e.g.,
        #          <----------------  Hsize[1] ----------------->
        #    ^     [ h(n+11*M)      0          0         ...    ]
        #    |     [ h(n+10*M)  h(n+11*M)      0         ...    ]
        #    |     [ h(n+ 9*M)  h(n+10*M)  h(n+11*M)     ...    ]
        # Hsize[0] [    ...        ...        ...        ...    ]
        #    |     [    ...     h(n+ 0*M)  h(n+ 1*M)  h(n+ 2*M) ]
        #    |     [    ...         0      h(n+ 0*M)  h(n+ 1*M) ]
        #    v     [    ...         0          0      h(n+ 0*M) ]
        H = np.zeros(Hsize)
        Mmat = np.zeros(Hsize)
        for i in range(Hsize[0]):
            for j in range(Hsize[1]):
                m = (j - i) - (Href[1] - Href[0])
                #Mmat[i,j] = m
                H[i,j] = h(n + m*M)

        # Now, we want to solve the matrix equation HF = D for F,
        # where D is the "Kronecker delta"
        D = np.zeros((Hsize[0], 1))
        D[(Hsize[0] + 1)//2 - 1] = 1

        if filter_size == "Mirror":
            Fhat[:,n] = [h(n + m*M) for m in range(Ftaps)] # Set synth filter = analysis filter
        elif filter_size == "fft":
            Fhat[:,n] = 1/Ftaps
        else:
            # H is non-square
            HT        = np.transpose(H)
            HTH       = np.matmul(HT, H)
            HTD       = np.matmul(HT, D)
            invHTH    = np.linalg.pinv(HTH)
            Fhat[:,n:n+1] = np.matmul(invHTH, HTD)

        # Let's see how well the solution does
        Dhat[:,n:n+1] = np.matmul(H, Fhat[:,n:n+1])
        #impulse_response[:,n:n+1] = np.matmul...

    F = Fhat.flatten()
    Fns = np.arange(F.size)
    impulse_response = Dhat.flatten()**2 # Units of power (i.e. not dB)

    # Plot figures
    # Plots for the filter coefficients and the S/N in dB as a function of tap position
    ns = np.arange(len(impulse_response)) - (Ftaps-1)*M/2
    if a in plottable:
        if isinstance(filter_size, int):
            label = "Least squares filter"
        elif filter_size == "Mirror":
            label = "Mirror filter"
        else:
            label = "(undefined)"
        label += ", {0} taps".format(Ftaps)

        plot_no = plottable.index(a)
        axs[plot_no].plot(Fns - F.size//2, F, 'k-', label=label)
        axs[plot_no].legend()
        axs[plot_no].set_yticks([0,1])
        axs[plot_no].set_ylabel(" ")

        xlims = [-18.5,18.5]
        axs_imp[plot_no].plot(xlims,[-25,-25], color='gray', linestyle='dashed')
        axs_imp[plot_no].plot(ns/M - 5.5, 10*np.log10(impulse_response), 'k-', label=label)
        axs_imp[plot_no].legend()
        axs_imp[plot_no].set_yticks([-50,-25,0])
        axs_imp[plot_no].set_ylabel(" ")
        axs_imp[plot_no].set_ylim([-70,5])
        axs_imp[plot_no].set_xlim(xlims)

        plt.figure("Impulse response (zoom)")
        plt.plot(ns - 5.5*M, 10*np.log10(impulse_response), dashtypes[a], label=label)

    # Write out the impulse response to file
    header =  "Impulse response for the {0}-tap filter\n".format(filter_size)
    header += "Generated by\n"
    header += "python" + " ".join(sys.argv)
    np.savetxt("impulse_response_{0}.txt".format(filter_size), impulse_response, header=header)

    # Finally, calculate the relative S/N for each n
    NSR[a,:] = 1 / (Dhat[(Hsize[0] + 1)//2 - 1])**2 - 1 # noise-to-signal ratio

    # Write out the spectrum to rinv.txt
    header =  "Synthesis filter coefficients for the {0}-tap filter\n".format(filter_size)
    header += "Generated by\n"
    header += "python" + " ".join(sys.argv)
    np.savetxt("rinv_{0}.txt".format(filter_size), F, header=header)



    # Plot the ACTUAL S/N (i.e. measured using real data)
    # Only for MIRROR filter and LSQ12
    ms = np.arange(M)
    if Ftaps == 12:
        plt.figure("Measured S/N")
        if isinstance(filter_size, int):
            label = "LSQ12"
            measured = np.loadtxt("sample_error_LSQ12.txt")
            linecolor = 'k'
        elif filter_size == "Mirror":
            label = "MIRROR"
            measured = np.loadtxt("sample_error_MIRROR.txt")
            linecolor = 'r'
        else:
            label = "(undefined)"
        NSR_measured = (measured[:,1]/measured[:,0])**2
        plt.plot(ms, -10*np.log10(1 + NSR_measured), linecolor+'-', label=label)
        plt.plot(ms, -10*np.log10(1 + (NSR_measured - NSR[a,:])), linecolor+'--', label=label+" (filter-subtracted)")

axs_imp[-1].set_xlabel("Tap number")
fig_imp.text(0.01, 0.5, 'Power (dB)', va='center', rotation='vertical')
#fig_imp.set_size_inches([12.8,4.8])
fig_imp.tight_layout()
fig_imp.savefig("impulse_response.eps")

plt.figure("Impulse response (zoom)")
plt.xlabel("Position within tap, $n$")
plt.ylabel("Power (dB)")
plt.xlim([0,127])
plt.ylim([-1.65,0.05])
plt.tight_layout()
plt.legend()
plt.savefig("snr.eps")

plt.figure("Measured S/N")
plt.xlabel("Position within tap, $n$")
plt.ylabel("Reconstructed S/N (dB)")
plt.ylim([-0.68,0])
plt.tight_layout()
plt.legend()
plt.savefig("snr_measured.eps")

axs[-1].set_xlabel("Sample number, $n$")
fig.text(0.02, 0.5, 'Normalised filter coefficients', va='center', rotation='vertical')
fig.tight_layout()
fig.savefig("inverse_filters.eps")

'''
# Show figures
plt.legend()
plt.show()
'''

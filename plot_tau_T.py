"""Plots tau vs. T for various N"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt

newparams = {'axes.labelsize': 15,
             'axes.linewidth': 1,
             'lines.linewidth': 1.5, 
             'figure.figsize': (9, 8),
             'ytick.labelsize': 15,
             'xtick.labelsize': 15,
             'ytick.major.pad': 5,
             'xtick.major.pad': 5,
             'legend.fontsize': 15,
             'legend.frameon': True, 
             'legend.handlelength': 1.5,
             'axes.titlesize': 20,
             'mathtext.fontset': 'stix',
             'font.family': 'STIXGeneral'}
plt.rcParams.update(newparams)

Tc = 2.26919  # Onsager's critical temp

def onsager(T):
    return 2 - T*Tc*np.log(1/np.tanh(1/(T*Tc)))

options, files = getopt.gnu_getopt(sys.argv[1:], 'N:t:', ["title=", "save="])

print('\n'.join(files))

Ns = []
title = ""
save = False
for opt, arg in options:
    if opt == "-N":
        Ns = arg
    elif opt in ("-t", "--title"):
        title = arg
    elif opt in ("--save"):
        save = True
        figname = arg
if not Ns:
    Ns = input("Comma spearated list of N: ")

Ns = Ns.split(",")

markers = ["D", "s"]
for i, file in enumerate(files):
    T, tau, tau_std = np.loadtxt(file)
    plt.errorbar(T, tau, yerr=tau_std, label=f"N={Ns[i]}", fmt=markers[i])

T = np.linspace(0, 1, 100)
plt.plot(T, onsager(T), "--", label="Onsager")

if title:
    plt.title(title)
plt.xlabel("$T/T_c$")
plt.ylabel("$\\tau$")
plt.axhline(y=0, linestyle="dashed", color="gray", linewidth=0.2)

plt.axvline(x=1, linestyle="dashed", color="gray", linewidth=0.2)
plt.legend()
if save:
    plt.savefig(figname)
plt.show()

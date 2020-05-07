"""Plots Ntau vs. N at Tc"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt

# newparams = {'axes.labelsize': 15,
#              'axes.linewidth': 1,
#              'lines.linewidth': 1.5, 
#              'figure.figsize': (16, 8),
#              'ytick.labelsize': 15,
#              'xtick.labelsize': 15,
#              'ytick.major.pad': 5,
#              'xtick.major.pad': 5,
#              'legend.fontsize': 15,
#              'legend.frameon': True, 
#              'legend.handlelength': 1.5,
#              'axes.titlesize': 20,
#              'mathtext.fontset': 'stix',
#              'font.family': 'STIXGeneral'}
# plt.rcParams.update(newparams)
plt.style.use("bmh")

options, files = getopt.gnu_getopt(sys.argv[1:], 'f:t:', ["title=", "file=", "log", "save="])

title = ""
log = False
save = False
for opt, arg in options:
    if opt in ("-t", "--title"):
        title = arg
    elif opt in ("-f", "--file"):
        file = arg
    elif opt == "--log":
        log = True
    elif opt == "--save":
        save = True
        figname = arg

try:
    file
except NameError as e:
    print("File not defined!")
    print("Usage: [(-t|--title)=<title>]  [(-f|--file)=<file>] [--log] [--save=<figfile>]")
    exit()

N, tau, tau_std = np.loadtxt(file)
plt.errorbar(N, tau, yerr=tau_std, fmt="o", capsize=2)

if title:
    plt.title(title)
if log:
    plt.xscale("log")
    plt.yscale("log")
plt.xlabel("$N$")
plt.ylabel("$N\\tau$")
if save:
    plt.savefig(figname)
plt.show()

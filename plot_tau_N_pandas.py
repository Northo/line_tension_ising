import matplotlib.pyplot as plt
import numpy as np
import pandas
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

def usage():
    print("Usage: ")
    print("--data=<datafile> [--title=<title>] [--save=<figfile>] [--logy] [--Nmax/Nmin=<Nmax/Nmin>]")
    exit()

options, rest = getopt.gnu_getopt(sys.argv[1:], 't:', [
    "title=",
    "save=",
    "logy",
    "data=",
    "Nmin=",
    "Nmax=",
])


title = ""
ts = ""
save = False
N_min = -1
N_max = -1
filename = ""
for opt, arg in options:
    if opt == "--title":
        title = arg
    elif opt in ("--save"):
        save = True
        figname = arg
    elif opt == "--logy":
        plt.yscale("log")
    elif opt == "--data":
        filename = arg
    elif opt == "--Nmin":
        N_min = int(arg)
    elif opt == "--Nmax":
        N_max = int(arg)

if not filename:
    print("Missing filename! ")
    usage()


df = pandas.read_csv(filename, delimiter="\t")
# Remove values over critica

if N_min != -1:
    df = df[df["N_sweeps"] >= N_min]
if N_max != -1:
    df = df[df["N_sweeps"] <= N_max]


systems = df.groupby(["system"])
for system_name, system in systems:
    plt.errorbar(
        system["N"],
        system["N"]*system["tau"],
        yerr=system["N"]*system["tau_std"],
        fmt="o",
        label=system_name,
    )

if title:
    plt.title(title)
plt.legend()
plt.xlabel("N")
plt.ylabel("$N\\tau$")

if save:
    plt.savefig(figname)
plt.show()

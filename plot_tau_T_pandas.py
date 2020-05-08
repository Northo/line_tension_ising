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
    print("--data=<datafile> [--title=<title>] [--save=<figfile>] [[--system=<torus|pp>]]")
    exit()

options, rest = getopt.gnu_getopt(sys.argv[1:], 't:', [
    "title=",
    "save=",
    "data=",
    "system=",
])


title = ""
ts = ""
save = False
N_min = -1
N_max = -1
filename = ""
ylim = ""
system = ""
for opt, arg in options:
    if opt == "--title":
        title = arg
    elif opt in ("--save"):
        save = True
        figname = arg
    elif opt == "--system":
        system = arg
    elif opt == "--data":
        filename = arg

if not filename:
    print("Missing filename! ")
    usage()

df = pandas.read_csv(filename, delimiter="\t")
# Remove values over critica
if system:
    df = df[df["system"] == system]

lattice_numbers = df.groupby(["N"])
i = 0
markers = ["D", "s", "^", "v", "o"]
for lattice_number, df_lattice in lattice_numbers:
    plt.errorbar(
        df_lattice["T"]/Tc,
        df_lattice["tau"],
        df_lattice["tau_std"],
        fmt=markers[i],
        label=f"N = {lattice_number}",
    )
    i += 1

if title:
    plt.title(title)

T = np.linspace(0, 1, 100)
plt.plot(T, onsager(T), "--", label="Onsager")
plt.legend()
plt.xlabel("$T/T_c$")
plt.ylabel("$\\tau$")

if save:
    plt.savefig(figname)
plt.show()

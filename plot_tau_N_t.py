"""Plots tau vs. T for various N"""

import numpy as np
import pandas
import matplotlib.pyplot as plt
import sys
import getopt
from scipy.optimize import curve_fit
from scipy import stats

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
filename = "datadir/T_N_fresh_run.dat"

def onsager(T):
    return 2 - T*Tc*np.log(1/np.tanh(1/(T*Tc)))

def t(T):
    return 1 - T/Tc

options, rest = getopt.gnu_getopt(sys.argv[1:], 't:', [
    "title=",
    "save=",
    "logy",
    "data=",
    "Nmin=",
    "Nmax=",
    "system=",
])


title = ""
ts = ""
save = False
N_min = -1
N_max = -1
system = ""
for opt, arg in options:
    if opt == "--title":
        title = arg
    elif opt in ("--save"):
        save = True
        figname = arg
    elif opt == "--logy":
        plt.yscale("log")
    elif opt == "-t":
        ts = arg
    elif opt == "--data":
        filename = arg
    elif opt == "--Nmin":
        N_min = int(arg)
    elif opt == "--Nmax":
        N_max = int(arg)
    elif opt == "--system":
        system = arg


df = pandas.read_csv(filename, delimiter="\t")
df["t"] = t(df["T"])
# Remove values over critica
df = df[df["t"]>0]
if system:
    df = df[df["system"]==system]

if N_min != -1:
    df = df[df["N_sweeps"] >= N_min]
if N_max != -1:
    df = df[df["N_sweeps"] <= N_max]

def plot(df, marker="s", label=""):
    y = df["tau"]/df["t"]
    yerr = df["tau_std"]/df["t"]
    x = 1/(df["N"]*df["t"])
    plt.errorbar(x, y, yerr=yerr, label=label, fmt=marker)

# If given ts, group
i = 0  # For markers
markers = ["D", "o", "v", "s", "^", "x", "+", "<", ">"]
if ts:
    ts = ts.split(",")
    ts = [float(t) for t in ts]
    groups = df.groupby(["t"])
    newDf = pandas.DataFrame()
    for t, group in groups:
        # Hack to avoid round off trouble
        t = float(f"{t:.4f}")
        if t in ts:
            newDf = newDf.append(group)
            plot(group, marker=markers[i], label=t)
            i+=1
    df = newDf
else:
    plot(df)

if title:
    plt.title(title)

x_data = 1/(df["N"]*df["t"])
y_data = df["tau"]/df["t"]

slope, intercept, r, p, std_err = stats.linregress(x_data, y_data)

curve_x = np.array([0, np.max(x_data)])
plt.plot(curve_x, curve_x*slope + intercept, label=f"$\\tau_0 = {intercept:.2f}\\pm{std_err:.2f}$")
plt.xlabel("$1/Nt$")
plt.ylabel("$\\frac{\\tau}{t}$")
plt.axhline(y=0, linestyle="dashed", color="gray", linewidth=0.2)

plt.axvline(x=0, linestyle="dashed", color="gray", linewidth=0.2)
plt.legend()
if save:
    plt.savefig(figname)
plt.show()

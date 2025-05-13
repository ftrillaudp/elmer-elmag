##################
### Python 3.x ###
##################

# AUTHOR: F. Trillaud <ftrillaudp@gmail.com>
# DATE: 01/25/2024

# Libraries:
import math
import itertools
import csv
import glob
import seaborn as sns
import scipy as sp
import numpy as np
# Functions
import matplotlib.pyplot as plt
import re
from sklearn.metrics import r2_score
from termcolor import colored

# For windows:
myPath_resu = './resu'
file = myPath_resu+'/circuit.dat'

### Figures #######
# ~ plt.rcParams.update({"text.usetex": True, "font.family": "DejaVu Sans", "font.sans-serif": ["Helvetica"]})
#plt.rcParams.update({"text.usetex": True, "font.family": "DejaVu Sans", "font.sans-serif": ["Avant Garde"]})
#plt.rcParams.update({"text.usetex": True, "font.family": "DejaVu Sans", "font.sans-serif": ["Computer Modern Sans serif"]})
## for Palatino and other serif fonts use:
#plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Palatino"],})
#plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],})

plt.rc('xtick', labelsize = 12)    # fontsize of the tick labels
plt.rc('ytick', labelsize = 12)    # fontsize of the tick labels
plt.grid(False)

markers = itertools.cycle(('.', '+', 's', 'o', '*', 'v', 'H', 'D'))
palette = itertools.cycle(sns.color_palette())

# plot with various axes scales
#plt.figure(1, figsize = (3,6))
fig, axs = plt.subplots(1, figsize=(8, 8))
fig.tight_layout()
axs.set_xlabel(r'$t$ (s)', fontsize = 16)
axs.set_ylabel(r'$i$ (A)', fontsize = 16)

# ~ plt.ylim((2.9, 3.1))

X, Y = list(), list()
### Open external file ###
with open(file) as myfile:
    mydata = csv.reader(myfile, delimiter = ' ')
    for row in mydata:
        row = list(filter(None, row))
        #print(row)
        X.append(float(row[2]))
        Y.append(abs(float(row[5])))


axs.plot(X, Y, color = 'red', linestyle='-', linewidth=2)
# ~ axs.legend(fontsize=12, loc='upper right')


plt.savefig('./resu/figure.png', format = 'png', dpi = 75, bbox_inches = 'tight')
# ~ plt.show()

from __future__ import division, print_function
from sfincs_simulation import Sfincs_simulation

from pathlib import Path # python 3.4+

import matplotlib.pyplot as plt

def findTrailingInt(s):
    N = len(s)
    for i in range(N):
        if s[N-1-i] not in "01234567890":
            if i>0:
                return (s[:N-i],int(s[N-i:]))
            else:
                return None
            

p = Path('./')
subdirs  = sorted([str(f) for f in p.iterdir() if f.is_dir()])

d = {}

basename = "baseCase"
base = Sfincs_simulation(basename,load_geometry=False)
print(subdirs)
maxNspecies = base.Nspecies
baseyval = base.n1Hat2
for subdir in subdirs:

    if subdir == basename:
        # skip base case
        continue
    simul = Sfincs_simulation(subdir,load_geometry=False)
    try:
        yval = simul.n1Hat2
    except KeyError:
        # skip this value
        print("Output not found for '" + subdir + "'.")
        continue
    except AttributeError:
        print("Simulation not (yet?) run for '"+ subdir + "'.")
        continue
    param, xval = findTrailingInt(subdir)
    print(param)
    
    
    if param not in d:
        d[param] = [(getattr(base,param),baseyval)]
        
    d[param].append((xval,yval))

# todo maxNspecies
Nrows = maxNspecies * 1
Ncols = len(d)
fig,axes = plt.subplots(Nrows,Ncols,sharey=True,squeeze=False)
for i,k in enumerate(d):
    for i_s in range(maxNspecies):
        for j in range(len(d[k])):
            axes[i_s,i].plot(d[k][j][0],d[k][j][1][i_s],'*b')
            axes[i_s,i].text(s="{:0.2f}".format(d[k][j][1][i_s]/baseyval[i_s]),x=d[k][j][0],y=d[k][j][1][i_s])
        axes[i_s,i].set_ylabel("n1Hat2_" + str(i_s))
        axes[i_s,i].set_xlabel(k)
        axes[i_s,i].axhline(baseyval[i_s],color='r')
        
            
plt.show()

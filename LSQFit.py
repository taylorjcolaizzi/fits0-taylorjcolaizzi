import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

xmin=1.0
xmax=20.0
npoints=12
sigma=0.2
lx=np.zeros(npoints)
ly=np.zeros(npoints)
ley=np.zeros(npoints)
pars=[0.5,1.3,0.5]

from math import log
def f(x,par):
    return par[0]+par[1]*log(x)+par[2]*log(x)*log(x)

from random import gauss
# this function is identital to np.arange()
def getX(x):  # x = array-like
    step=(xmax-xmin)/npoints
    for i in range(npoints):
        x[i]=xmin+i*step
        
def getY(x,y,ey):  # x,y,ey = array-like
    for i in range(npoints):
        y[i]=f(x[i],pars)+gauss(0,sigma)
        ey[i]=sigma

# get a random sampling of the (x,y) data points, rerun to generate different data sets for the plot below
getX(lx)
getY(lx,ly,ley)

fig, ax = plt.subplots()
ax.errorbar(lx, ly, yerr=ley)
ax.set_title("Pseudoexperiment")
fig.show()


# *** modify and add your code here ***
nexperiments = 1000  # for example


# perform many least squares fits on different pseudo experiments here
# fill histograms w/ required data

par_a = np.random.rand(1000)   # simple placeholders for making the plot example
par_b = np.random.rand(1000)   # these need to be filled using results from your fits
par_c = np.random.rand(1000)
chi2_reduced = np.random.rand(1000)

fig, axs = plt.subplots(2, 2)
plt.tight_layout()

# careful, the automated binning may not be optimal for displaying your results!
axs[0, 0].hist2d(par_a, par_b)
axs[0, 0].set_title('Parameter b vs a')

axs[0, 1].hist2d(par_a, par_c)
axs[0, 1].set_title('Parameter c vs a')

axs[1, 0].hist2d(par_b, par_c)
axs[1, 0].set_title('Parameter c vs b')

axs[1, 1].hist(chi2_reduced)
axs[1, 1].set_title('Reduce chi^2 distribution')

fig.show()

# **************************************
  

input("hit Enter to exit")

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

xmin=1.0
xmax=20.0
npoints=12*5*5
sigma=0.2
lx=np.zeros(npoints)
ly=np.zeros(npoints)
ley=np.zeros(npoints)
pars=[0.5,1.3,0.5]

from math import log
def f(x,par):
    return par[0]+par[1]*log(x)+par[2]*log(x)*log(x)

from random import gauss
# this function is identital to np.arange(i, f, n)
def getX(x):  # x = array-like
    step=(xmax-xmin)/npoints
    for i in range(npoints):
        x[i]=xmin+i*step
        
# just f(x) plus a random gaussian error
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
ax.set_xlabel("x")
ax.set_ylabel("y")
# fig.show() # no need to show plot for now.


# *** modify and add your code here ***
nexperiments = 1000*100  # for example
nPar = 3

# perform many least squares fits on different pseudo experiments here
# fill histograms w/ required data

par_a = np.random.rand(nexperiments)   # simple placeholders for making the plot example
par_b = np.random.rand(nexperiments)   # these need to be filled using results from your fits
par_c = np.random.rand(nexperiments)
chi2 = np.random.rand(nexperiments)

# for each fit, do the following:
for fit in range(nexperiments):
    getX(lx)
    getY(lx, ly, ley)
    function_points = []
    for i in range(npoints):
        function_points.append(f(lx[i], pars)) #this is the function without gaussian errors
    residuals = ly - function_points
    residuals_squared = residuals * residuals
    chi_square = sum(residuals_squared / (sigma * sigma)) # since all points have same sigma

    ax = np.array(lx, dtype = 'd')
    ax2 = lx*lx
    ay = np.array(ly, dtype = 'd')
    ae = np.array(ley, dtype = 'd')
    nPnts = len(ax)

    A = np.matrix(np.zeros((nPnts, nPar)))

    # fill matrix

    for nr in range(nPnts):
        for nc in range(nPar):
            if nc == 0: A[nr, nc] = 1
            if nc == 1: A[nr, nc] = ax[nr]
            if nc == 2: A[nr, nc] = ax2[nr]
    for i in range(nPnts):
        A[i] = A[i] / ae[i] # weight with errors
    yw = (ay/ae).reshape(nPnts, 1)
    theta = inv(np.transpose(A).dot(A)).dot(np.transpose(A)).dot(yw)

    # transform theta from a jank 1x1 numpy matrix into a sensible python 3 element list.

    theta_list = list(theta.A.flatten())
    # make list from theta as flattened array

    # now fill the placeholders.
    par_a[fit] = theta_list[0]
    par_b[fit] = theta_list[1]
    par_c[fit] = theta_list[2]
    chi2[fit] = chi_square
chi2_reduced = chi2 / (nPnts - nPar)


# to save as pdfs, just save the figures and add them to a 
# word doc later on in life

# now, plot the histograms of the fit parameters

fig, axs = plt.subplots(2, 2)
plt.tight_layout()

axs[0, 0].hist(par_a, bins = 100)
axs[0, 0].set_title('Parameter a')
axs[0, 1].hist(par_b, bins = 100)
axs[0, 1].set_title('Parameter b')
axs[1, 0].hist(par_c, bins = 100)
axs[1, 0].set_title('Parameter c')
axs[1, 1].hist(chi2, bins = 100)
axs[1, 1].set_title('Chi-Square')

plt.tight_layout()
fig.savefig('python_1d_histograms.pdf')
fig.show()

##################################

fig, axs = plt.subplots(2, 2)

# careful, the automated binning may not be optimal for displaying your results!
axs[0, 0].hist2d(par_a, par_b, bins = 100)
axs[0, 0].set_title('Parameter b vs a')

axs[0, 1].hist2d(par_a, par_c, bins = 100)
axs[0, 1].set_title('Parameter c vs a')

axs[1, 0].hist2d(par_b, par_c, bins = 100)
axs[1, 0].set_title('Parameter c vs b')

# to calculate reduced chi^2, all we need to do is to divide the chi_square by the number of degrees of freedom.
# the degrees of freedom is just the number of data points minus the number of fit parameters.
# in this case, reduces chi^2 is equal to
# chi^2 ?=/ (nexperiments - nPar)

axs[1, 1].hist(chi2_reduced, bins = 100)
axs[1, 1].set_title('Reduced Chi^2')

plt.tight_layout()
fig.savefig('python_2d_histograms.pdf')
fig.show()

# **************************************
  

input("hit Enter to exit")
print(np.mean(chi2), 'is chi-2 mean')
print(np.mean(chi2_reduced), 'is reduced chi-2 mean')
print(np.std(chi2), 'is chi-2 stdev')
print(np.std(chi2_reduced), 'is reduced chi-2 stdev')
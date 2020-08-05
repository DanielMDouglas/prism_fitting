import numpy as np
import scipy.stats as st
from scipy.linalg import block_diag
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from ROOT import *

def root_to_array(infileName, branchName, binEdges = []):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    con = np.ndarray(shape)
    err = np.ndarray(shape)
    for i in xrange(shape[0]):
        for j in xrange(shape[1]):
            for k in xrange(shape[2]):
                con[i,j,k] = TH.GetBinContent(i+1, j+1, k+1)
                err[i,j,k] = TH.GetBinError(i+1, j+1, k+1)
    infile.Close()
    con = con.squeeze()
    err = err.squeeze()

    if list(binEdges):
        oldBins = list(root_to_axes(infileName, branchName))

        for axis, theseBinEdges in enumerate(binEdges):
            con = average_by_bin_edge(con, oldBins[axis], theseBinEdges, axis = axis)
            err = average_by_bin_edge(err, oldBins[axis], theseBinEdges, axis = axis)

    return con, err


def tgraph_to_array(infileName, branchName, bins):
    infile = TFile(infileName)
    TG = infile.Get(branchName)
    con = np.ndarray(bins.shape)
    for i, xi in enumerate(bins):
        con[i] = TG.Eval(xi)
    infile.Close()
    return con


def root_to_axes(infileName, branchName, where = 'mid'):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    axes = [TH.GetXaxis(),
            TH.GetYaxis(),
            TH.GetZaxis()]
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    if where == 'mid':
        xBins = np.array([axes[0].GetBinCenter(i+1) for i in range(shape[0])])
        yBins = np.array([axes[1].GetBinCenter(i+1) for i in range(shape[1])])
        zBins = np.array([axes[2].GetBinCenter(i+1) for i in range(shape[2])])
    elif where == 'pre':
        xBins = np.array([axes[0].GetBinLowEdge(i+1) for i in range(shape[0])])
        yBins = np.array([axes[1].GetBinLowEdge(i+1) for i in range(shape[1])])
        zBins = np.array([axes[2].GetBinLowEdge(i+1) for i in range(shape[2])])
    elif where == 'post':
        xBins = np.array([axes[0].GetBinUpEdge(i+1) for i in range(shape[0])])
        yBins = np.array([axes[1].GetBinUpEdge(i+1) for i in range(shape[1])])
        zBins = np.array([axes[2].GetBinUpEdge(i+1) for i in range(shape[2])])
    else:
        print("options are 'mid', 'pre', and 'post'")
        return
    
    return (xBins, yBins, zBins)
    

def solution_norm(flux_matrix, target, lamb):
    nBinsOA = flux_matrix.shape[1]
    A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
    c = coefficients(flux_matrix, target_vect, lamb)
    return np.sqrt((np.sum(np.power(np.dot(A, c),2))))


def residual_norm(flux_matrix, error_array, target, lamb):
    c = coefficients(flux_matrix, target, lamb)
    solution_vec = np.dot(flux_matrix, c).flatten()
    residual_vec = np.array(solution_vec - target).squeeze()
    # var_vec = np.array(np.dot(np.power(error_array, 2), np.power(c, 2))).squeeze()
    var_vec = np.array(np.dot(np.power(error_array, 2), np.power(c, 2)).flatten() + np.power(0.1*target, 2)).squeeze()

    return np.sqrt(np.sum(np.power(residual_vec, 2)/var_vec))


def optimize_reg(flux_matrix, err_array, target):
    l_space = np.logspace(-10, -1, 1000)
    sol_norm = np.array([solution_norm(flux_matrix, target, l)
                         for l in l_space])
    res_norm = np.array([residual_norm(flux_matrix, err_array, target, l)
                         for l in l_space])

    dl = np.diff(l_space)
    xi = np.log(sol_norm)
    rho = np.log(res_norm)

    plt.plot(xi, rho)
    plt.show()
    
    xi_prime = np.diff(xi)/dl
    rho_prime = np.diff(rho)/dl

    xi_prime_prime = np.diff(xi_prime)/dl[:-1]
    rho_prime_prime = np.diff(rho_prime)/dl[:-1]

    curv = 2*(rho_prime[:-1]*xi_prime_prime - rho_prime_prime*xi_prime[:-1])/np.power(np.power(rho_prime[:-1], 2) + np.power(xi_prime[:-1], 2), 3./2)

    plt.plot(l_space[1:-1], curv)
    plt.semilogx()
    plt.show()

    opt_lambda = l_space[1:-1][curv==np.max(curv)][0]

    return opt_lambda


def resize_hist_1(oldHist, oldBins, newBins):
    nBinsOld = len(oldBins)
    nBinsNew = len(newBins)
    newHist = np.zeros_like(newBins)
    for i in range(nBinsNew):
        newHist[i] = oldHist[int(float(i*nBinsOld/nBinsNew))]

    newHist *= np.sum(oldHist)/np.sum(newHist)
    return newHist
        

def resize_hist_2(oldHist, oldBinsX, newBinsX, oldBinsY, newBinsY):
    nBinsOldX = len(oldBinsX)
    nBinsNewX = len(newBinsX)
    
    nBinsOldY = len(oldBinsY)
    nBinsNewY = len(newBinsY)

    newHist = np.zeros(shape = (nBinsNewX, nBinsNewY))
    for i in range(nBinsNewX):
        for j in range(nBinsNewY):
            newHist[i][j] = oldHist[int(float(i*nBinsOldX/nBinsNewX))][int(float(j*nBinsOldY/nBinsNewY))]

    return newHist


def rebin(oldHist, rebinF, axis = 0):
    oldShape = oldHist.shape
    newShape = oldShape[:axis] + (oldShape[axis]/rebinF, rebinF) + oldShape[axis+1:]
    newHist = np.sum(oldHist.reshape(newShape), axis = axis + 1)
    return newHist


def rebin_by_bin_edge(oldHist, oldBinCenters, newBinEdges, axis = 0):
    # WARNING: Only works for axis = 0 or 1 for now!
    oldShape = oldHist.shape
    newShape = oldShape[:axis] + tuple((newBinEdges.size-1,)) + oldShape[axis+1:]
    newHist = np.ndarray(newShape)
    if axis == 0:
        for i, (leftEdge, rightEdge) in enumerate(zip(newBinEdges[:-1], newBinEdges[1:])):
            newHist[i] = np.sum(oldHist[np.logical_and(leftEdge < oldBinCenters,
                                                       oldBinCenters <= rightEdge)],
                                axis = axis)
        return newHist
    elif axis == 1:
        for i in range(newHist.shape[0]):
            newHist[i] = rebin_by_bin_edge(oldHist[i], oldBinCenters, newBinEdges)
        return newHist


def average(oldHist, rebinF, **kwargs):
    return rebin(oldHist, rebinF, **kwargs)/float(rebinF)


def average_by_bin_edge(oldHist, oldBinCenters, newBinEdges, axis = 0):
    rebinned = rebin_by_bin_edge(oldHist, oldBinCenters, newBinEdges, axis = axis)
    # binWidths = np.ndarray((rebinned.shape[axis]))
    for i, (leftEdge, rightEdge) in enumerate(zip(newBinEdges[:-1], newBinEdges[1:])):
        nInside = np.sum(np.logical_and(leftEdge < oldBinCenters,
                                        oldBinCenters <= rightEdge),
                         dtype = float)
        rebinned[i] /= nInside
    return rebinned


def cut_arrays(ND, FD, FDunosc, Ebins, OAbins, Emax = 4, OAmax = None):
    # peak finding
    threshold = np.mean(FD)
    peaks = Ebins[np.pad(np.diff(FD, n=2), 1, 'constant') < -threshold]

    end = peaks[-1]
    # ax1.axvline(x = end, ls = '--')

    size = 0.15
    norm = np.sqrt(2*np.pi)*size*FD[Ebins == end]
    
    FD = np.where(Ebins > end, FD, norm*st.norm.pdf(Ebins, loc = end, scale = size))

    # cut off at 4 GeV
    FD = FD[Ebins < Emax]
    FDunosc = FDunosc[Ebins < Emax]
    if OAmax:
        ND = ND[Ebins < Emax].T[OAbins <= OAmax].T
        OAbins = OAbins[OAbins <= OAmax]
    else:
        ND = ND[Ebins < Emax]
    Ebins = Ebins[Ebins < Emax]
        
    return ND, FD, FDunosc, Ebins, OAbins


def plot_with_bands(x, y_coll, ax = plt, *args, **kwargs):
    """
    Plot the median and +/- 1 sigma values as a line and band 
    """
    quantiles = [0.16, 0.5, 0.84]
    lastAxis = len(y_coll.shape) - 1
    lower, med, upper = np.quantile(y_coll, quantiles, axis = lastAxis)
    
    band = ax.fill_between(x, lower, upper, alpha = 0.5, *args, **kwargs)
    line, = ax.plot(x, med, *args, **kwargs)
    return band, line


def float_to_sci(thisFloat, digits = 2):
    raw_float_string = str(thisFloat)
    if 'e' in raw_float_string:
        roundFactor = -int(raw_float_string.split('e')[-1]) + digits
        float_string = str(round(thisFloat, roundFactor)).replace('+', '')
        return r'$'+float_string.replace('e', r'\times 10^{')+'}$'
    else:
        if thisFloat < 1:
            shift = -5
        else:
            shift = 5
        adjFloat = thisFloat*(10**shift)
        newStr = float_to_sci(adjFloat, digits)
        expBeg = newStr.index('{')+1
        expEnd = newStr.index('}')
        exp = int(newStr[expBeg:expEnd])
        newExp = exp - shift
        return newStr[:expBeg] + str(newExp) + newStr[expEnd:]
        

from matplotlib.widgets import Slider


class Sliderlog(Slider):

    """Logarithmic slider.

    Takes in every method and function of the matplotlib's slider.

    Set slider to *val* visually so the slider still is linear but display 10**val next to the slider.

    Return 10**val to the update function (func)"""

    def set_val(self, val):

        xy = self.poly.xy
        # if self.orientation == 'vertical':
        #     xy[1] = 0, val
        #     xy[2] = 1, val
        # else:
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        # self.valtext.set_text(self.valfmt % 10**val)   # Modified to display 10**val instead of val
        # self.valtext.set_text(float_to_sci(10**val))   # Modified to display 10**val instead of val
        self.valtext.set_text("")   # Modified to display 10**val instead of val
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        # self.val = val
        # if not self.eventson:
        #     return
        # for cid, func in self.observers.items():
        #         func(10**val)
        self.val = 10**val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
                func(10**val)

import numpy as np
import scipy.stats as st
from ROOT import *

def root_to_array(infileName, branchName):
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
    return con.squeeze(), err.squeeze()


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
    else:
        print "options are 'mid' and 'pre'"
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

    print sol_norm.shape, res_norm.shape
    
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
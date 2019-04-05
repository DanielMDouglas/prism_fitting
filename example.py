import numpy as np
import matplotlib.pyplot as plt
from ROOT import *

def root_to_array(infileName, branchName):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    con = np.empty(shape)
    err = np.empty(shape)
    for i in xrange(shape[0]):
        for j in xrange(shape[1]):
            for k in xrange(shape[2]):
                con[i,j,k] = TH.GetBinContent(i+1, j+1, k+1)
                err[i,j,k] = TH.GetBinError(i+1, j+1, k+1)
    infile.Close()
    return con.squeeze(), err.squeeze()

def root_to_axes(infileName, branchName):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    mins = [TH.GetXaxis().GetXmin(),
            TH.GetYaxis().GetXmin(),
            TH.GetZaxis().GetXmin()]
    maxs = [TH.GetXaxis().GetXmax(),
            TH.GetYaxis().GetXmax(),
            TH.GetZaxis().GetXmax()]
    axes = [np.linspace(min_i, max_i, nBins_i)
            for min_i, max_i, nBins_i
            in zip(mins, maxs, shape)]
    return tuple(axes)

def coefficients(flux_matrix, target, lamb):
    nBinsOA = flux_matrix.shape[1]
    A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
    Gamma = lamb*A
    return np.dot(np.matmul((np.matmul(flux_matrix.T,flux_matrix) + np.matmul(Gamma.T, Gamma)).I, flux_matrix.T), target).T


if __name__ == '__main__':

    # Load FD
    FDoscillatedFileName = "../flux/oscillated/nominal.FHC.Oscillated.dm230.0026.sin2theta230.5.uniform.root"
    FDoscillatedBranchName = "flux"
    FDoscillatedFlux, FDoscillatedFluxStat = root_to_array(FDoscillatedFileName, FDoscillatedBranchName)

    # Load ND
    NDFileName = "../flux/nominal.FHC.uniform.root"
    NDBranchName = "numu_flux_2D"
    NDFlux, NDFluxStat = root_to_array(NDFileName, NDBranchName)
    NDMatrix = np.matrix(NDFlux)

    # Load axes (use ND, they are the same for these files,
    # otherwise, you have to rebin them so that they are the same)
    Ebins, OAbins, junk = root_to_axes(NDFileName, NDBranchName)
    print "N E bins:", Ebins.size, "from ", Ebins[0], "to ", Ebins[-1]
    print "N OA bins:", OAbins.size, "from ", OAbins[0], "to ", OAbins[-1]

    fig, (fluxAx, coeffAx) = plt.subplots(1, 2)
    fig.set_figheight(4.8)
    fig.set_figwidth(12.8)
    
    fluxAx.plot(Ebins,
                FDoscillatedFlux,
                color = 'red',
                label = r'Target')

    # set the regularization factor
    # l = 0 is no regularization and
    # coefficients reduce to a simple matrix inversion
    l = 0
    c = coefficients(NDMatrix, FDoscillatedFlux, l)
    fluxAx.plot(Ebins,
                np.dot(NDMatrix, c),
                color = 'blue',
                label = r'Fit, no regularization')
    coeffAx.plot(OAbins,
                 c,
                 color = 'blue')
    
    # Now, with regularization
    # the right amount depends entirely on the normalization of the target
    # and the normalization of the ND flux
    # for now, l = 1.e-10 will give us an idea...
    l = 1.e-10
    c = coefficients(NDMatrix, FDoscillatedFlux, l)
    fluxAx.plot(Ebins,
                np.dot(NDMatrix, c),
                color = 'green',
                label = r'Fit, $\lambda = 10^{-10}$')
    coeffAx.plot(OAbins,
                 c,
                 color = 'green')

    # plot labels and whatnot
    fluxAx.set_title(r'Flux Fit')
    fluxAx.set_xlabel(r'$E$ [GeV]')
    fluxAx.set_ylabel(r'$\Phi$')
    fluxAx.set_xlim([0, 4])
    fluxAx.legend()
    fluxAx.grid()

    coeffAx.set_title(r'Coefficients')
    coeffAx.set_xlabel(r'$D_{OA}$ [m]')
    coeffAx.grid()
    
    plt.show()

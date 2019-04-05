import pickle
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# from ROOT import TH1D, TH2D, TH1F, TH2F, TH1, TH2, TFile, TGraph
from ROOT import *
from scipy.optimize import fmin_l_bfgs_b as minimizer
import scipy.stats as st

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


def plot_fit(flux_matrix, target, unosc, lamb, Ebins, OAbins, coeffErr, fitErr, compErr, title, label, fig, ((ax1, ax2), (ax3, ax4))):
    c = np.array(coefficients(flux_matrix, target, lamb))

    # fig = plt.figure()
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    fig.set_figheight(4.8)
    fig.set_figwidth(12.8)

    # ax1 = fig.add_subplot(221)
    # ax1.axvline(x = cent_en, ls = '--', color = 'red')
    fit = np.array(np.dot(flux_matrix, c)).squeeze()
    ax1.fill_between(Ebins,
                     fit+fitErr,
                     fit-fitErr,
                     alpha = 0.5)
    p1 = ax1.plot(Ebins, fit)
    p2 = ax1.plot(Ebins, target)
    p3 = ax1.fill(np.NaN, np.NaN, color = '#1f77b4', alpha = 0.5)
    ax1.set_xlim(Ebins[0], Ebins[-1])
    # ax1.set_xlabel(r'$E$ (GeV)')
    # ax1.set_ylim(0, 4)
    ax1.set_ylabel(r'$\Phi$')
    ax1.grid()
    ax1.legend([(p3[0], p1[0]), p2[0]],
               [label, r'FD Target'],
               title = title)
    # ax1.legend(title=r'Sum of squared residuals per bin: [something else]')
    ax1.set_title("Flux fits")

    
    # def chi2(args):
    #     return np.sum(np.power(args[2]*st.norm.pdf(Ebins, loc = args[0], scale = args[1]) - target, 2)[np.logical_and(0.45 < Ebins, Ebins < 0.55)])

    # print minimizer(chi2,
    #                 [0.5, 0.1, 2.e-15],
    #                 approx_grad = True,
    #                 factr = 1.e30)
    # ax1.plot(Ebins, st.
    
    ax5 = fig.add_subplot(122)
    
    ax5.fill_between(OAbins,
                     np.array(c).squeeze() + coeffErr,
                     np.array(c).squeeze() - coeffErr,
                     alpha = 0.5)
    ax5.plot(OAbins, c)
    ax5.set_xlabel(r'$D_{OA}$ (m)')
    ax5.set_xlim(OAbins[0], OAbins[-1])
    # ylim = 3.e-4
    # ax5.set_ylim(-ylim, ylim)
    ax5.grid()
    # ax5.legend(title = label)
    ax5.set_title('Coefficients')

    # ax3 = fig.add_subplot(223)
    # ax3.axvline(x = cent_en, ls = '--', color = 'red')
    comp = (np.array(np.dot(flux_matrix, c)).flatten() - target)/unosc
    compErrUp = (fit + fitErr - target)/unosc
    compErrLow = (fit - fitErr - target)/unosc
    # ax3.fill_between(Ebins,
    #                  comp - compErr,
    #                  comp + compErr,
    #                  alpha = 0.5)
    ax3.fill_between(Ebins,
                     compErrUp,
                     compErrLow,
                     alpha = 0.5)
    ax3.axhline(y = 0, linestyle='--', color = '#ff7f0e')
    ax3.plot(Ebins, comp)
    ax3.set_xlim(0, 4)
    ylim = 0.1
    ax3.set_ylim(-ylim, ylim)
    ax3.set_xlabel(r'$E$ (GeV)')
    ax3.set_ylabel(r'$\frac{ND - FD}{FD_{unosc}}$')
    ax3.grid()

    ax2.axis('off')
    ax4.axis('off')

    fig.subplots_adjust(hspace=0)
    
    # plt.show()
    # plt.savefig('foo.png')
    # plt.clf()
    # plt.close()


def resize_hist_1(oldHist, oldBins, newBins):
    nBinsOld = len(oldBins)
    nBinsNew = len(newBins)
    newHist = np.zeros_like(newBins)
    for i in range(nBinsNew):
        newHist[i] = oldHist[int(float(i*nBinsOld/nBinsNew))]

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


if __name__ == '__main__':

    filenames = ["../flux/oscillated/nominal.FHC.Oscillated.dm230.0022.sin2theta230.4.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.0022.sin2theta230.5.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.0022.sin2theta230.6.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.0026.sin2theta230.4.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.0026.sin2theta230.5.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.0026.sin2theta230.6.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.003.sin2theta230.4.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.003.sin2theta230.5.uniform.root",
                 "../flux/oscillated/nominal.FHC.Oscillated.dm230.003.sin2theta230.6.uniform.root"]
    labels = [r'$\Delta m^2_{32} = 2.2 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.4$',
              r'$\Delta m^2_{32} = 2.2 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.5$',
              r'$\Delta m^2_{32} = 2.2 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.6$',
              r'$\Delta m^2_{32} = 2.6 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.4$',
              r'$\Delta m^2_{32} = 2.6 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.5$',
              r'$\Delta m^2_{32} = 2.6 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.6$',
              r'$\Delta m^2_{32} = 3.0 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.4$',
              r'$\Delta m^2_{32} = 3.0 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.5$',
              r'$\Delta m^2_{32} = 3.0 (10^{-3} eV^2), \sin^2 \theta_{23} = 0.6$']

    # OAmaxes = [33, 27, 20]

    OAmax = 37.5
    
    for oscFDFileName, oscLabel in zip(filenames, labels):
        print oscFDFileName
        # 400, [0,10]
        # oscFDFileName         = "../flux/oscillated/nominal.FHC.Oscillated.dm230.0025.sin2theta230.5.uniform.root"
        FDBranchName          = "flux"
        FDnominal, targetStat = root_to_array(oscFDFileName, FDBranchName)

        FDUnoscBranchName = "LBNF_numu_flux"
        FDunoscNominal, usoscStat = root_to_array(oscFDFileName, FDUnoscBranchName)
    
        # 400x390, [0,10], [-1.5,37.5]
        NDFileName        = "../flux/nominal.FHC.uniform.root" 
        NDBranchName      = "numu_flux_2D"
        NDnominal, NDstat = root_to_array(NDFileName, NDBranchName)
    
        Ebins, OAbins, zbins = root_to_axes(NDFileName, NDBranchName)

        NDnominalMatrix = np.matrix(NDnominal)
        l = 1.e-8

        figStuff = plt.subplots(2, 2, sharex='col')
        
        # for OAmax in OAmaxes:
        #     print OAmax
            
        NDcutNominalMatrix, FDcutNominal, FDcutUnoscNominal, cutEbins, cutOAbins = cut_arrays(NDnominalMatrix,
                                                                                              FDnominal,
                                                                                              FDunoscNominal,
                                                                                              Ebins,
                                                                                              OAbins,
                                                                                              OAmax = OAmax)

            
        # 40x81, [0, 20], [-0.25,40.25]
        Ntrials = 1000
        Nparams = 20

        systFileName = "../flux/syst/FluxErrors_40mOffAxis_Total_BothBeamModes_AllSpecies.root"

        NDshifts = [root_to_array(systFileName,
                                  "EffectiveFluxParameters/param_"+paramName+"/ND_nu_numu")
                    for paramName in [str(i) for i in range(Nparams)]]
        FDshifts = [root_to_array(systFileName,
                                  "EffectiveFluxParameters/param_"+paramName+"/FD_nu_numu")
                    for paramName in [str(i) for i in range(Nparams)]]

        shiftEbins, shiftOAbins, shiftZbins = root_to_axes(systFileName,
                                                           "EffectiveFluxParameters/param_0/ND_nu_numu")


        NDshiftsResized = np.array([resize_hist_2(NDshift, shiftEbins, Ebins, shiftOAbins, OAbins)
                                    for NDshift, NDshiftErr in NDshifts])
        FDshiftsResized = np.array([resize_hist_1(FDshift, shiftEbins, Ebins)
                                    for FDshift, FDshiftErr in FDshifts])

        coeffs_array = np.empty((Ntrials, cutOAbins.shape[0]))
        fits_array = np.empty((Ntrials, cutEbins.shape[0]))
        comp_array = np.empty((Ntrials, cutEbins.shape[0]))
    
        for trial in range(Ntrials):
            shiftMag = np.random.normal(size = Nparams)

            NDshift = np.dot(NDshiftsResized.T, shiftMag).T
            FDshift = np.dot(FDshiftsResized.T, shiftMag)
        
            NDshifted = NDnominal*(1 + NDshift)
            FDshifted = FDnominal*(1 + FDshift)
            FDunoscShifted = FDunoscNominal*(1 + FDshift)
            
            NDshiftedMatrix = np.matrix(NDshifted)
            NDcutShiftedMatrix, FDcutShifted, FDcutUnoscShifted, foo, bar = cut_arrays(NDshiftedMatrix,
                                                                                       FDshifted,
                                                                                       FDunoscShifted,
                                                                                       Ebins,
                                                                                       OAbins,
                                                                                       OAmax = OAmax)

    
            c = coefficients(NDcutShiftedMatrix, FDcutShifted, l).flatten()
            
            coeffs_array[trial] = c

            fit = np.dot(NDcutShiftedMatrix, c.T).flatten()
            
            fits_array[trial] = fit

            # comp = np.array(fit - FDcutShifted)/FDcutUnoscShifted

            # comp_array[trial] = comp
        
            # plot_fit(NDcutShiftedMatrix,
                #          FDcutShifted,
                #          FDcutUnoscShifted,
                #          l,
                #          cutEbins,
                #          cutOAbins,
                #          "param "+paramName,
                #          *figStuff)

        coeffs_std = np.std(coeffs_array, axis=0)
        fits_std = np.std(fits_array, axis=0)
        comp_std = np.std(comp_array, axis=0)
        comp_bounds = np.percentile(comp_array, [16., 84.], axis=0)
        
        plot_fit(NDcutNominalMatrix,
                 FDcutNominal,
                 FDcutUnoscNominal,
                 l,
                 cutEbins,
                 cutOAbins,
                 coeffs_std,
                 fits_std,
                 comp_bounds,
                 oscLabel,
                 r'ND Combined Flux $\pm 1 \sigma$',
                 # r'' + str(OAmax) + r'm $\pm 1 \sigma$',
                 *figStuff)
    
        # plt.show()
        plt.savefig(oscFDFileName.replace("../flux/oscillated/", "oscPlots/pdf/"+str(OAmax)+"m/").replace(".root", "_"+str(OAmax)+"m.pdf"))

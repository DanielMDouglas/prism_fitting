import pickle
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ROOT import TH1D, TH2D, TFile, TGraph

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
    return np.dot(np.dot((np.dot(flux_matrix.T,flux_matrix) + np.dot(Gamma.T, Gamma)).I, flux_matrix.T), target).T


def solution_norm(flux_matrix, target, lamb):
    nBinsOA = flux_matrix.shape[1]
    A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
    c = coefficients(flux_matrix, target_vect, lamb)
    return np.sqrt((np.sum(np.power(np.dot(A, c),2))))


def residual_norm(flux_matrix, target, target_norm, lamb):
    c = coefficients(flux_matrix, target, lamb)
    solution_vec = np.dot(flux_matrix, c).flatten()
    residual_vec = np.array(solution_vec - target).squeeze()

    # var_vec = np.where(target,
    #                    np.power(solution_vec, 2),
    #                    np.power(target_norm, 2))
    # var_vec = np.power(0.0001*target, 2) + np.power(0.00005*target_norm, 2)
    var_vec = np.ones_like(target)*np.power(target_norm, 2)*0.05
    
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
        

def plot_fit(flux_matrix, target, lamb, coeff_err, fit_err, OAbins, cent_en, name, ssr):
    c = coefficients(flux_matrix, target, lamb)

    # fig = plt.figure()
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_figheight(4.8)
    fig.set_figwidth(12.8)
    
    # ax1 = fig.add_subplot(221)
    fit = np.array(np.dot(flux_matrix, c)).squeeze()

    ax1.axvline(x = cent_en, ls = '--', color = '#d62728')
    ax1.fill_between(Ebins,
                     fit + fit_err,
                     fit - fit_err,
                     alpha = 0.5)
    p2 = ax1.plot(Ebins, target, color = '#ff7f0e')
    p1 = ax1.plot(Ebins, fit, color = '#1f77b4')
    p3 = ax1.fill(np.NaN,
                  np.NaN,
                  color = '#1f77b4',
                  alpha = 0.5)
    # ax1.plot(Ebins, target, label = r'Target, $\langle E \rangle$ = '+str(float(int(100*cent_en))/100) + r' GeV', color = '#ff7f0e')
    # ax1.plot(Ebins, np.dot(flux_matrix, c), label = r'ND Combined Flux', color = '#1f77b4')
    ax1.set_xlim(Ebins[0], Ebins[-1])
    ax1.set_xlabel(r'$E$ (GeV)')
    ax1.set_xlim(0, 4)
    ax1.set_ylabel(r'$\Phi$')
    ax1.grid()
    # ax1.legend(title=r'$\chi^2$ per bin: $10^{'+ ssr + r'}$')
    ax1.legend([(p3[0], p1[0]), p2[0]],
               [r'ND Combined Flux $\pm 1 \sigma$', r'Target, $\langle E \rangle$ = '+str(float(int(100*cent_en))/100) + r' GeV'])
    ax1.set_title("Flux fits")

    # ax2 = fig.add_subplot(122)
    ax2.fill_between(OAbins,
                     np.array(c).squeeze() + coeff_err,
                     np.array(c).squeeze() - coeff_err,
                     alpha = 0.5)
    ax2.plot(OAbins, c)
    ax2.set_xlabel(r'$D_{OA}$ (m)')
    ax2.set_xlim(OAbins[0], OAbins[-1])
    # ax2.set_ylim(-5.e9, 5.e9)
    ax2.grid()
    ax2.set_title('Coefficients')

    fig.subplots_adjust(hspace=0)
    
    # plt.show()
    plt.savefig(name+'.png')
    plt.clf()
    plt.close()
    

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
        
    
if __name__ == '__main__':
    gauss_center_space = np.linspace(0.3, 1.5, 50)
    E_prob_space = np.linspace(0.3, 1.5, 1000)
    inv_gauss_center_space = np.linspace(0.5, 3.5, 50)
    inv_E_prob_space = np.linspace(0.5, 3.5, 1000)
    max_OA_pos_space = []

    # pickleFileName = "e.pick"
    
    resid = []

    Ebins, OAbins, zbins = root_to_axes("../flux/OptimizedEngineeredNov2017Review/ND_numu_numode_50cm_0_40m.root",
                                        "LBNF_numu_flux")
    target_vect = st.norm.pdf(Ebins, loc = 1, scale = 0.1)

    ND, err = root_to_array("../flux/OptimizedEngineeredNov2017Review/ND_numu_numode_50cm_0_40m.root",
                            "LBNF_numu_flux")

    Ntrials = 1000
    Nparams = 20
    
    systFileName = "../flux/syst/FluxErrors_40mOffAxis_Total_BothBeamModes_AllSpecies.root"
    NDshifts = [root_to_array(systFileName,
                              "EffectiveFluxParameters/param_"+paramName+"/ND_nu_numu")
                for paramName in [str(i) for i in range(Nparams)]]
    shiftEbins, shiftOAbins, shiftZbins = root_to_axes(systFileName,
                                                       "EffectiveFluxParameters/param_0/ND_nu_numu")
    NDshiftsResized = np.array([resize_hist_2(NDshift, shiftEbins, Ebins, shiftOAbins, OAbins)
                                for NDshift, NDshiftErr in NDshifts])
    
    peak_inds = [list(slice).index(max(slice)) for slice in ND.T]
    peak_E = [Ebins[i] for i in peak_inds]
    peak_heights = [max(slice) for slice in ND.T]
    full_ND_matrix = np.matrix(ND)
    # Ecut_ND_matrix = full_ND_matrix[:100,:]
    # full_ND_matrix = np.matrix(root_to_array("../flux/OptimizedEngineeredNov2017Review/ND_numu_numode_50cm_0_40m.root",
    #                                          "LBNF_numu_flux"))
    # l = optimize_reg(Ecut_ND_matrix, target_vect)
    # l = optimize_reg(full_ND_matrix, err, target_vect)
    l = 5.e-10
    # l = 0
    # print l

    # optimal_values = {}

    # coeffs = pickle.load(open('coeffs.pick'))
    
    for j, g_cent in enumerate(gauss_center_space):
    # for j, inv_g_cent in enumerate(inv_gauss_center_space):
        # g_cent = 1./inv_g_cent
        resid.append([])
        # targ_norm = np.power(10, 8.295)*peak_heights[[abs(peak_E_i-g_cent)
        #                                for peak_E_i in peak_E].index(sorted([abs(peak_E_i-g_cent)
        #                                                                      for peak_E_i in peak_E])[0])]
        targ_norm = peak_heights[[abs(peak_E_i-g_cent)
                                  for peak_E_i in peak_E].index(sorted([abs(peak_E_i-g_cent)
                                                                        for peak_E_i in peak_E])[0])]
        
        target_vect = targ_norm*st.norm.pdf(Ebins, loc = g_cent, scale = 0.1*g_cent)

        for i, max_OA in list(enumerate(OAbins[:30:-1]))[1:]:
        # for i, max_OA in list(enumerate(OAbins[::-1]))[1:]:

            # trunc_ND_matrix = Ecut_ND_matrix[:,:-i]
            trunc_ND_matrix = full_ND_matrix[:,:-i]
            trunc_err = err[:,:-i]
            
            # opt_l = optimize_reg(trunc_ND_matrix, target_vect)

            # optimal_values[g_cent, max_OA] = opt_l

            # l = coeffs[g_cent, max_OA]

            ssr = np.power(residual_norm(trunc_ND_matrix,
                                         target_vect, targ_norm, l), 2)/trunc_ND_matrix.shape[0]
            
            # print g_cent
            # print max_OA
            # # print l
            # print ssr
            # # print trunc_ND_matrix.shape[0]
            # print
            
            # if max_OA == 33.1625 and g_cent == 0.7408163265306122:

            if 0.20 < ssr < 0.21: 
            # if max_OA == 33.1625:
                print g_cent
                print max_OA
                print ssr
                print
                # if g_cent == 0.4959183673469387 and max_OA == 30.631249999999998:

                coeffs_array = np.empty((Ntrials, OAbins[:-i].shape[0]))
                fits_array = np.empty((Ntrials, Ebins.shape[0]))
                for trial in range(Ntrials):
                    shiftMag = np.random.normal(size = Nparams)

                    NDshift = np.dot(NDshiftsResized.T, shiftMag).T
                    NDshifted = ND*(1 + NDshift)

                    NDshiftedMatrix = np.matrix(NDshifted)[:,:-i]


                    c = coefficients(NDshiftedMatrix, target_vect, l).flatten()

                    coeffs_array[trial] = c

                    fit = np.dot(NDshiftedMatrix, c.T).flatten()

                    fits_array[trial] = fit

                coeffs_std = np.std(coeffs_array, axis = 0)
                fits_std = np.std(fits_array, axis = 0)
                    
                plot_fit(trunc_ND_matrix,
                         target_vect,
                         l,
                         coeffs_std,
                         fits_std,
                         OAbins[:-i],
                         g_cent,
                         "good_g_cent_"+str(g_cent)+"_max_OA_"+str(max_OA),
                         str(float(int(100*np.log10(ssr)))/100))
                # plot_fit(trunc_ND_matrix, target_vect, opt_l, OAbins[:-i], g_cent, j, str(float(int(1))))
            if 1. < ssr < 1.05: 
            # if max_OA == 33.1625:
                print g_cent
                print max_OA
                print ssr
                print
                # if g_cent == 0.4959183673469387 and max_OA == 30.631249999999998:

                coeffs_array = np.empty((Ntrials, OAbins[:-i].shape[0]))
                fits_array = np.empty((Ntrials, Ebins.shape[0]))
                for trial in range(Ntrials):
                    shiftMag = np.random.normal(size = Nparams)

                    NDshift = np.dot(NDshiftsResized.T, shiftMag).T
                    NDshifted = ND*(1 + NDshift)

                    NDshiftedMatrix = np.matrix(NDshifted)[:,:-i]


                    c = coefficients(NDshiftedMatrix, target_vect, l).flatten()

                    coeffs_array[trial] = c

                    fit = np.dot(NDshiftedMatrix, c.T).flatten()

                    fits_array[trial] = fit

                coeffs_std = np.std(coeffs_array, axis = 0)
                fits_std = np.std(fits_array, axis = 0)
                    
                plot_fit(trunc_ND_matrix,
                         target_vect,
                         l,
                         coeffs_std,
                         fits_std,
                         OAbins[:-i],
                         g_cent,
                         "marginal_g_cent_"+str(g_cent)+"_max_OA_"+str(max_OA),
                         str(float(int(100*np.log10(ssr)))/100))

                
            max_OA_pos_space.append(max_OA)
            resid[-1].append(ssr)

    # pickle.dump(optimal_values, open(pickleFileName, 'w'))
    max_OA_pos_space = np.unique(max_OA_pos_space)
    resid = np.array(resid)

    fig, ax1 = plt.subplots()
    
    cmesh = ax1.pcolormesh(gauss_center_space,
                           max_OA_pos_space,
                           np.flip(resid, 1).T,
                           # vmin = 0.5,
                           # vmax = 10,
                           norm=LogNorm(vmin=1.e-2, vmax=10.),
                           cmap='afmhot_r')
    # cmesh = ax1.pcolormesh(inv_gauss_center_space,
    #                        max_OA_pos_space,
    #                        np.flip(resid, 1).T,
    #                        norm=LogNorm(vmin=2.e-2, vmax=2.e1),
    #                        cmap='afmhot_r')
    plt.colorbar(cmesh, pad=0.12)


    ax2 = ax1.twinx()
    oscFile = TFile("../osc_prob/Osc_flux.root")
    TG = oscFile.Get("osc_extr/POsc")
    # TG = oscFile.Get("osc_extrrecip/POsc")
    osc = np.array([TG.Eval(Ei) for Ei in E_prob_space])
    # osc = np.array([TG.Eval(iEi) for iEi in inv_E_prob_space])
    ax2.plot(E_prob_space, osc, color='r')
    # ax2.plot(inv_E_prob_space, osc, color='r')
    ax2.tick_params('y', colors='r')
    ax2.yaxis.label.set_color('red')
    ax2.set_ylabel(r'$\nu_e$ Appearance Probability')

    cs = ax1.contour(gauss_center_space,
                     max_OA_pos_space,
                     np.flip(resid, 1).T,
                     # [np.power(10, -2.5), np.power(10, -1.5), np.power(10, -0.5)], 
                     # [0.5, 1., 3.],
                     [0.2, 1.],
                     colors = 'black')

    # cs = plt.contour(max_OA_pos_space, inv_gauss_center_space, np.flip(resid, 1), [5.e-3, 1.e-2, 1.e-1], colors = 'black')
    # ax = plt.gca()
    ax1.clabel(cs, cs.levels, inline=True)

    
    ax1.set_xlabel(r'Gaussian Peak Position (GeV)')
    ax1.set_ylabel(r'Max Off-axis Position (m)')
    # ax1.set_xlabel(r'1/Gaussian Peak Position (GeV$^{-1}$)')
    ax1.set_title(r'Sum of Squared Residuals per Bin')
    ax1.set_xlim(0.3, 1.5)
    # ax1.set_xlim(0.5, 3.5)
    plt.tight_layout()
    plt.savefig('energy.png')
    # plt.savefig('inv_energy.png')
    # plt.show()

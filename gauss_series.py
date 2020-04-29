from flux_fitter import *
from plots import *
import matplotlib as mpl

if __name__ == '__main__':

    # fine spacing for plotting the target s
    fine_E_space = np.linspace(0, 1.5, 1000)

    # these are good colors
    # colors = ['#20473D', '#92AC4C', '#6D005F', '#ED8227', '#CEDD42']
    colors = [DUNEblue, DUNElightOrange, DUNEgreen, DUNEyellow, DUNEdarkOrange]
    gauss_peak_locs = [0.3, 0.4, 0.5, 0.75, 1.0]

    maxOA = 33.
    histfig = plt.figure()
    histax = histfig.gca()
    histbins = np.linspace(0, 8e-8, 30)
    
    
    fitfig = plt.figure()
    fitax = fitfig.gca()
    
    fitter = flux_fitter("nu", "numu", "numu", "numu")

    # load systs.  This makes fitter.ND_univs available
    fitter.load_ppfx_systs()
    # set maximum off-axis distance to use in the fit
    fitter.set_maxOA(maxOA)
    # only fit inside of the plot window.  Sneaky...
    fitter.set_fit_region([0, 1.5])
    
    # not pretty, but get the normalization of the target from the height of the nearest OA flux
    peak_inds = [list(slice).index(max(slice)) for slice in fitter.ND.T]
    peak_E = [Ebins[i] for i in peak_inds]
    peak_heights = [max(slice) for slice in fitter.ND.T]
    for color, loc in zip(colors, gauss_peak_locs):
        target_norm = np.sqrt(2*np.pi*(0.1*loc)**2)*peak_heights[[abs(peak_E_i-loc) for peak_E_i in peak_E].index(sorted([abs(peak_E_i-loc) for peak_E_i in peak_E])[0])]

        fitter.target = target_norm*st.norm.pdf(Ebins, loc = loc, scale = 0.1*loc)
        fine_target = target_norm*st.norm.pdf(fine_E_space, loc = loc, scale = 0.1*loc)

        # play with this number!
        reg = 1.e-9
        fitter.calc_coeffs(reg, reg)

        fitax.plot(fine_E_space, fine_target, ls = '--', color = color)
        plot_with_bands(Ebins,
                        np.dot(fitter.ND_ppfx_univs, fitter.c).T,
                        ax = fitax,
                        color = color)
        resids = np.sqrt(np.sum(np.power(np.matmul(fitter.P,
                                                   (np.dot(fitter.ND_ppfx_univs,
                                                           fitter.c) - fitter.target).T),
                                         2),
                                axis = 0))
        histax.hist(resids,
                    color = color,
                    # histtype = 'step',
                    weights = 1./len(resids)*np.ones_like(resids),
                    alpha = 0.5,
                    bins = histbins)
        histax.axvline(x = fitter.residual_norm(),
                       color = color)

    fitax.text(0.75, 4.e-8, r'Fluxes up to '+str(maxOA)+'m')
        
    fitax.set_xlim(0, 1.5)
    fitax.set_ylim(-0.5e-8, 4.5e-8)
    
    fitax.set_xlabel(r'$E_\nu$ (GeV)')
    fitax.set_ylabel(r'$\Phi_\nu$ (A.U.)')
    fitfig.tight_layout()
    
    histax.set_xlim(histbins[0], histbins[-1])
    histax.set_xlabel(r'$|P(\Phi_{ND} \vec{c} - \vec{\Phi}_{FD})|$')
    histfig.tight_layout()
    plt.show()

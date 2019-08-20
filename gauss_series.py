from flux_fitter import *
import matplotlib as mpl

if __name__ == '__main__':

    # set matplotlib for making pretty plots
    mpl.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
    mpl.rc('text', usetex = True)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    # fine spacing for plotting the target s
    fine_E_space = np.linspace(0, 1.5, 1000)

    # these are good colors
    colors = ['#20473D', '#92AC4C', '#6D005F', '#ED8227', '#CEDD42']
    gauss_peak_locs = [0.3, 0.4, 0.5, 0.75, 1.0]

    maxOA = 30
    
    for color, loc in zip(colors, gauss_peak_locs):
        fitter = flux_fitter("nu", "numu", "numu", "numu")

        # not pretty, but get the normalization of the target from the height of the nearest OA flux
        peak_inds = [list(slice).index(max(slice)) for slice in fitter.ND.T]
        peak_E = [Ebins[i] for i in peak_inds]
        peak_heights = [max(slice) for slice in fitter.ND.T]
        target_norm = np.sqrt(2*np.pi*(0.1*loc)**2)*peak_heights[[abs(peak_E_i-loc) for peak_E_i in peak_E].index(sorted([abs(peak_E_i-loc) for peak_E_i in peak_E])[0])]

        fitter.target = target_norm*st.norm.pdf(Ebins, loc = loc, scale = 0.1*loc)
        fine_target = target_norm*st.norm.pdf(fine_E_space, loc = loc, scale = 0.1*loc)

        # load systs.  This makes fitter.ND_univs available
        fitter.load_systs()
        # set maximum off-axis distance to use in the fit
        fitter.set_maxOA(maxOA)
        # only fit inside of the plot window.  Sneaky...
        fitter.set_fit_region([0, 1.5])
        # play with this number!
        fitter.calc_coeffs(8.e-10)

        plt.plot(fine_E_space, fine_target, ls = '--', color = color)
        plot_with_bands(Ebins, np.dot(fitter.ND_univs, fitter.c).T, color = color)

    plt.text(0.75, 4.e-8, r'Fluxes up to '+str(maxOA)+'m')
        
    plt.xlim(0, 1.5)
    plt.ylim(-0.5e-8, 4.5e-8)
    
    plt.xlabel(r'$E_\nu$ (GeV)')
    plt.ylabel(r'$\Phi_\nu$ (A.U.)')
    plt.show()

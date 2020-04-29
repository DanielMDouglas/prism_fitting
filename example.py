from flux_fitter import *
from plots import *

if __name__ == '__main__':

    # make a new fitter object
    # use FHC
    # use numu -> numu at FD
    # use numu at ND
    fitter = flux_fitter("nu", "numu", "numu", "numu", useHC = True, rebin = 5)
    fitter.use_currents([])
    fitter.set_maxOA(33.)
    
    # use the following oscillation parameters
    dcp = -np.pi/2
    s23 = 0.53
    dm32 = 2.46e-3
    # a nicely formatted string for plotting
    dm32Pretty = r'$2.4 \times 10^{-3}$'

    # get the numu -> numu probability for these parameters
    osc_hyp = oscProb("numu", "numu", s23 = s23, dm32 = dm32, dcp = dcp)
    # load them into the fitter
    fitter.set_oscHypothesis(osc_hyp.load(fitter.Ebins))

    # we want to fit between the first and fourth highest energy peaks
    fitter.set_fit_region(energies = [0.4, 3.865])
    # with a 5% effect on the low-E side and no weight on the high-E side
    fitter.set_OOR([0.8, 0])
    # find the coefficients with a regularization factor of 8.e-9
    reg = 5.5e-9
    fitter.calc_coeffs(reg, reg, fluxTimesE = False)
    
    # a nicely formatted title 
    title = r'$\sin^2 \theta_{23} = $'+str(s23)+r'$, \Delta m^2_{32} = $'+dm32Pretty
    # plot the fit and the (fit - target)/target_unosc
    fp = FD_flux_plot(fitter, label = r'$\nu_\mu \rightarrow \nu_\mu$')
    ND_flux_slice_plot(fitter, slices = [0, 10, 20, 30, 40, 50])
    fitPlot = fit_and_ratio_plot()
    coeffPlot = coeff_plot(HC = True)

    fitPlot.add_target(fitter, label = r'FD $\nu_\mu \rightarrow \nu_\mu$', Ebounds = False)
    fitPlot.add(fitter, label = r'Off-axis Only')
    coeffPlot.add(fitter, label = r'Off-axis Only')
    
    fitter.use_currents([280])
    fitter.Ebounds[-1] = 10.
    fitter.calc_coeffs(reg, reg, fluxTimesE = False)
    fitPlot.add(fitter, label = r'280 kA')
    fp.add_fit(fitter)
    coeffPlot.add(fitter, label = r'280 kA')
    
    L_curv_curv_plot = L_curve_curvature_plot(fitter,
                                              regRange = np.logspace(-10, -6, 1000),
                                              fluxTimesE = True)
    opt_l = L_curv_curv_plot.opt_l[0]

    L_curve_plot(fitter,
                 regRange = np.logspace(-10, -6, 1000),
                 # ND = fitter.ND_full/fitter.FD_unoscillated.reshape(-1, 1),
                 # target = fitter.target/fitter.FD_unoscillated,
                 highlight = [opt_l])
    fitter.load_ppfx_systs()
    mike_plot(fitter, varKey = 18)
    plt.show()

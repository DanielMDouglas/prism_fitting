from flux_fitter import *
from plots import *

if __name__ == '__main__':

    # use the following oscillation parameters
    dcp = -np.pi/2
    s23 = 0.53
    dm32 = 2.46e-3
    # get the numu -> numu probability for these parameters
    osc_hyp = oscProb("numu", "nue", s23 = s23, dm32 = dm32, dcp = dcp)

    # make a new fitter object
    # use FHC
    # use numu -> numu at FD
    # use numu at ND
    fitter = flux_fitter("nu",
                         "numu", "nue",
                         "numu",
                         oscParam = osc_hyp,
                         useHC = False,
                         Erebin = 5)
    fitter.set_maxOA(33.)

    eventRate = fitter.ND_rate
    bw = fitter.EbinWidths
    print ("@4m:", sum(sum(bw*eventRate[:,i] for i in range(4, 12))))
    print ("@8m:", sum(sum(bw*eventRate[:,i] for i in range(12, 20))))
    print ("@12m:", sum(sum(bw*eventRate[:,i] for i in range(20, 28))))
    print ("@16m:", sum(sum(bw*eventRate[:,i] for i in range(28, 36))))
    print ("@20m:", sum(sum(bw*eventRate[:,i] for i in range(36, 44))))
    print ("@24m:", sum(sum(bw*eventRate[:,i] for i in range(44, 52))))
    print ("@28m:", sum(sum(bw*eventRate[:,i] for i in range(52, 60))))
    print ("@30.5m:", sum(sum(bw*eventRate[:,i] for i in range(57, 65))))
    
    # we want to fit between the first and fourth highest energy peaks
    fitter.set_fit_region(energies = [0.4, 3.865])
    # with a 5% effect on the low-E side and no weight on the high-E side
    fitter.set_OOR([0.05, 0])
    # find the coefficients with a regularization factor of 8.e-9
    reg = 5.e-9
    # reg = 0
    fitter.calc_coeffs(reg, reg, fluxTimesE = False)
    
    # a nicely formatted title 
    title = r'$\sin^2 \theta_{23} = $'+str(s23)+r'$, \Delta m^2_{32} = $'+float_to_sci(dm32)

    # plot the fit and the target
    fp = FD_flux_plot(fitter,
                      label = r'FD $\nu_\mu \rightarrow \nu_e$',
                      style = 'plot',
                      xlim = (0, 6))

    rp = FD_rate_plot(fitter, style = 'errorbandstep', label = "10 years")
    # rp.add_fit(fitter)
    
    # plot some ND fluxes at various off-axis positions
    # These are specified by the index (50cm windows from 0 to 33m) 
    ND_flux_slice_plot(fitter,
                       slices = [0, 10, 20, 30, 40, 50, 60],
                       title = r'ND $\nu_e$ Flux').save('nue_flux.png', dpi = 300)
    ND_ER_slice_plot(fitter,
                     style = 'errorbandstep',
                     slices = [0, 10, 20, 30, 40, 50, 60],
                     title = r'ND $\nu_e$ Event Rate').save('nue_rate.png', dpi = 300)
    
    coeffPlot = coeff_plot(HC = False)
    # add the fitter's coefficients
    # coeffPlot.add(fitter,
    #               label = r'Off-axis only')

    # # create a new fit and ratio plot
    fitPlot = fit_and_ratio_plot(xlim = (0, 6))
    # add the fitter's target to this plot
    fitPlot.add_target(fitter,
                       label = r'FD $\nu_\mu \rightarrow \nu_e$',
                       Ebounds = False)
    # add the fit to this plot
    # fitPlot.add(fitter,
    #             label = r'Off-axis only')
    
    fitter.use_currents([280])
    fitter.set_fit_region(energies = [0.4, 10])
    fitter.calc_coeffs(reg, reg, fluxTimesE = False)
    # fitPlot.add(fitter,
    #             label = r'Off-axis + horn current')

    # coeffPlot.add(fitter,
    #               label = r'Off-axis + horn current')
    coeffPlot.add(fitter,
                  label = None)
    fitPlot.add(fitter,
                label = r'ND Flux Match')

    fp.save("nue_app_spect.png", dpi = 300)
    
    fp.add_fit(fitter, label = r'ND Flux Match')
    rp.add_fit(fitter)
    
    # # fitPlot.save("hc_fit.pdf")
    # # coeffPlot.save("hc_coeffs.pdf")
    
    # # plt.show()
    
    # # # create a new coefficient plot
    # # coeffPlot = coeff_plot(HC = True)
    # # # add the fitter's coefficients
    # # coeffPlot.add(fitter,
    # #               label = r'Off-axis Only')

    # # make a new fitter that makes use of an alternative horn current
    # HCfitter = flux_fitter("nu",
    #                        "numu", "numu",
    #                        "numu",
    #                        oscParam = osc_hyp,
    #                        useHC = True,
    #                        Erebin = 5)
    # # configure it similarly to the first fitter
    # HCfitter.set_maxOA(33.)
    # HCfitter.use_currents([280])
    # HCfitter.set_fit_region(energies = [0.4, 10])
    # HCfitter.set_OOR([0, 0])
    # HCfitter.calc_coeffs(reg, reg,
    #                      fluxTimesE = False)

    # # add this fitter to the previous plots
    # # fitPlot.add(HCfitter,
    # #             label = r'280 kA')
    # fp.add_fit(HCfitter, label = r'ND Flux Match')
    # coeffPlot.add(HCfitter,
    #               label = r'280 kA')

    # fp.save("FluxMatch.png", dpi = 300)
    # fitPlot.save("FluxMatch_ratio.png", dpi = 300)
    # coeffPlot.save("coeffs.png", dpi = 300)
    # fp.save("FluxMatch_noreg.png", dpi = 300)
    # fitPlot.save("FluxMatch_ratio_noreg.png", dpi = 300)
    # coeffPlot.save("coeffs_noreg.png", dpi = 300)
    fitPlot.save("nue_app_fluxmatch_ratio.png", dpi = 300)
    coeffPlot.save("nue_app_coeffs.png", dpi = 300)
    
    # # to optimize the regularization, we will scan a large range of values
    # regRange = np.logspace(-10, -6, 1000)
    # # this plot will scan over these parameters and find the place
    # # of maximum curvature along the "L-curve"
    # # see (http://people.compute.dtu.dk/pcha/DIP/chap5.pdf)
    # L_curv_curv_plot = L_curve_curvature_plot(HCfitter,
    #                                           regRange = regRange)
    # # the optimal value is saved as attribute, opt_l
    # opt_l = L_curv_curv_plot.opt_l

    # # this plot will show the "L-curve" itself
    # L_curve_plot(HCfitter,
    #              regRange = regRange,
    #              highlight = [opt_l])

    # # load some hadron production throws
    # HCfitter.load_ppfx_systs(nUniv = 20)
    # # produce a "mike" plot (we need a better name)
    # mike_plot(HCfitter, varKey = 18)

    # show all of the plots we've made
    plt.show()

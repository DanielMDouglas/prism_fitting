from flux_fitter import *
from plots import *
import scipy.stats as st

newOAbins = OAbinEdges[OAbinEdges<=33.25]
    
rebinArgs = dict(Erebin = 20,
                 OArebin = newOAbins)

# fitLimLo = 0.7
fitLimLo = 0
# fitLimHi = 1
fitLimHi = 5

# reg = 8e-11
reg = 1.e-10

calc_coeffs_kwargs = {'fluxTimesE': False}

nueFitter = flux_fitter("nu",
                        "numu", "nue",
                        "nue",
                        useHC = True,
                        oscParam = oscProb("numu", "nue"),
                        **rebinArgs)

nueFitter.use_currents([280])
nueFitter.set_fit_region(energies = [fitLimLo, fitLimHi])
nueFitter.set_OOR([0, 0])

numuFitter = flux_fitter("nu",
                         "numu", "numu",
                         "numu",
                         useHC = True,
                         **rebinArgs)

numuFitter.use_currents([280])
numuFitter.set_fit_region(energies = [fitLimLo, fitLimHi])
numuFitter.set_OOR([0, 0])

xs_true = xSec["nue"]["CCInc"].load([numuFitter.EbinEdges])/xSec["numu"]["CCInc"].load([numuFitter.EbinEdges])

class xs_ratio_plot (plot):
    def __init__(self):
        self.fig = plt.figure()

        gs = GridSpec(2, 1,
                      figure = self.fig,
                      height_ratios = [0.7, 0.3],
                      hspace = 0)
        self.axUp = self.fig.add_subplot(gs[0, :])
        self.axLo = self.fig.add_subplot(gs[1, :])
        
        # ax = fig.gca()
            
        self.axUp.set_title(r'Cross Section Ratio')
            
        self.axUp.axvline(x = numuFitter.Ebounds[0],
                          ls = '--',
                          color = 'red')
        self.axUp.axvline(x = numuFitter.Ebounds[1],
                          ls = '--',
                          color = 'red')

        self.axLo.axvline(x = numuFitter.Ebounds[0],
                          ls = '--',
                          color = 'red')
        self.axLo.axvline(x = numuFitter.Ebounds[1],
                          ls = '--',
                          color = 'red')
        
        self.axLo.set_xlabel(r'$E_\nu$ [GeV]')
        self.axUp.set_ylabel(r'$\sigma(\nu_e)/\sigma(\nu_\mu)$')
        self.axLo.set_ylabel(r'Measured / True')
        
        self.axUp.semilogy()
        self.axUp.set_ylim(5.e-1, 1.e1)
        # self.axLo.set_ylim(0.4, 1.6)
        self.axLo.set_ylim(0.89, 1.11)

        # self.axUp.set_xlim(0, 1)
        # self.axLo.set_xlim(0, 1)
        
        self.axUp.step(numuFitter.Ebins,
                       xs_true,
                       where = 'mid',
                       # label = 'Genie (Input) Ratio',
                       label = 'Genie Predicted Ratio',
        )

        self.axUp.legend()
        
        plt.tight_layout()

        
    def add(self, xs_rat, xs_rat_statErr,
            label = 'ND Ratio Measurement'):
        lines = self.axUp.step(numuFitter.Ebins,
                               xs_rat,
                               where = 'mid',
                               label = label)

        color = lines[0].get_color()
        
        self.axUp.fill_between(numuFitter.Ebins,
                               xs_rat - xs_rat_statErr,
                               xs_rat + xs_rat_statErr,
                               step = 'mid',
                               alpha = 0.5,
                               color = color)

        self.axLo.step(numuFitter.Ebins,
                       xs_rat/xs_true,
                       where = 'mid',
                       color = color)

        self.axLo.fill_between(numuFitter.Ebins,
                               (xs_rat - xs_rat_statErr)/xs_true,
                               (xs_rat + xs_rat_statErr)/xs_true,
                               step = 'mid',
                               alpha = 0.5,
                               color = color)

        self.axUp.legend()
        
        plt.tight_layout()


class xs_ratio_plot_simple (plot):
    def __init__(self):
        self.fig = plt.figure()

        self.ax = self.fig.gca()

        self.ax.set_title(r'Cross Section Ratio')
                    
        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\sigma(\nu_e)/\sigma(\nu_\mu)$')
        
        # self.ax.semilogy()
        self.ax.axhline(y = 1, ls = '--', color = 'gray')
        # self.ax.set_ylim(5.e-1, 1.e1)
        # self.axLo.set_ylim(0.4, 1.6)

        # self.axUp.set_xlim(0, 1)
        # self.axLo.set_xlim(0, 1)
        
        self.ax.step(numuFitter.Ebins,
                     xs_true,
                     where = 'mid',
                     # label = 'Genie (Input) Ratio',
                     label = 'Genie Predicted Ratio',
        )

        self.ax.set_xlim(0.2, 5)
        self.ax.set_ylim(0.8, 1.5)

        self.ax.legend()
        
        plt.tight_layout()

def fit_injected_target(target,
                        injTarOutFile = None,
                        secTarOutFile = None,
                        xsRatioOutFile = None,
                        secFitOutFile = None,
                        secRatesOutFile = None):

    # inj_target_fit_ratio_plot = fit_and_ratio_plot(style = 'step',
    #                                                full_flux = False,
    #                                                title = 'Secondary Target') # flux + ratio plot for the injected target
    # inj_target_rate_ratio_plot = fit_and_ratio_rate_plot(style = 'errorbandstep',
    #                                                      full_flux = False,
    #                                                      title = 'Secondary Target') # rate + ratio plot for the injected target
    inj_target_fit_plot = FD_flux_plot(style = 'step',
                                       full_flux = False,
                                       ltitle = 'Target Definition') # flux + ratio plot for the injected target
    sec_tar_fit_plot = FD_flux_plot(style = 'step',
                                    full_flux = False,
                                    ltitle = 'Secondary Target') # flux + ratio plot for the injected target
    # inj_target_rate_plot = FD_rate_plot(style = 'errorbandstep',
    #                                     full_flux = False,
    #                                     title = 'Secondary Target') # rate + ratio plot for the injected target
    
    nueFitter.target = target
    nueFitter.FD_rate = target
    nueFitter.FD_rate_statErr = np.sqrt(nueFitter.FD_rate)
 
    # inj_target_fit_ratio_plot.add_target(nueFitter,
    #                                      label = "Injected Target",
    #                                      full_flux = False) # add the injected target flux
    # inj_target_rate_ratio_plot.add_target(nueFitter,
    #                                       label = "Injected Target",
    #                                       full_flux = False,
    #                                       scatter = True) # add the injected target rate
    inj_target_fit_plot.add(nueFitter,
                            label = "Injected Target",
                            color = 'black',
                            Ebounds = True)
    sec_tar_fit_plot.add(nueFitter,
                         label = "Injected Target",
                         color = 'black',
                         Ebounds = True)

    if injTarOutFile:
        inj_target_fit_plot.save(injTarOutFile, dpi = 300)

    # inj_target_rate_plot.add(nueFitter,
    #                          label = "Injected Target",
    #                          color = 'black',
    #                          Ebounds = True)

    print ('fitting to injected target')
    nueFitter.calc_coeffs(reg, reg, **calc_coeffs_kwargs)

    # inj_target_fit_ratio_plot.add(nueFitter,
    #                               label = r'ND '+LaTeXflavor['nue']) # add the fitted ND nue flux + ratio
    # inj_target_rate_ratio_plot.add(nueFitter,
    #                                label = r'ND '+LaTeXflavor['nue']) # add the fitted ND nue rate + ratio

    sec_tar_fit_plot.add_fit(nueFitter,
                             label = r'ND '+LaTeXflavor['nue'])
    # inj_target_rate_plot.add_fit(nueFitter,
    #                              label = r'ND '+LaTeXflavor['nue'])

    if secTarOutFile:
        sec_tar_fit_plot.save(secTarOutFile, dpi = 300)
    
    ############################
    # now, fit it with a ND numu!
    ############################
    
    # fitPlot = fit_and_ratio_plot(style = 'step', full_flux = False)
    ratePlot = fit_and_ratio_rate_plot(style = 'errorbandstep', full_flux = False)
    
    numuFitter.target = nueFitter.fluxPred
    numuFitter.FD_rate = nueFitter.ratePred
    numuFitter.FD_rate_statErr = nueFitter.ratePred_statErr

    print ('fitting to secondary target')
    numuFitter.calc_coeffs(reg, reg, **calc_coeffs_kwargs)

    # fitPlot.add_target(numuFitter,
    #                    label = r'ND '+LaTeXflavor['nue']+r' target',
    #                    full_flux = False)
    ratePlot.add_target(numuFitter,
                        label = r'ND '+LaTeXflavor['nue']+r' target',
                        full_flux = False,
                        scatter = True)
    
    # fitPlot.add(numuFitter, label = r'ND '+LaTeXflavor['numu'])
    ratePlot.add(numuFitter, label = r'ND '+LaTeXflavor['numu'])

    # rp = FD_rate_plot(numuFitter,
    #                   label = r'FD $\nu_e$',
    #                   title = '',
    #                   style = 'errorbandstep',
    #                   figSize = (8.5, 5),
    #                   Ebounds = True,
    #                   scatter = True)
    # rp.add_fit(numuFitter,
    #            label = r'ND $\nu_e$')

    fp = FD_flux_plot(numuFitter,
                      label = r'ND $\nu_e$',
                      title = 'ND Flux Match',
                      style = 'step',
                      figSize = (8.5, 5),
                      Ebounds = True,
                      color = 'black')
    fp.add_fit(numuFitter,
               label = r'ND $\nu_\mu$',
               color = DUNEblue)

    # fitPlot.save('NDnumu_fittotarget_flux_withratio.png', dpi = 300)
    if secRatesOutFile:
        # ratePlot.save('NDnumu_fittotarget_rate_withratio.png', dpi = 300)
        ratePlot.save(secRatesOutFile, dpi = 300)
    # fp.save('NDnumu_fittotarget_flux_noratio.png', dpi = 300)
    if secFitOutFile:
        fp.save(secFitOutFile, dpi = 300)
    # rp.save('NDnumu_fittotarget_rate_noratio.png', dpi = 300)

    xs_rat = numuFitter.FD_rate/numuFitter.ratePred
    xs_rat_statErr = xs_rat*np.sqrt(np.power(numuFitter.FD_rate_statErr/numuFitter.FD_rate, 2) +
                                    np.power(numuFitter.ratePred_statErr/numuFitter.ratePred, 2))

    plt.figure()
    plt.title(r'Cross Section Ratio')
    plt.axhline(y = 1,
                ls = '--',
                color = 'gray')
    plt.fill_between(numuFitter.Ebins,
                     xs_rat - xs_rat_statErr,
                     xs_rat + xs_rat_statErr,
                     step = 'mid',
                     alpha = 0.5)
    plt.step(numuFitter.Ebins,
             xs_rat,
             where = 'mid',
             label = 'ND Ratio Measurement')
    plt.step(numuFitter.Ebins,
             xSec["nue"]["CCInc"].load([numuFitter.EbinEdges])/xSec["numu"]["CCInc"].load([numuFitter.EbinEdges]),
             where = 'mid',
             label = 'Genie (Input) Ratio')
    plt.axvline(x = numuFitter.Ebounds[0],
                ls = '--',
                color = 'red')
    plt.axvline(x = numuFitter.Ebounds[1],
                ls = '--',
                color = 'red')
    plt.legend()
    
    # plt.ylim(0, 10)
    plt.ylim(0, 2)
    plt.xlabel(r'$E_\nu$ [GeV]')
    plt.ylabel(r'$\sigma(\nu_e)/\sigma(\nu_\mu)$')

    plt.tight_layout()

    if xsRatioOutFile:
        plt.savefig(xsRatioOutFile, dpi = 300)

    ######################
    # now, apply the ratio to an FD nue prediction
    ######################

    # FDfitPlot = fit_and_ratio_plot(style = 'step', full_flux = False)
    # FDratePlot = fit_and_ratio_rate_plot(style = 'errorbandstep', full_flux = False, title = r'FD Rate Match')
    
    # FDfitter = flux_fitter("nu",
    #                        "numu", "nue",
    #                        "numu",
    #                        useHC = True,
    #                        **rebinArgs)

    # FDfitter.use_currents([280])
    # FDfitter.set_fit_region(energies = [fitLimLo, fitLimHi])
    # FDfitter.set_OOR([0, 0])

    # FDfitter.target = nueFitter.target
    # FDfitter.FD_rate = nueFitter.FD_rate
    # FDfitter.FD_rate_statErr = nueFitter.FD_rate_statErr

    # FDfitter.calc_coeffs(reg, reg, **calc_coeffs_kwargs)

    # FDfitPlot.add_target(FDfitter,
    #                      label = r'FD '+LaTeXflavor['nue'],
    #                      full_flux = False)
    # FDratePlot.add_target(FDfitter,
    #                       label = r'FD '+LaTeXflavor['nue'],
    #                       full_flux = False,
    #                       scatter = True)
    
    # FDfitPlot.add(FDfitter,
    #               label = r'ND '+LaTeXflavor['numu'])
    # FDratePlot.add(FDfitter,
    #                label = r'ND '+LaTeXflavor['numu'])

    # FDfitPlot.save('NDnumu_FDnue_flux_nocorr.png', dpi = 300)
    # FDratePlot.save('NDnumu_FDnue_rate_nocorr.png', dpi = 300)

    # FDfitter.ratePred *= xs_rat
    # FDfitter.ratePred_statErr = FDfitter.ratePred*np.sqrt(np.power(FDfitter.ratePred_statErr/FDfitter.ratePred, 2) +
    #                                                       np.power(xs_rat_statErr/xs_rat, 2))
    # FDratePlot.add(FDfitter,
    #                label = r'ND '+LaTeXflavor['numu']+', Corrected',
    #                color = DUNEdarkOrange)

    # FDfitPlot.save('NDnumu_FDnue_flux_corr.png', dpi = 300)
    # FDratePlot.save('NDnumu_FDnue_rate_corr.png', dpi = 300)
    
    # plt.show()

    return (xs_rat, xs_rat_statErr)

        
def fit_injected_target_noplots(target):

    nueFitter.target = target
    nueFitter.FD_rate = target
    nueFitter.FD_rate_statErr = np.sqrt(nueFitter.FD_rate)
    
    nueFitter.calc_coeffs(reg, reg, **calc_coeffs_kwargs)

    ############################
    # now, fit it with a ND numu!
    ############################
    
    numuFitter.target = nueFitter.fluxPred
    numuFitter.FD_rate = nueFitter.ratePred
    numuFitter.FD_rate_statErr = nueFitter.ratePred_statErr

    numuFitter.calc_coeffs(reg, reg, **calc_coeffs_kwargs)

    xs_rat = numuFitter.FD_rate/numuFitter.ratePred
    xs_rat_statErr = xs_rat*np.sqrt(np.power(numuFitter.FD_rate_statErr/numuFitter.FD_rate, 2) +
                                    np.power(numuFitter.ratePred_statErr/numuFitter.ratePred, 2))

    return (xs_rat, xs_rat_statErr)


if __name__ == '__main__':

    xrps = xs_ratio_plot_simple()
    
    plt.show()
    
    # good fit locations
    # 0.4 0.4
    # 0.7 0.07
    # 0.7 0.28 (wide)
    # 0.8 0.48
    # 1.0 0.1
    # 1.2 0.12
    #
    #
    #
    

    
    # centers = [0.5, 1, 2, 3, 4, 5]
    # widths = [0.1, 0.2, 0.4, 1, 1, 1]
    centers = np.array([])
    widths = np.array([])

    # nom_center_space = np.linspace(0.1, 2, 50)
    nom_center_space = np.linspace(0.4, 1.2, 9)
    # width_ratios = [0.1, 0.3, 0.7, 1, 2]
    width_ratios = np.linspace(0.1, 1, 10)
    # widths = np.linspace(0.01, 2, 50)

    # nom_center_space = np.linspace(0, 1.5, 51)
    # widths = np.linspace(0, 2, 51)

    # nom_center_space = np.linspace(0, 2, 41)
    # widths = np.linspace(0, 3, 51)[1:]
    # widths = np.linspace(0, 3, 61)[1:]

    # nom_center_space = np.array([2.5])
    # widths = np.array([1])
    
    for wr in width_ratios:
        centers = np.concatenate((centers, nom_center_space)) 
        widths = np.concatenate((widths, wr*nom_center_space))

    # for w in widths:
    #     centers = np.concatenate((centers, nom_center_space)) 
    #     widths = np.concatenate((widths, w*np.ones_like(nom_center_space)))
    
    # centers = [0.75, 1.2, 0.2, 0.9]
    # widths = [0.375, 1.6, 1.0, 0.375]
    
    # centers = [0.759]
    # widths = [0.3755]
    
    fig = None

    xs = []
    xsErr = []
    
    minErrs = 1.e10*np.ones_like(nueFitter.Ebins)
    minErrs_val = np.zeros_like(nueFitter.Ebins)
    bestLoc = np.zeros_like(nueFitter.Ebins)
    bestWid = np.zeros_like(nueFitter.Ebins)
    
    for center, width in zip(centers, widths):
        # iterate over targets
        print (round(center, 2), round(width, 2))

        target = 1.e-16*st.norm.pdf(nueFitter.Ebins, loc = center, scale = width)

        # if (center, width) in [(0.75, 0.36),
        #                        (1.29, 1.6400000000000001)]:
        xs_rat, xs_rat_statErr = fit_injected_target(target,
                                                     injTarOutFile = 'plots/xs_inj_target_plotDump/inj_target_'+str(round(center, 3))+'_'+str(round(width, 3))+'.png',
                                                     secTarOutFile = 'plots/xs_secondary_target_plotDump/secondary_target_'+str(round(center, 3))+'_'+str(round(width, 3))+'.png',
                                                     xsRatioOutFile = 'plots/xs_ratio_measured_plotDump/xs_ratio_'+str(round(center, 3))+'_'+str(round(width, 3))+'.png',
                                                     secFitOutFile = 'plots/xs_secondary_fit_plotDump/sec_fit__'+str(round(center, 3))+'_'+str(round(width, 3))+'.png',
                                                     secRatesOutFile = 'plots/xs_secondary_rates_plotDump/sec_rates_'+str(round(center, 3))+'_'+str(round(width, 3))+'.png',
                                                     )

        plt.close('all')
        # else:
        # xs_rat, xs_rat_statErr = fit_injected_target_noplots(target)

        # print (xs_rat_statErr)

        allNans = np.all(np.isnan(xs_rat_statErr))
        hasZero = np.any(xs_rat_statErr == 0)
        if not (allNans or hasZero): # is the result meaningful? sometimes they are all nan or zero
            for i in range(len(nueFitter.Ebins)):
                # iterate over energy bins

                if abs(xs_rat_statErr[i]) < minErrs[i]:
                    minErrs[i] = abs(xs_rat_statErr[i])
                    minErrs_val[i] = xs_rat[i]

                    bestLoc[i] = center
                    bestWid[i] = width

                    print ("this fit is the current min for bin", i, "with statErr of", xs_rat_statErr[i])
                
            # if the abs(xs_error) in this bin is smaller than the current minimum, replace the current minimum
            # and save some other info
            
    for i in range(len(nueFitter.Ebins)):
        print (nueFitter.Ebins[i], bestLoc[i], bestWid[i])
        # print (minErrs)
        # print (minErrs_val)
        # print (bestLoc)
        # print (bestWid)
        
    # for xs_i, xsE_i, E_i, truth_i in zip(xs.T, xsErr.T, nueFitter.Ebins, xs_true):
    #     minResid = 1.e10
    #     minErr = None
    #     minErr_val = None
    #     for xs_j, xsE_j in zip(xs_i, xsE_i):
    #         if abs(xs_j - truth_i) < minResid:
    #             minResid = abs(xs_j - truth_i)
    #             minErr = abs(xsE_j)
    #             minErr_val = xs_j
    #     print ("Ebins: ", E_i, " most confident measure: ", minErr_val, "+/-", minErr)
    #     print ("True value: ", )
    #     if minErr < 0:
    #         print(xsE_i)
    #     xs_report.append(minErr_val)
    #     xsErr_report.append(minErr)

    # xs_report = np.array(xs_report)
    # xsErr_report = np.array(xsErr_report)
        
    # print(np.mean(xsErr_report[nueFitter.Ebins < 5]))
    # print(np.mean((np.abs(xs_report - xs_true)/xs_true)[nueFitter.Ebins < 5]))

    xrp = xs_ratio_plot()
    xrp.add(np.array(minErrs_val), np.array(minErrs), label = 'Combined Measurement')

    plt.figure()
    plt.scatter(nueFitter.Ebins, bestLoc)
    plt.title(r'Best Location')
    plt.xlabel(r'Measurement Location')
    plt.ylabel(r'Injected Target Center')

    plt.figure()
    plt.scatter(nueFitter.Ebins, bestWid)
    plt.title(r'Best Width')
    plt.xlabel(r'Measurement Location')
    plt.ylabel(r'Injected Target Width')
    
    plt.show()
        

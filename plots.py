import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os

from utils import *
from fluxes import *

DUNEblue = '#7FAED5'
DUNElightOrange = '#F19E54'
DUNEdarkOrange = '#F0652B'
DUNEgreen = '#8ACA6F'
DUNEgray = '#626466'
DUNEyellow = '#FBD03F'
DUNEpurple = '#5B3069'
DUNElightPurple = '#8C6E96'
DUNEcyan = '#42C2A8'
DUNEpink = '#F7AEC2'

DUNEcolors = [DUNEblue,
              DUNElightOrange,
              DUNEgreen,
              DUNEdarkOrange,
              DUNEyellow,
              DUNEpink,
              DUNEpurple,
              DUNEcyan,
              DUNEgray,
              DUNElightPurple,
]

LaTeXflavor = {"numu": r'$\nu_\mu$',
               "numubar": r'$\bar{\nu}_\mu$',
               "nue": r'$\nu_e$',
               "nuebar": r'$\bar{\nu}_e$',
               "nutau": r'$\nu_\tau$',
               "nutaubar": r'$\bar{\nu}_\tau$'}

matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
matplotlib.rc('text', usetex = True)
matplotlib.rc('axes', prop_cycle = matplotlib.cycler(color = DUNEcolors))

class plot (object):
    def __init__(self, style = "plot", *args, **kwargs):
        self.style = style
        self.legLineList = []
        self.legLabelList = []
    def show(self, *args, **kwargs):
        self.fig.show(*args, **kwargs)

    def save(self, fileName, remake = True, *args, **kwargs):
        if os.path.exists and not remake:
            pass
        else:
            self.fig.savefig(fileName, *args, **kwargs)
    def plot(self, ax, *args, **kwargs):
        if self.style == "plot":

            return ax.plot(*args, **kwargs)
        elif self.style == "step":
            return ax.step(*args, where = 'mid', **kwargs)
        elif self.style == "errorbar":
            return ax.errorbar(*args, **kwargs), 
        elif self.style == "errorband":
            yerr = kwargs.pop('yerr')
            lower = args[1] - yerr
            upper = args[1] + yerr
            band = ax.fill_between(args[0], lower, upper, alpha = 0.1)
            return ax.plot(*args, **kwargs)
        elif self.style == "errorbandstep":
            yerr = kwargs.pop('yerr')
            lower = args[1] - yerr
            upper = args[1] + yerr
            line =  ax.step(*args, where = 'mid', **kwargs)
            if 'color' in kwargs:
                kwargs.pop('color')
            band = ax.fill_between(args[0], lower, upper, step = 'mid', alpha = 0.5, color = line[0].get_color(), **kwargs)
            return line
 
class fit_and_ratio_plot (plot):
    def __init__(self, fitter = None, useTarget = True,
                 useFit = True, title = None, Ebounds = True,
                 xlim = (0, 10),
                 *args, **kwargs):
        super(fit_and_ratio_plot, self).__init__(*args, **kwargs)
        bandBounds = (0.16, 0.84)

        self.fig = plt.figure(figsize = (8.5, 5))
        gs = GridSpec(2, 1,
                      figure = self.fig,
                      height_ratios = [0.7, 0.3],
                      hspace = 0)
        self.axUp = self.fig.add_subplot(gs[0, :])
        self.axLo = self.fig.add_subplot(gs[1, :])

        self.legArgs = {"frameon": True,
                        "loc": "upper right"}
        if title:
            self.legArgs.update({"title": title})
            
        if fitter:
            if useTarget:
                self.add_target(fitter, Ebounds)
            if useFit:
                self.add(fitter, **kwargs)
            
        self.axUp.set_xlim(*xlim)
        self.axLo.set_xlim(*xlim)

        self.axUp.set_xticklabels([])
        self.axUp.grid(True, which = 'both')
        self.axUp.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]', labelpad = 15)

        self.axLo.grid(True)
        self.axLo.set_xlabel(r'$E_\nu$ [GeV]')
        if "ylabel" in kwargs:
            self.axLo.set_ylabel(kwargs["ylabel"], labelpad = 5)
        else:
            self.axLo.set_ylabel(r'$\frac{ND - \mathrm{target}}{\mathrm{target}}$', labelpad = 5)
        
        # self.axLo.set_ylim(-0.06, 0.06)
        # self.axLo.set_yticks([-0.04, -0.02, 0, 0.02, 0.04])
        # self.axLo.set_yticklabels(["-4\%", "-2\%", "0\%", "2\%", "4\%"])

        # self.axLo.set_ylim(-0.6, 0.6)
        # self.axLo.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
        # self.axLo.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axLo.set_ylim(-0.25, 0.25)
        self.axLo.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
        self.axLo.set_yticklabels(["-20\%", "-10\%", "0\%", "10\%", "20\%"])

        self.axLo.grid(True, which = 'both')

        self.fig.tight_layout()

        self.full_flux = True

    def add_target(self, fitter, label = None, full_flux = True,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        targetNomLine, = self.plot(self.axUp,
                                   fitter.Ebins,
                                   self.target,
                                   color = color,
                                   **kwargs)
        if not label:
            label = ''.join([r'FD ',
                             LaTeXflavor[fitter.FDfromFlavor],
                             r' $\rightarrow$ ',
                             LaTeXflavor[fitter.FDtoFlavor]])
        if label:
            self.legLineList.append(targetNomLine)
            self.legLabelList.append(label)

        if Ebounds:
            self.axUp.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.axUp.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')
            self.axUp.arrow(fitter.Ebounds[0], 3.85e-16,
                            0.15, 0,
                            width = 2.e-18,
                            head_length = 0.05,
                            color = 'red')
            self.axUp.arrow(fitter.Ebounds[1], 0.35e-16,
                            -0.15, 0,
                            width = 2.e-18,
                            head_length = 0.05,
                            color = 'red')
        
            self.axLo.axvline(x = fitter.Ebounds[0],
                          ls = '--',
                              color = 'red')
            self.axLo.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')

        self.axUp.set_ylim(-0.2*np.max(self.target),
                           1.2*np.max(self.target))
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)
        
    def add(self, fitter, label = None, color = None, **kwargs):
        NDNomLine, = self.plot(self.axUp,
                               fitter.Ebins,
                               fitter.fluxPred,
                               color = color,
                               **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(None)
        else:
            self.legLabelList.append(label)

        if self.full_flux:
            target = fitter.FD_oscillated
            # denom = fitter.FD_oscillated
            denom = fitter.FD_unoscillated
        else:
            target = self.target
            denom = self.target
        self.plot(self.axLo,
                  fitter.Ebins,
                  (fitter.fluxPred - target)/denom,
                  color = NDNomLine.get_color())
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)

class fit_and_ratio_rate_plot (plot):
    def __init__(self, fitter = None, useTarget = True,
                 useFit = True, title = None, Ebounds = True,
                 xlim = (0, 10),
                 *args, **kwargs):
        super(fit_and_ratio_rate_plot, self).__init__(*args, **kwargs)
        bandBounds = (0.16, 0.84)

        self.fig = plt.figure(figsize = (8.5, 5))
        gs = GridSpec(2, 1,
                      figure = self.fig,
                      height_ratios = [0.7, 0.3],
                      hspace = 0)
        self.axUp = self.fig.add_subplot(gs[0, :])
        self.axLo = self.fig.add_subplot(gs[1, :])

        self.legArgs = {"frameon": True,
                        "loc": "upper right"}
        if title:
            self.legArgs.update({"title": title})
            
        if fitter:
            if useTarget:
                self.add_target(fitter, Ebounds)
            if useFit:
                self.add(fitter, **kwargs)

        self.axUp.set_xlim(*xlim)
        self.axLo.set_xlim(*xlim)
            
        self.axUp.set_xticklabels([])
        self.axUp.grid(True, which = 'both')
        self.axUp.set_ylabel(r'Event Rate [events per bin]', labelpad = 15)

        self.axLo.grid(True)
        self.axLo.set_xlabel(r'$E_\nu$ [GeV]')
        if "ylabel" in kwargs:
            self.axLo.set_ylabel(kwargs["ylabel"], labelpad = 5)
        else:
            self.axLo.set_ylabel(r'$\frac{ND - \mathrm{target}}{\mathrm{target}}$', labelpad = 5)
        
        # self.axLo.set_ylim(-0.60, 0.60)
        # self.axLo.set_yticks([-0.50, -0.25, 0, 0.25, 0.50])
        # self.axLo.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axLo.set_ylim(-0.25, 0.25)
        self.axLo.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
        self.axLo.set_yticklabels(["-20\%", "-10\%", "0\%", "10\%", "20\%"])
        self.axLo.grid(True, which = 'both')

        self.fig.tight_layout()

        self.full_flux = True

    def add_target(self, fitter, label = None, full_flux = True, scatter = False,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        if self.style in ['errorbar', 'errorband', 'errorbandstep']:
            yerr = fitter.FD_rate_statErr
        else:
            yerr = None
        if scatter:
            targetNomLine = self.axUp.errorbar(fitter.Ebins,
                                               fitter.FD_rate,
                                               yerr = yerr,
                                               xerr = 0.5*fitter.EbinWidths,
                                               color = color,
                                               ls = 'none',
                                               **kwargs)
            self.axLo.errorbar(fitter.Ebins,
                               np.zeros_like(fitter.Ebins),
                               yerr = yerr/fitter.FD_rate,
                               xerr = 0.5*fitter.EbinWidths,
                               color = color,
                               ls = 'none',
                               **kwargs)
        else:
            targetNomLine, = self.plot(self.axUp,
                                       fitter.Ebins,
                                       fitter.FD_rate,
                                       yerr = yerr,
                                       color = color,
                                       **kwargs)
        
        if not label:
            label = ''.join([r'FD ',
                             LaTeXflavor[fitter.FDfromFlavor],
                             r' $\rightarrow$ ',
                             LaTeXflavor[fitter.FDtoFlavor]])
        if label:
            self.legLineList.append(targetNomLine)
            self.legLabelList.append(label)

        if Ebounds:
            self.axUp.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.axUp.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')
            self.axUp.arrow(fitter.Ebounds[0], 0.5*self.axUp.get_ylim()[-1],
                            0.15, 0,
                            width = 0.01*self.axUp.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
            self.axUp.arrow(fitter.Ebounds[1], 0.5*self.axUp.get_ylim()[-1],
                            -0.15, 0,
                            width = 0.01*self.axUp.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
        
            self.axLo.axvline(x = fitter.Ebounds[0],
                          ls = '--',
                              color = 'red')
            self.axLo.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')

        if yerr is not None:
            self.axUp.set_ylim(-0.2*np.max(fitter.FD_rate+yerr),
                               1.2*np.max(fitter.FD_rate+yerr))
        else:
            self.axUp.set_ylim(-0.2*np.max(fitter.FD_rate),
                               1.2*np.max(fitter.FD_rate))
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)
        
    def add(self, fitter, label = None, color = None, **kwargs):
        if self.style in ['errorbar', 'errorband', 'errorbandstep']:
            yerr = fitter.ratePred_statErr
        else:
            yerr = None
        NDNomLine, = self.plot(self.axUp,
                               fitter.Ebins,
                               fitter.ratePred,
                               yerr = yerr,
                               color = color,
                               **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(None)
        else:
            self.legLabelList.append(label)

        if self.full_flux:
            target = fitter.FD_oscillated
            print ("thing!")
            denom = fitter.FD_oscillated
            # denom = fitter.FD_unoscillated
        else:
            print ("wrong thing :(")
            target = self.target
            denom = self.target
        self.plot(self.axLo,
                  fitter.Ebins,
                  (fitter.ratePred - fitter.FD_rate)/fitter.FD_rate,
                  yerr = fitter.ratePred_statErr/fitter.FD_rate,
                  color = NDNomLine.get_color())
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)

class fit_and_ratio_plot_with_sliders (plot):
    def __init__(self, fitter = None, useTarget = True, useFit = True, title = None, Ebounds = True, regLims = [-12, -3], calcKwargs = dict(), *args, **kwargs):
        super(fit_and_ratio_plot_with_sliders, self).__init__(**kwargs)
        from matplotlib.widgets import Slider, Button, RadioButtons
        
        bandBounds = (0.16, 0.84)

        self.fig = plt.figure(figsize = (8.5, 6))
        gs1 = GridSpec(2, 1,
                       figure = self.fig,
                       height_ratios = [0.7, 0.3],
                       hspace = 0,
                       top = 0.95,
                       bottom = 0.22)
        gs2 = GridSpec(3, 1,
                       figure = self.fig,
                       top = 0.15,
                       bottom = 0.05)
        self.axUp = self.fig.add_subplot(gs1[0, :])
        self.axLo = self.fig.add_subplot(gs1[1, :])
        self.axLoBoundSlider = self.fig.add_subplot(gs2[0, :])
        self.axHiBoundSlider = self.fig.add_subplot(gs2[1, :])
        self.axRegSlider = self.fig.add_subplot(gs2[2, :])
        
        self.legLineList = []
        self.ratioLineList = []
        self.legLabelList = []

        self.legArgs = {"frameon": True,
                        "loc": "upper right"}
        if title:
            self.legArgs.update({"title": title})

        self.fitters = []
        if fitter:
            if useTarget:
                self.add_target(fitter, Ebounds)
            if useFit:
                self.add(fitter, **kwargs)
            
        self.axUp.set_xlim(0, 10)
        self.axUp.set_xticklabels([])
        self.axUp.grid(True, which = 'both')
        self.axUp.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]', labelpad = 15)

        self.axLo.grid(True)
        self.axLo.set_xlim(0, 10)
        self.axLo.set_xlabel(r'$E_\nu$ [GeV]')
        if "ylabel" in kwargs:
            self.axLo.set_ylabel(kwargs["ylabel"], labelpad = 5)
        else:
            self.axLo.set_ylabel(r'$\frac{ND - FD (osc.)}{FD (unosc.)}$', labelpad = 5)
        # self.axLo.set_ylim(-0.6, 0.6)
        # self.axLo.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
        # self.axLo.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axLo.set_ylim(-0.12, 0.12)
        self.axLo.set_yticks([-0.5, -0.25, 0, 0.25, 0.50])
        self.axLo.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axLo.grid(True, which = 'both')

        # self.axRegSlider = plt.axes([0.1, 0.1, 0.65, 0.03])
        self.regLims = regLims
        self.sReg = Sliderlog(self.axRegSlider,
                              r'$\lambda$',
                              regLims[0], regLims[1],
                              valinit = regLims[0])
        ticks = range(regLims[0], regLims[1]+1)[::2]
        self.axRegSlider.set_xticks(ticks)
        self.axRegSlider.set_xticklabels([r'$10^{'+str(i)+r'}$' for i in ticks])

        self.sLoBound = Slider(self.axLoBoundSlider,
                               r'$E_{min}$',
                               0, 10, valinit = 0)
        self.sHiBound = Slider(self.axHiBoundSlider,
                               r'$E_{max}$',
                               0, 10, valinit = 10)

        self.sReg.on_changed(self.update)
        self.sLoBound.on_changed(self.update)
        self.sHiBound.on_changed(self.update)

        # self.axLoBoundSlider.set

        self.calcKwargs = calcKwargs
        
        self.fig.tight_layout()

        self.full_flux = True

    def add_target(self, fitter, label = None, full_flux = True,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        targetNomLine, = self.plot(self.axUp,
                                   fitter.Ebins,
                                   self.target,
                                   color = color,
                                   **kwargs)
        # if not label:
        #     label = ''.join([r'FD ',
        #                      LaTeXflavor[fitter.FDfromFlavor],
        #                      r' $\rightarrow$ ',
        #                      LaTeXflavor[fitter.FDtoFlavor]])
        if label:
            self.legLineList.append(targetNomLine)
            self.legLabelList.append(label)

        if Ebounds:
            self.loBoundLineUp = self.axUp.axvline(x = fitter.Ebounds[0],
                                                   ls = '--',
                                                   color = 'red')
            self.hiBoundLineUp = self.axUp.axvline(x = fitter.Ebounds[1],
                                                   ls = '--',
                                                   color = 'red')
            # self.loBoundArrowUp = self.axUp.arrow(fitter.Ebounds[0], 3.85e-16,
            #                                       0.15, 0,
            #                                       width = 2.e-18,
            #                                       head_length = 0.05,
            #                                       color = 'red')
            # self.hiBoundArrowUp = self.axUp.arrow(fitter.Ebounds[1], 0.35e-16,
            #                                       -0.15, 0,
            #                                       width = 2.e-18,
            #                                       head_length = 0.05,
            #                                       color = 'red')
        
            self.loBoundLineDown = self.axLo.axvline(x = fitter.Ebounds[0],
                                                     ls = '--',
                                                     color = 'red')
            self.hiBoundLineDown = self.axLo.axvline(x = fitter.Ebounds[1],
                                                     ls = '--',
                                                     color = 'red')

        self.axUp.set_ylim(-0.2*np.max(self.target),
                           1.2*np.max(self.target))
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)
        
    def add(self, fitter, label = None, **kwargs):
        self.fitters.append(fitter)
        
        NDNomLine, = self.plot(self.axUp,
                               fitter.Ebins,
                               fitter.fluxPred,
                               **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(None)
        else:
            self.legLabelList.append(label)

        if self.full_flux:
            target = fitter.FD_oscillated
            # denom = fitter.FD_oscillated
            denom = fitter.FD_unoscillated
        else:
            target = self.target
            denom = self.target
        RatioLine, = self.plot(self.axLo,
                               fitter.Ebins,
                               (fitter.fluxPred - target)/denom,
                               **kwargs)
        self.ratioLineList.append(RatioLine)
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)

    def update(self, val):
        reg = self.sReg.val
        loBound = self.sLoBound.val
        hiBound = self.sHiBound.val
        for i, fitter in enumerate(self.fitters):
            if self.full_flux:
                target = fitter.FD_oscillated
                # denom = fitter.FD_oscillated
                denom = fitter.FD_unoscillated
            else:
                target = self.target
                denom = self.target

            fitter.set_fit_region(energies = (loBound, hiBound))
            fitter.calc_coeffs(reg, reg, **self.calcKwargs)
            self.legLineList[i+1].set_ydata(fitter.fluxPred)
            self.ratioLineList[i].set_ydata((fitter.fluxPred - target)/denom)
        # self.loBoundLineUp.set_xdata([loBound, loBound])
        # self.loBoundLineDown.set_xdata([loBound, loBound])
        # self.hiBoundLineUp.set_xdata([hiBound, hiBound])
        # self.hiBoundLineDown.set_xdata([hiBound, hiBound])
            
        self.fig.canvas.draw_idle()

class slider_super_plot (plot):
    def __init__(self, fitter = None, useTarget = True, useFit = True, title = None, Ebounds = True, regLims = [-12, -3], calcKwargs = dict(), **kwargs):
        super(slider_super_plot, self).__init__(**kwargs)
        from matplotlib.widgets import Slider, Button, RadioButtons
        
        bandBounds = (0.16, 0.84)

        self.fig = plt.figure(figsize = (11, 8.5))
        gsFitAndRatio = GridSpec(2, 1,
                                 figure = self.fig,
                                 height_ratios = [0.7, 0.3],
                                 hspace = 0,
                                 top = 0.95,
                                 bottom = 0.22,
                                 left = 0.05,
                                 right = 0.45)
        gsESliders = GridSpec(2, 1,
                              figure = self.fig,
                              top = 0.15,
                              bottom = 0.05,
                              left = 0.05,
                              right = 0.45)
        gsCoeffs = GridSpec(1, 2,
                            figure = self.fig,
                            wspace = 0,
                            width_ratios = [0.6, 0.4],
                            top = 0.95,
                            bottom = 0.65,
                            left = 0.55,
                            right = 0.95)
        gsLcurve = GridSpec(2, 1,
                            figure = self.fig,
                            top = 0.55,
                            bottom = 0.22,
                            left = 0.55,
                            right = 0.95)
        gsRegSliders = GridSpec(1, 1,
                                figure = self.fig,
                                top = 0.15,
                                bottom = 0.05,
                                left = 0.55,
                                right = 0.95)

        
        self.axFit = self.fig.add_subplot(gsFitAndRatio[0, :])
        self.axRatio = self.fig.add_subplot(gsFitAndRatio[1, :])

        self.axOAcoeff = self.fig.add_subplot(gsCoeffs[:,0])
        self.axHCcoeff = self.fig.add_subplot(gsCoeffs[:,1])

        self.axLcurve = self.fig.add_subplot(gsLcurve[0,:])
        self.axCurvature = self.fig.add_subplot(gsLcurve[1,:])
        
        self.axLoBoundSlider = self.fig.add_subplot(gsESliders[0, :])
        self.axHiBoundSlider = self.fig.add_subplot(gsESliders[1, :])

        self.axRegSlider = self.fig.add_subplot(gsRegSliders[0, :])
        
        self.legLineList = []
        self.ratioLineList = []
        self.OAcoeffLineList = []
        self.HCcoeffLineList = []
        self.LcurveLineList = []
        self.curvatureLineList = []
        self.legLabelList = []

        self.legArgs = {"frameon": True,
                        "loc": "upper right"}
        if title:
            self.legArgs.update({"title": title})

        self.fitters = []
        # if fitter:
        #     if useTarget:
        #         self.add_target(fitter, Ebounds)
        #     if useFit:
        #         self.add(fitter, **kwargs)
            
        self.axFit.set_xlim(0, 10)
        self.axFit.set_xticklabels([])
        self.axFit.grid(True, which = 'both')
        self.axFit.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]', labelpad = 15)

        self.axRatio.grid(True)
        self.axRatio.set_xlim(0, 10)
        self.axRatio.set_xlabel(r'$E_\nu$ [GeV]')
        if "ylabel" in kwargs:
            self.axRatio.set_ylabel(kwargs["ylabel"], labelpad = 5)
        else:
            self.axRatio.set_ylabel(r'$\frac{ND - FD (osc.)}{FD (unosc.)}$', labelpad = 5)
        # self.axRatio.set_ylim(-0.6, 0.6)
        # self.axRatio.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
        # self.axRatio.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axRatio.set_ylim(-0.12, 0.12)
        self.axRatio.set_yticks([-0.5, -0.25, 0, 0.25, 0.50])
        self.axRatio.set_yticklabels(["-50\%", "-25\%", "0\%", "25\%", "50\%"])
        self.axRatio.grid(True, which = 'both')

        coeffYlim = 2.5e-7
        self.axOAcoeff.set_ylim(-coeffYlim, coeffYlim)
        self.axOAcoeff.grid(True)
        self.axOAcoeff.set_xlabel(r'$D_{OA}$ [m]')
        self.axOAcoeff.set_ylabel(r'$c_i$')
        self.axOAcoeff.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))

        self.axHCcoeff.set_ylim(-coeffYlim, coeffYlim)
        self.axHCcoeff.grid(True)
        self.axHCcoeff.set_xlabel(r'Horn Current [kA]')
        self.axHCcoeff.set_yticklabels([])
        self.axHCcoeff.ticklabel_format(axis = 'x', style = 'plain', useOffset = False)


        # self.axRegSlider = plt.axes([0.1, 0.1, 0.65, 0.03])
        self.regRange = np.logspace(-12, -4, 1000)
        self.regLims = [self.regRange[0], self.regRange[-1]]
        self.sReg = Sliderlog(self.axRegSlider,
                              r'$\lambda$',
                              regLims[0], regLims[1],
                              valinit = regLims[0])
        ticks = range(regLims[0], regLims[1]+1)[::2]
        self.axRegSlider.set_xticks(ticks)
        self.axRegSlider.set_xticklabels([r'$10^{'+str(i)+r'}$' for i in ticks])

        self.sLoBound = Slider(self.axLoBoundSlider,
                               r'$E_{min}$',
                               0, 10, valinit = 0)
        self.sHiBound = Slider(self.axHiBoundSlider,
                               r'$E_{max}$',
                               0, 10, valinit = 10)

        self.sReg.on_changed(self.updateEsliders)
        self.sLoBound.on_changed(self.updateEsliders)
        self.sHiBound.on_changed(self.updateRegSliders)

        # self.axLoBoundSlider.set

        self.fig.tight_layout()
        self.updateEsliders('dumb')

        self.full_flux = True

        self.calcKwargs = calcKwargs

        if fitter:
            if useTarget:
                self.add_target(fitter, Ebounds)
            if useFit:
                self.add(fitter, **kwargs)
                

    def add_target(self, fitter, label = None, full_flux = True,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        targetNomLine, = self.plot(self.axFit,
                                   fitter.Ebins,
                                   self.target,
                                   color = color,
                                   **kwargs)
        # if not label:
        #     label = ''.join([r'FD ',
        #                      LaTeXflavor[fitter.FDfromFlavor],
        #                      r' $\rightarrow$ ',
        #                      LaTeXflavor[fitter.FDtoFlavor]])
        if label:
            self.legLineList.append(targetNomLine)
            self.legLabelList.append(label)

        if Ebounds:
            self.loBoundLineUp = self.axFit.axvline(x = fitter.Ebounds[0],
                                                    ls = '--',
                                                    color = 'red')
            self.hiBoundLineUp = self.axFit.axvline(x = fitter.Ebounds[1],
                                                    ls = '--',
                                                    color = 'red')
            # self.loBoundArrowUp = self.axFit.arrow(fitter.Ebounds[0], 3.85e-16,
            #                                       0.15, 0,
            #                                       width = 2.e-18,
            #                                       head_length = 0.05,
            #                                       color = 'red')
            # self.hiBoundArrowUp = self.axFit.arrow(fitter.Ebounds[1], 0.35e-16,
            #                                       -0.15, 0,
            #                                       width = 2.e-18,
            #                                       head_length = 0.05,
            #                                       color = 'red')
        
            self.loBoundLineDown = self.axRatio.axvline(x = fitter.Ebounds[0],
                                                     ls = '--',
                                                     color = 'red')
            self.hiBoundLineDown = self.axRatio.axvline(x = fitter.Ebounds[1],
                                                     ls = '--',
                                                     color = 'red')

        self.axFit.set_ylim(-0.2*np.max(self.target),
                           1.2*np.max(self.target))
        
        self.axFit.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)
        
    def add(self, fitter, label = None, color = None, **kwargs):
        self.fitters.append(fitter)

        self.addFit(fitter, color, label)
        self.addRatio(fitter, color)
        self.addCoeffs(fitter, color)
        self.addLcurves(fitter, color)

    def addFit(self, fitter, color, label, **kwargs):
        NDNomLine, = self.plot(self.axFit,
                               fitter.Ebins,
                               fitter.fluxPred,
                               color = color,
                               **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(None)
        else:
            self.legLabelList.append(label)

    def addRatio(self, fitter, color):
        if self.full_flux:
            target = fitter.FD_oscillated
            # denom = fitter.FD_oscillated
            denom = fitter.FD_unoscillated
        else:
            target = self.target
            denom = self.target
        RatioLine, = self.plot(self.axRatio,
                               fitter.Ebins,
                               (fitter.fluxPred - target)/denom,
                               color = color)
        self.ratioLineList.append(RatioLine)
        
        self.axFit.legend(self.legLineList,
                          self.legLabelList,
                          **self.legArgs)

    def addCoeffs(self, fitter, color):
        OAcoeffLine, = self.axOAcoeff.step(fitter.OAbins[fitter.OAbins <= fitter.maxOA],
                                           fitter.cOA,
                                           color = color,
                                           where = 'mid')
        self.OAcoeffLineList.append(OAcoeffLine)

        HCcoeffLine = self.axHCcoeff.scatter(fitter.HCbins,
                                             fitter.cHC,
                                             color = color,
                                             marker = "+",
                                             s = 500)
        self.HCcoeffLineList.append(HCcoeffLine)
        
    def addLcurves(self, fitter, color, **kwargs):
        LcurveLine, = self.axLcurve.plot(np.zeros_like(self.regRange), np.zeros_like(self.regRange))
        self.LcurveLineList.append(LcurveLine)
        
        curvatureLine, = self.axCurvature.plot(self.regRange[1:-1], np.zeros_like(self.regRange[1:-1]))
        self.curvatureLineList.append(curvatureLine)

    def updateLcurves(self, fitter, **kwargs):
        res = []
        sol = []
        for reg in self.regRange:
            if not "HCreg" in kwargs:
                HCreg = reg
            else:
                HCreg = kwargs["HCreg"]

            fitter.calc_coeffs(reg, HCreg, **self.calcKwargs)
            res.append(fitter.residual_norm(**kwargs))
            sol.append(fitter.solution_norm(**kwargs))

        self.axLcurve.plot(res, sol)
            
        res = np.array(res)
        sol = np.array(sol)
        dl = np.diff(self.regRange)
        xi = np.log(sol)
        rho = np.log(res)
        xi_prime = np.diff(xi)/dl
        rho_prime = np.diff(rho)/dl
        xi_prime_prime = np.diff(xi_prime)/dl[:-1]
        rho_prime_prime = np.diff(rho_prime)/dl[:-1]

        curv = 2*(rho_prime[:-1]*xi_prime_prime - rho_prime_prime*xi_prime[:-1])/np.power(np.power(rho_prime[:-1], 2) + np.power(xi_prime[:-1], 2), 3./2)
    
        curvatureLine, = self.axCurvature.plot(self.regRange[1:-1], curv)
        self.curvatureLineList.append(curvatureLine)
 
    def updateRegSliders(self, val, **kwargs):
        reg = self.sReg.val
        loBound = self.sLoBound.val
        hiBound = self.sHiBound.val
        for i, fitter in enumerate(self.fitters):
            if self.full_flux:
                target = fitter.FD_oscillated
                # denom = fitter.FD_oscillated
                denom = fitter.FD_unoscillated
            else:
                target = self.target
                denom = self.target

            fitter.set_fit_region(energies = (loBound, hiBound))
            fitter.calc_coeffs(reg, reg, **self.calcKwargs)
            self.legLineList[i+1].set_ydata(fitter.fluxPred)
            self.ratioLineList[i].set_ydata((fitter.fluxPred - target)/denom)
            self.OAcoeffLineList[i].set_ydata(fitter.cOA)
            self.HCcoeffLineList[i].set_offsets(np.array([fitter.HCbins, fitter.cHC]).T)

        # self.loBoundLineUp.set_xdata([loBound, loBound])
        # self.loBoundLineDown.set_xdata([loBound, loBound])
        # self.hiBoundLineUp.set_xdata([hiBound, hiBound])
        # self.hiBoundLineDown.set_xdata([hiBound, hiBound])
            
        self.fig.canvas.draw_idle()

    def updateEsliders(self, val):
        reg = self.sReg.val
        loBound = self.sLoBound.val
        hiBound = self.sHiBound.val
        for i, fitter in enumerate(self.fitters):
            if self.full_flux:
                target = fitter.FD_oscillated
                # denom = fitter.FD_oscillated
                denom = fitter.FD_unoscillated
            else:
                target = self.target
                denom = self.target

            fitter.set_fit_region(energies = (loBound, hiBound))
            fitter.calc_coeffs(reg, reg, **self.calcKwargs)
            self.legLineList[i+1].set_ydata(fitter.fluxPred)
            self.ratioLineList[i].set_ydata((fitter.fluxPred - target)/denom)
            self.OAcoeffLineList[i].set_ydata(fitter.cOA)
            self.HCcoeffLineList[i].set_offsets(np.array([fitter.HCbins, fitter.cHC]).T)

        # self.loBoundLineUp.set_xdata([loBound, loBound])
        # self.loBoundLineDown.set_xdata([loBound, loBound])
        # self.hiBoundLineUp.set_xdata([hiBound, hiBound])
        # self.hiBoundLineDown.set_xdata([hiBound, hiBound])
            
        self.fig.canvas.draw_idle()

class coeff_plot (plot):
    def __init__(self, fitter = None, title = "Coefficients", HC = True, legend = True, style = "step", **kwargs):
        super(coeff_plot, self).__init__(style = style, **kwargs)
        self.fig = plt.figure()

        self.HC = HC
        if self.HC:
            gs = GridSpec(1, 2,
                          figure = self.fig,
                          wspace = 0,
                          width_ratios = [0.6, 0.4])
            self.OAax = self.fig.add_subplot(gs[:, 0])
            self.HCax = self.fig.add_subplot(gs[:, 1])

        else:
            self.OAax = self.fig.gca()

        if "ylim" in kwargs:
            self.ylim = kwargs["ylim"]
        else:
            self.ylim = 2.5e-7

        self.OAax.set_ylim(-self.ylim, self.ylim)
        self.OAax.grid(True)
        self.OAax.set_xlabel(r'$D_{OA}$ [m]')
        self.OAax.set_ylabel(r'$c_i$')
        self.OAax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))

        if self.HC:
            self.HCax.set_ylim(-self.ylim, self.ylim)
            self.HCax.grid(True)
            self.HCax.set_xlabel(r'Horn Current [kA]')
            self.HCax.set_yticklabels([])
            self.HCax.ticklabel_format(axis = 'x', style = 'plain', useOffset = False)

        self.fig.suptitle(title)
        if "legendLoc" in kwargs:
            self.legendLoc = kwargs["legendLoc"]
        else:
            self.legendLoc = "left"

        self.legArgs = {}
        if "fontSize" in kwargs:
            self.legArgs.update({"prop": {"size": kwargs["fontSize"]}})
        self.legend = legend
            
        self.fig.tight_layout()

        self.coeffs = []

        if fitter:
            self.add(fitter, **kwargs)

    def add(self, fitter, label = None, color = None, HCplot = False, **kwargs):

        if self.HC:
            self.HCax.set_xticks(fitter.HCbins)
            self.HCax.set_xticklabels([str(HC) for HC in fitter.HCbins])

        self.plot(self.OAax,
                  fitter.OAbins,
                  fitter.cOA,
                  color = color,
                  label = label,
                  **kwargs)
        
        self.coeffs.append(fitter.cOA[fitter.OAbins <= fitter.maxOA])

        if self.HC:
            if HCplot:
                self.HCax.plot(fitter.HCbins,
                               fitter.cHC,
                               color = color,
                               label = label)
            else:
                self.HCax.scatter(fitter.HCbins,
                                  fitter.cHC,
                                  color = color,
                                  marker = "+",
                                  s = 500,
                                  label = label)
            self.coeffs.append(fitter.cHC)


        self.ylim = 1.2*max([np.max(np.abs(c)) for c in self.coeffs if np.any(c)])
        self.OAax.set_ylim(-self.ylim, self.ylim)
        if self.HC:
            self.HCax.set_ylim(-self.ylim, self.ylim)
        
        if self.legend:
            if label:
                if self.legendLoc == "left":
                    self.OAax.legend(**self.legArgs)
                elif self.legendLoc == "right":
                    self.HCax.legend(**self.legArgs)

class ND_flux_plot (plot):
    def __init__(self, fitter, title = "ND Flux", **kwargs):
        super(ND_flux_plot, self).__init__(**kwargs)
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        
        self.cmesh = self.ax.pcolormesh(fitter.Ebins,
                                        fitter.OAbins[fitter.OAbins <= fitter.maxOA],
                                        fitter.ND_OA.T)

        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$D_{OA}$ [m]')

        self.cb = plt.colorbar(mappable = self.cmesh)
        self.cb.set_label(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        self.ax.set_title(title)
        self.fig.tight_layout()

class ND_flux_slice_plot (plot):
    def __init__(self, fitter, slices = [0, 10, 20, 30, 40, 50, 60], **kwargs):
        super(ND_flux_slice_plot, self).__init__(**kwargs)
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        
        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        meanLines = []
        errLines = []
        labels = []
        for sliceInd in slices:
            if "ND_univs" in dir(fitter):
                thisSlice = fitter.ND_univs[:,:,sliceInd]
                lower, upper = np.quantile(thisSlice, [0.16, 0.84], axis = 0)
                errLine = self.ax.fill_between(fitter.Ebins,
                                               lower,
                                               upper,
                                               alpha = 0.5)
                errLines.append(errLine)
            # meanLine, = self.ax.plot(fitter.Ebins,
            #                          fitter.ND[:,sliceInd])
            meanLine, = self.plot(self.ax,
                                  fitter.Ebins,
                                  fitter.ND[:,sliceInd])
            labels.append(str(fitter.OAbins[sliceInd]) + " m")
            meanLines.append(meanLine)

        if errLines:
            self.ax.legend([(errLine, meanLine) for errLine, meanLine
                            in zip(errLines, meanLines)],
                           labels,
                           frameon = False)
        else:
            self.ax.legend(meanLines,
                           labels,
                           frameon = False)
            
        # self.ax.set_xlim(np.min(fitter.Ebins), np.max(fitter.Ebins))
        self.ax.set_xlim(0, 5)

        if "title" in kwargs:
            self.ax.set_title(kwargs["title"])
        else:
            self.ax.set_title("ND Flux")
        self.fig.tight_layout()

class ND_ER_slice_plot (plot):
    def __init__(self, fitter, slices = [0, 10, 20, 30, 40, 50, 60], **kwargs):
        super(ND_ER_slice_plot, self).__init__(**kwargs)
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        
        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'Event Rate [events per bin]')

        meanLines = []
        errLines = []
        labels = []
        for sliceInd in slices:
            # if "ND_univs" in dir(fitter):
            #     thisSlice = fitter.ND_univs[:,:,sliceInd]
            #     lower, upper = np.quantile(thisSlice, [0.16, 0.84], axis = 0)
            #     errLine = self.ax.fill_between(fitter.Ebins,
            #                                    lower,
            #                                    upper,
            #                                    alpha = 0.5)
            #     errLines.append(errLine)
            meanLine, = self.plot(self.ax,
                                  fitter.Ebins,
                                  fitter.ND_rate[:,sliceInd],
                                  yerr = fitter.ND_rate_statErr[:,sliceInd])
            labels.append(str(fitter.OAbins[sliceInd]) + " m")
            meanLines.append(meanLine)

        if errLines:
            self.ax.legend([(errLine, meanLine) for errLine, meanLine
                            in zip(errLines, meanLines)],
                           labels,
                           frameon = False)
        else:
            self.ax.legend(meanLines,
                           labels,
                           frameon = False,
                           ncol = 3)
            
        # self.ax.set_xlim(np.min(fitter.Ebins), np.max(fitter.Ebins))
        self.ax.set_xlim(0, 5)
        self.ax.semilogy()
        
        if "title" in kwargs:
            self.ax.set_title(kwargs["title"])
        else:
            self.ax.set_title("ND Flux")
        self.fig.tight_layout()

class FD_flux_plot (plot):
    def __init__(self, fitter = None, title = "FD Oscillated Flux",
                 aspect = None, figSize = (6.4, 4.8), legendOn = True,
                 xlim = (0, 10), **kwargs):
        super(FD_flux_plot, self).__init__(**kwargs)
        self.fig = plt.figure(figsize = figSize)
        self.ax = self.fig.gca()

        self.legendOn = legendOn

        self.legLabelList = []
        self.legLineList = []

        self.legArgs = {"frameon": False}
        if "legCols" in kwargs:
            self.legArgs.update({"ncol": kwargs["legCols"]})
        if "xlim" in kwargs:
            self.ax.set_xlim(*kwargs["xlim"])
        else:
            self.ax.set_xlim(0, 10)
        if "ymax" in kwargs:
            self.ymax = kwargs["ymax"]
        else:
            self.ymax = None

        if fitter:
             self.add(fitter, **kwargs)

        self.ax.set_xlim(*xlim)

        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        self.ax.set_title(title)

        # self.ax.grid()
        self.ax.set_title(title)
        self.fig.tight_layout()
     
    def add(self, fitter, label = None, color = None, Ebounds = False, **kwargs):
        # if not label:
        #     label = r'FD ' + LaTeXflavor[fitter.FDfromFlavor] + r' $\rightarrow$ ' + LaTeXflavor[fitter.FDtoFlavor]
        self.legLabelList.append(label)
        
        if "FD_unosc_univs" in dir(fitter):
            oscillated = fitter.FD_unosc_univs*self.Posc
            lower, upper = np.quantile(oscillated, [0.16, 0.84], axis = 0)
            self.ax.fill_between(fitter.Ebins,
                                 lower,
                                 upper,
                                 alpha = 0.5,
                                 color = color)
            
        line, = self.plot(self.ax,
                          fitter.Ebins,
                          fitter.target,
                          color = color)
        
        self.legLineList.append(line)

        if Ebounds:
            self.ax.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.ax.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')
            self.ax.arrow(fitter.Ebounds[0], 0.5*self.ax.get_ylim()[-1],
                            0.15, 0,
                            width = 0.01*self.ax.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
            self.ax.arrow(fitter.Ebounds[1], 0.5*self.ax.get_ylim()[-1],
                            -0.15, 0,
                            width = 0.01*self.ax.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
     
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

        if self.ymax:
            self.ymax = max(self.ymax, 1.2*np.max(fitter.target))
            self.ax.set_ylim(0, self.ymax)

        return line

    def add_unosc(self, fitter, label = r'FD Unoscillated', color = None):
        self.legLabelList.append(label)
        
        if "FD_unosc_univs" in dir(fitter):
            lower, upper = np.quantile(fitter.FD_unosc_univs, [0.16, 0.84], axis = 0)
            self.ax.fill_between(fitter.Ebins,
                                 lower,
                                 upper,
                                 alpha = 0.5,
                                 color = color)
            
        line, = self.plot(self.ax,
                          fitter.Ebins,
                          fitter.FD_unoscillated,
                          color = color)
        
        self.legLineList.append(line)
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)
        
        self.ymax = max(self.ymax, 1.2*np.max(fitter.FD_unoscillated))
        self.ax.set_ylim(0, self.ymax)

        return line
         
    def add_fit(self, fitter, color = None, linestyle = None, label = r'ND Flux Match'):
        line, = self.plot(self.ax,
                          fitter.Ebins,
                          fitter.fluxPred,
                          color = color,
                          linestyle = linestyle)
        self.legLineList.append(line)
        self.legLabelList.append(label)
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

        return line
            
class FD_rate_plot (plot):
    def __init__(self, fitter = None, title = "FD Event Rate", aspect = None,
                 figSize = (6.4, 4.8), legendOn = True, style = "errorbandstep",
                 xlim = (0, 10), **kwargs):
        super(FD_rate_plot, self).__init__(style = style, **kwargs)

        self.fig = plt.figure(figsize = figSize)
        self.ax = self.fig.gca()

        self.legendOn = legendOn

        self.legLabelList = []
        self.legLineList = []

        self.legArgs = {"frameon": False}
        if "legCols" in kwargs:
            self.legArgs.update({"ncol": kwargs["legCols"]})
        if "xlim" in kwargs:
            self.ax.set_xlim(*kwargs["xlim"])
        else:
            self.ax.set_xlim(0, 10)
        if "ymax" in kwargs:
            self.ymax = kwargs["ymax"]
        else:
            self.ymax = 0

        if fitter:
             self.add(fitter, **kwargs)

        self.ax.set_xlim(*xlim)
             
        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'Event Rate [events per bin per year]')

        self.ax.set_title(title)

        # self.ax.grid()
        self.ax.set_title(title)
        self.fig.tight_layout()
     
    def add(self, fitter, label = "FD Event Rate", scatter = False, color = 'black', Ebounds = False, **kwargs):
        self.legLabelList.append(label)

        if self.style in ['errorbar', 'errorband', 'errorbandstep']:
            yerr = fitter.FD_rate_statErr
        else:
            yerr = None
        if scatter:
            line = self.ax.errorbar(fitter.Ebins,
                                     fitter.FD_rate,
                                     xerr = 0.5*fitter.EbinWidths,
                                     yerr = yerr,
                                     color = color,
                                     ls = 'none',
                                     **kwargs)
        else:
            line, = self.plot(self.ax,
                              fitter.Ebins,
                              fitter.FD_rate,
                              yerr = yerr,
                              color = color,
                              **kwargs)
        
        self.legLineList.append(line)
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

        if Ebounds:
            self.ax.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.ax.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')
            self.ax.arrow(fitter.Ebounds[0], 0.5*self.ax.get_ylim()[-1],
                            0.15, 0,
                            width = 0.01*self.ax.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
            self.ax.arrow(fitter.Ebounds[1], 0.5*self.ax.get_ylim()[-1],
                            -0.15, 0,
                            width = 0.01*self.ax.get_ylim()[-1],
                            head_length = 0.05,
                            color = 'red')
        
        if yerr is not None:
            print (np.max(fitter.FD_rate), np.max(fitter.FD_rate+yerr))
            self.ymax = max(self.ymax, 1.2*np.max(fitter.FD_rate+yerr))
        else:
            self.ymax = max(self.ymax, 1.2*np.max(fitter.FD_rate))

        self.ax.set_ylim(0, self.ymax)

        return line

    def add_fit(self, fitter, label = r'ND Match', **kwargs):
        if self.style in ['errorbar', 'errorband', 'errorbandstep']:
            yerr = fitter.ratePred_statErr
        else:
            yerr = None
        line, = self.plot(self.ax,
                          fitter.Ebins,
                          fitter.ratePred,
                          yerr = yerr,
                          **kwargs)
        self.legLineList.append(line)
        self.legLabelList.append(label)
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

        return line
            
class FD_flux_osc_and_unosc_plot (plot):
    def __init__(self, fitter, title = "FD Flux", inset_text = None, logScale = False, figSize = (6.4, 4.8), **kwargs):
        super(FD_flux_osc_and_unosc_plot, self).__init__(**kwargs)

        self.fig = plt.figure(figsize = figSize)
        self.ax = plt.gca()

        self.legLabelList = [r'Oscillated',
                             r'Unoscillated']

        oscLine, = self.ax.plot(fitter.Ebins,
                                fitter.FD_oscillated)

        unoscLine, = self.ax.plot(fitter.Ebins,
                                  fitter.FD_unoscillated)
        
        if "FD_unosc_univs" in dir(fitter):
            oscillated = fitter.FD_unosc_univs*fitter.Posc
            lower, upper = np.quantile(oscillated,
                                       [0.16, 0.84],
                                       axis = 0)
            oscBand = self.ax.fill_between(fitter.Ebins,
                                           lower,
                                           upper,
                                           alpha = 0.5,
                                           color = DUNEblue)

            lower, upper = np.quantile(fitter.FD_unosc_univs,
                                       [0.16, 0.84],
                                       axis = 0)
            unoscBand = self.ax.fill_between(fitter.Ebins,
                                             lower,
                                             upper,
                                             alpha = 0.5,
                                             color = DUNElightOrange)

            self.legLineList = [(oscBand, oscLine),
                                (unoscBand, unoscLine)]

        else:
            self.legLineList = [oscLine,
                                unoscLine]

        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        self.ax.legend(self.legLineList,
                       self.legLabelList,
                       frameon = False)
        
        self.ax.set_xlim(np.min(fitter.Ebins), np.max(fitter.Ebins))
        self.ax.set_ylim(0, 1.2*np.max(fitter.FD_unoscillated))
        
        if logScale:
            self.ax.semilogx()
            
        self.ax.set_title(title)

        if inset_text:
            top_line_y_loc = 1.05*np.max(fitter.FD_unoscillated)
            if logScale:
                top_line_x_loc = 2.e-2
            else:
                top_line_x_loc = 6.
                
            if type(inset_text) == list:
                line_spacing = 0.15*top_line_y_loc
                for i, line in enumerate(inset_text):
                    self.fig.text(top_line_x_loc,
                                  top_line_y_loc - i*line_spacing,
                                  line)

            else:
                self.fig.text(top_line_x_loc,
                              top_line_y_loc,
                              inset_text)
        self.fig.tight_layout()

    def add_fit(self, fitter, color = None):
        fitLine, = self.ax.plot(fitter.Ebins,
                                fitter.fluxPred,
                                color = color)
        self.legLineList.append(fitLine)
        self.legLabelList.append(r'ND Flux Match')
        
        self.ax.legend(self.legLineList,
                       self.legLabelList,
                       frameon = False)
 
class mike_plot (plot):
    def __init__(self, fitter = None, varType = "ppfx", varKey = 0, label = r'ND \& FD', title = None, binEdges = [], ylim = None, **kwargs):
        super(mike_plot, self).__init__(**kwargs)
        
        self.fig = plt.figure()
        gs = GridSpec(2, 1,
                      figure = self.fig,
                      height_ratios = [0.6, 0.4],
                      hspace = 0)
        self.axUp = self.fig.add_subplot(gs[0, :])
        self.axLo = self.fig.add_subplot(gs[1, :])

        self.legLineList = []
        self.legLabelList = []
        
        self.axUp.grid()
        if title:
            self.axUp.set_title(title)
        self.axUp.set_ylabel(r'$(\Phi_V - \Phi_N)/\Phi_{FD, unosc}$')
        self.axUp.set_xlim(0, 10)
        if ylim:
            self.axUp.set_ylim(-ylim, ylim)
        self.axUp.set_xticklabels([])

        self.axLo.grid()
        self.axLo.set_xlabel(r'$E_\nu$')
        self.axLo.set_ylabel(r'Difference')
        self.axLo.set_xlim(0, 10)
        self.axLo.set_ylim(-0.06, 0.06)
        self.axLo.set_yticks([-0.05, 0, 0.05])
        self.fig.tight_layout()

        self.binEdges = binEdges

        if fitter:
            self.add(fitter, varType = varType, varKey = varKey, label = label) 

    def add(self, fitter, varType = "ppfx", varKey = 0, label = r'ND \& FD', sameColor = True):
        if varType == "ppfx":
            var_target = fitter.FD_unosc_ppfx_univs[varKey]*fitter.Posc
            ND_univ = fitter.ND_ppfx_univs[varKey]
        elif varType == "other":
            var_target = fitter.FD_unosc_other_univs[varKey]*fitter.Posc
            ND_univ = fitter.ND_other_univs[varKey]

        ND = (np.dot(ND_univ, fitter.c) - fitter.fluxPred)/fitter.FD_unoscillated
        FD = (var_target - fitter.target)/fitter.FD_unoscillated
        if np.any(self.binEdges):
            binCenters = average_by_bin_edge(fitter.Ebins, fitter.Ebins, self.binEdges)
            ND = average_by_bin_edge(ND, fitter.Ebins, self.binEdges)
            FD = average_by_bin_edge(FD, fitter.Ebins, self.binEdges)
        else:
            binCenters = fitter.Ebins
        NDLine, = self.axUp.plot(binCenters,
                                 ND,
                                 color = DUNEblue)
        FDLine, = self.axUp.plot(binCenters,
                                 FD,
                                 color = DUNElightOrange)

        self.legLineList.append(NDLine)
        self.legLineList.append(FDLine)
        self.legLabelList.append("ND")
        self.legLabelList.append("FD")
        
        self.axUp.axvline(x = fitter.Ebounds[0],
                          ls = '--',
                          color = 'r')
        self.axUp.axvline(x = fitter.Ebounds[1],
                          ls = '--',
                          color = 'r')

        self.thisDiff = ND - FD
        self.axLo.plot(binCenters,
                       self.thisDiff,
                       color = DUNEdarkOrange)
        self.axLo.axvline(x = fitter.Ebounds[0],
                          ls = '--',
                          color = 'r')
        self.axLo.axvline(x = fitter.Ebounds[1],
                          ls = '--',
                          color = 'r')

        self.axUp.legend(self.legLineList,
                         self.legLabelList)
 
class mike_summary_plot (plot):
    def __init__(self, fitter, varType = "ppfx", title = None, makeIndividuals = False, ax = None, color = None, ylabels = True, binEdges = None, **kwargs):
        super(mike_summary_plot, self).__init__(**kwargs)

        if not ax:
            self.fig = plt.figure()
            self.ax = self.fig.gca()
            self.fig.tight_layout()
        else:
            self.ax = ax
            
        if title:
            self.ax.set_title(title)

        self.makeIndividuals = makeIndividuals

        self.ax.grid()
        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        if ylabels:
            self.ax.set_ylabel(r'$\left[\left(\Phi^{ND}_V - \Phi^{ND}_N\right) - \left(\Phi^{FD}_V - \Phi^FD_N\right)\right]/\Phi^{FD}_{\mathrm{unosc.}}$')
        self.ax.set_xlim(0, 10)
        self.ax.set_ylim(-0.12, 0.12)
        self.ax.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
        self.ax.set_yticklabels(["-10\%", "-5\%", "0\%", "5\%", "10\%"])

        if varType == "ppfx":
            varTargetsUnosc = fitter.FD_unosc_ppfx_univs
            varNDs = fitter.ND_ppfx_univs
            self.varKeys = range(100)
        elif varType == "other":
            varTargetsUnosc = fitter.FD_unosc_other_univs
            varNDs = fitter.ND_other_univs
            self.varKeys = systKeys

        if np.any(binEdges):
            binCenters = average_by_bin_edge(fitter.Ebins, fitter.Ebins, binEdges)
        else:
            binCenters = fitter.Ebins

        diffs = []
        self.subPlots = []
        for varKey in self.varKeys:
            thisMikePlot = mike_plot(fitter, varType, varKey, binEdges = binEdges)
            diffs.append(thisMikePlot.thisDiff)
            if not makeIndividuals:
                plt.close(thisMikePlot.fig)
            else:
                self.subPlots.append(thisMikePlot)
                
        self.ax.fill_between(binCenters,
                             *np.quantile(diffs, (0.16, 0.84), axis = 0),
                             color = color,
                             alpha = 0.5)
        self.line, = self.ax.plot(binCenters,
                                  np.median(diffs, axis = 0),
                                  color = color)
        self.ax.axvline(x = fitter.Ebounds[0],
                        ls = '--',
                        color = 'r')
        self.ax.axvline(x = fitter.Ebounds[1],
                        ls = '--',
                        color = 'r')

        plt.tight_layout()
    def save(self, fileName, remake = True, *args, **kwargs):
        plot.save(self, fileName, remake = remake, *args, **kwargs)
        if self.makeIndividuals:
            for subPlot, varKey in zip(self.subPlots, self.varKeys):
                subPlot.save(fileName.replace(".", "_"+str(varKey)+"."), remake = remake, **kwargs)
class L_curve_plot (plot):
    def __init__(self, fitter = None, **kwargs):
        super(L_curve_plot, self).__init__(**kwargs)

        self.fig = plt.figure()
        self.ax = self.fig.gca()

        self.legLabelList = []
        self.legLineList = []

        self.legArgs = {"frameon": False}
        
        if fitter:
            self.add(fitter, **kwargs)

        self.ax.loglog()
        self.ax.set_xlabel(r'Residual Norm = $|P(\hat{\Phi}_{ND} \vec{c} - \hat{\vec{\Phi}}_{FD})|$')
        # self.ax.set_ylabel(r'Solution Norm = $|A \vec{c}|$')
        self.ax.set_ylabel(r'Stat. Var. Norm = $|P \sigma_{ND}|$')

        self.fig.tight_layout()
        
    def add(self, fitter, regRange = np.logspace(-11, -5, 1000), highlight = [], **kwargs):
        if "ND" in kwargs: # user specified an alternative ND matrix
            ND = kwargs["ND"] 
        else:
            ND = [None]
        if "target" in kwargs: # user specified an alternative target array
            target = kwargs["target"]
        else:
            target = [None]

        res = []
        sol = []
        for reg in regRange:
            if not "HCreg" in kwargs:
                HCreg = reg
            else:
                HCreg = kwargs["HCreg"]

            fitter.calc_coeffs(reg, HCreg, ND = ND, target = target)
            res.append(fitter.residual_norm(**kwargs)) # residual norm (sum of squared residuals)
            sol.append(fitter.solution_norm(**kwargs)) # solution norm (UNWEIGHTED penalty matrix)

        line, = self.ax.plot(res, sol)

        if "label" in kwargs:
            self.legLabelList.append(kwargs["label"])
            self.legLineList.append(line)
            
        for reg in highlight: # highlight certain points with a scatter marker
            if not "HCreg" in kwargs:
                HCreg = reg
            else:
                HCreg = kwargs["HCreg"]

            fitter.calc_coeffs(reg, HCreg, ND = ND, target = target)
            point = self.ax.scatter(fitter.residual_norm(**kwargs),
                                    fitter.solution_norm(**kwargs))

            self.legLabelList.append(r'$\lambda_{OA} = $'+float_to_sci(reg))
            self.legLineList.append(point)

        if self.legLabelList:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

class L_curve_curvature_plot (plot):
    def __init__(self, fitter = None, **kwargs):
        super(L_curve_curvature_plot, self).__init__(**kwargs)

        self.fig = plt.figure()
        self.ax = self.fig.gca()

        self.legLabelList = []
        self.legLineList = []
        
        if fitter:
            self.opt = self.add(fitter, **kwargs)
            print (self.opt)
            
        self.ax.loglog()
        self.ax.set_xlabel(r'$\lambda$')
        self.ax.set_ylabel(r'Curvature')

        self.fig.tight_layout()
        
    def add(self, fitter, regRange = np.logspace(-11, -5, 1000), showOpt = True, highlight = [], **kwargs):
        if "ND" in kwargs:
            ND = kwargs["ND"]
        else:
            ND = [None]
        if "target" in kwargs:
            target = kwargs["target"]
        else:
            target = [None]

        res = []
        sol = []
        for reg in regRange:
            if not "HCreg" in kwargs:
                HCreg = reg
            else:
                HCreg = kwargs["HCreg"]

            fitter.calc_coeffs(reg, HCreg, ND = ND, target = target, **kwargs)
            res.append(fitter.residual_norm(**kwargs))
            # sol.append(fitter.solution_norm(**kwargs))
            sol.append(fitter.statvar_norm(**kwargs))

        fac = 1.
            
        res = np.array(res)
        sol = np.array(sol)
        dl = np.diff(regRange)
        xi = np.log(sol)
        rho = np.log(res)
        xi_prime = (1/fac)*np.diff(xi)/dl
        rho_prime = np.diff(rho)/dl
        xi_prime_prime = fac*np.diff(xi_prime)/dl[:-1]
        rho_prime_prime = np.diff(rho_prime)/dl[:-1]

        curv = 2*(rho_prime[:-1]*xi_prime_prime - rho_prime_prime*xi_prime[:-1])/np.power(np.power(rho_prime[:-1], 2) + np.power(xi_prime[:-1], 2), 3./2)

        line, = self.ax.plot(regRange[1:-1], curv)

        # maxCurv = np.max(np.abs(curv[~np.isnan(curv)]))
        # self.opt_l = regRange[1:-1][np.abs(curv) == maxCurv][0]

        maxCurv = np.max(curv[~np.isnan(curv)])
        self.opt_l = regRange[1:-1][curv == maxCurv][0]

        if showOpt:
            self.ax.axvline(x = self.opt_l,
                            label = r'$\lambda = $'+str(self.opt_l),
                            ls = '--',
                            color = 'red')

        if "label" in kwargs:
            self.legLabelList.append(kwargs["label"])
            self.legLineList.append(line)
            
        if self.legLabelList:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           frameon = False)

        
        return self.opt_l

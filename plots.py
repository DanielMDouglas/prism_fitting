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

class plot:
    def show(self, *args, **kwargs):
        self.fig.show(*args, **kwargs)

    def save(self, fileName, remake = True, *args, **kwargs):
        if os.path.exists and not remake:
            pass
        else:
            self.fig.savefig(fileName, *args, **kwargs)

class fit_and_ratio_plot (plot):
    def __init__(self, fitter = None, useTarget = True, useFit = True, title = None, Ebounds = True, **kwargs):
        
        bandBounds = (0.16, 0.84)

        self.fig = plt.figure(figsize = (8.5, 5))
        gs = GridSpec(2, 1,
                      figure = self.fig,
                      height_ratios = [0.7, 0.3],
                      hspace = 0)
        self.axUp = self.fig.add_subplot(gs[0, :])
        self.axLo = self.fig.add_subplot(gs[1, :])

        self.legLineList = []
        self.legLabelList = []

        self.legArgs = {"frameon": True,
                        "loc": "upper right"}
        if title:
            self.legArgs.update({"title": title})
            
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

        self.fig.tight_layout()

        self.full_flux = True

    def add_target(self, fitter, label = None, full_flux = True,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        targetNomLine, = self.axUp.plot(fitter.Ebins,
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
        NDNomLine, = self.axUp.plot(fitter.Ebins,
                                    np.dot(fitter.ND_full, fitter.c),
                                    color = color,
                                    **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(r'Fluxes up to ' + str(fitter.maxOA) + r'm')
        else:
            self.legLabelList.append(label)

        if self.full_flux:
            target = fitter.FD_oscillated
            # denom = fitter.FD_oscillated
            denom = fitter.FD_unoscillated
        else:
            target = self.target
            denom = self.target
        self.axLo.plot(fitter.Ebins,
                       (np.dot(fitter.ND_full, fitter.c) - target)/denom,
                       color = color)
        
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)

class fit_and_ratio_plot_with_sliders (plot):
    def __init__(self, fitter = None, useTarget = True, useFit = True, title = None, Ebounds = True, regLims = [-12, -3], **kwargs):
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

        self.fig.tight_layout()

        self.full_flux = True

    def add_target(self, fitter, label = None, full_flux = True,
                   Ebounds = True, color = 'black', **kwargs):
        self.full_flux = full_flux
        if full_flux:
            self.target = fitter.FD_oscillated
        else:
            self.target = fitter.target
        targetNomLine, = self.axUp.plot(fitter.Ebins,
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
        
    def add(self, fitter, label = None, color = None, **kwargs):
        self.fitters.append(fitter)
        
        NDNomLine, = self.axUp.plot(fitter.Ebins,
                                    np.dot(fitter.ND_full, fitter.c),
                                    color = color,
                                    **kwargs)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(r'Fluxes up to ' + str(fitter.maxOA) + r'm')
        else:
            self.legLabelList.append(label)

        if self.full_flux:
            target = fitter.FD_oscillated
            # denom = fitter.FD_oscillated
            denom = fitter.FD_unoscillated
        else:
            target = self.target
            denom = self.target
        RatioLine, = self.axLo.plot(fitter.Ebins,
                                    (np.dot(fitter.ND_full, fitter.c) - target)/denom,
                                    color = color)
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
            fitter.calc_coeffs(reg, reg)
            self.legLineList[i+1].set_ydata(np.dot(fitter.ND_full, fitter.c))
            self.ratioLineList[i].set_ydata((np.dot(fitter.ND_full, fitter.c) - target)/denom)
        # self.loBoundLineUp.set_xdata([loBound, loBound])
        # self.loBoundLineDown.set_xdata([loBound, loBound])
        # self.hiBoundLineUp.set_xdata([hiBound, hiBound])
        # self.hiBoundLineDown.set_xdata([hiBound, hiBound])
            
        self.fig.canvas.draw_idle()


class coeff_plot (plot):
    def __init__(self, fitter = None, title = "Coefficients", HC = True, legend = True, **kwargs):
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
            ylim = kwargs["ylim"]
        else:
            ylim = 2.5e-7

        self.OAax.set_ylim(-ylim, ylim)
        self.OAax.grid(True)
        self.OAax.set_xlabel(r'$D_{OA}$ [m]')
        self.OAax.set_ylabel(r'$c_i$')
        self.OAax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))

        if self.HC:
            self.HCax.set_ylim(-ylim, ylim)
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

        if fitter:
            self.add(fitter, **kwargs)

    def add(self, fitter, label = None, color = None, HCplot = False, **kwargs):
        # self.OAax.plot(fitter.OAbins[fitter.OAbins <= fitter.maxOA],
        #                fitter.cOA,
        #                color = color,
        #                label = label)
        self.OAax.step(fitter.OAbins[fitter.OAbins <= fitter.maxOA],
                       fitter.cOA,
                       color = color,
                       label = label,
                       where = 'mid')
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

        if self.legend:
            if label:
                if self.legendLoc == "left":
                    self.OAax.legend(**self.legArgs)
                elif self.legendLoc == "right":
                    self.HCax.legend(**self.legArgs)

class ND_flux_plot (plot):
    def __init__(self, fitter, title = "ND Flux"):
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
    def __init__(self, fitter, slices = [0, 10, 20, 30, 40, 50, 60]):
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
            meanLine, = self.ax.plot(fitter.Ebins,
                                     fitter.ND_full[:,sliceInd])
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
        
        self.ax.set_title("ND Flux")
        self.fig.tight_layout()

class FD_flux_plot (plot):
    def __init__(self, fitter = None, title = "FD Oscillated Flux", aspect = None, figSize = (6.4, 4.8), legendOn = True, **kwargs):
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

        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        self.ax.set_title(title)

        # self.ax.grid()
        self.ax.set_title(title)
        self.fig.tight_layout()

            
    def add(self, fitter, label = None, color = None, **kwargs):
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
            
        line = self.ax.plot(fitter.Ebins,
                            fitter.FD_oscillated,
                            color = color)
        
        self.legLineList.append(line[0])
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)
        
        self.ymax = max(self.ymax, 1.2*np.max(fitter.FD_oscillated))
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
            
        line = self.ax.plot(fitter.Ebins,
                            fitter.FD_unoscillated,
                            color = color)
        
        self.legLineList.append(line[0])
        
        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)
        
        self.ymax = max(self.ymax, 1.2*np.max(fitter.FD_unoscillated))
        self.ax.set_ylim(0, self.ymax)

        return line

                
    def add_fit(self, fitter, color = None, linestyle = None, label = r'ND Flux Match'):
        line, = self.ax.plot(fitter.Ebins,
                            np.dot(fitter.ND_full, fitter.c),
                            color = color,
                            linestyle = linestyle)
        self.legLineList.append(line)
        self.legLabelList.append(label)

        if self.legendOn:
            self.ax.legend(self.legLineList,
                           self.legLabelList,
                           **self.legArgs)

        return line
            
class FD_flux_osc_and_unosc_plot (plot):
    def __init__(self, fitter, title = "FD Flux", inset_text = None, logScale = False, figSize = (6.4, 4.8)):
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
                                np.dot(fitter.ND_full, fitter.c),
                                color = color)
        self.legLineList.append(fitLine)
        self.legLabelList.append(r'ND Flux Match')
        
        self.ax.legend(self.legLineList,
                       self.legLabelList,
                       frameon = False)
 
class mike_plot (plot):
    def __init__(self, fitter = None, varType = "ppfx", varKey = 0, label = r'ND \& FD', title = None, binEdges = [], ylim = None):
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

        ND = (np.dot(ND_univ, fitter.c) - np.dot(fitter.ND, fitter.c))/fitter.FD_unoscillated
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
    def __init__(self, fitter, varType = "ppfx", title = None, makeIndividuals = False, ax = None, color = None, ylabels = True, binEdges = None):
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
        self.fig = plt.figure()
        self.ax = self.fig.gca()

        self.legLabelList = []
        self.legLineList = []

        self.legArgs = {"frameon": False}
        
        if fitter:
            self.add(fitter, **kwargs)

        self.ax.loglog()
        self.ax.set_xlabel(r'Residual Norm = $|P(\hat{\Phi}_{ND} \vec{c} - \hat{\vec{\Phi}}_{FD})|$')
        self.ax.set_ylabel(r'Solution Norm = $|A \vec{c}|$')

        self.fig.tight_layout()
        
    def add(self, fitter, regRange = np.logspace(-11, -5, 1000), highlight = [], **kwargs):
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

            fitter.calc_coeffs(reg, HCreg, ND = ND, target = target)
            res.append(fitter.residual_norm(**kwargs))
            sol.append(fitter.solution_norm(**kwargs))

        line, = self.ax.plot(res, sol)

        if "label" in kwargs:
            self.legLabelList.append(kwargs["label"])
            self.legLineList.append(line)
            
        for reg in highlight:
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
        self.fig = plt.figure()
        self.ax = self.fig.gca()

        self.legLabelList = []
        self.legLineList = []
        
        if fitter:
            self.add(fitter, **kwargs)

        self.ax.semilogx()
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
            sol.append(fitter.solution_norm(**kwargs))

        res = np.array(res)
        sol = np.array(sol)
        dl = np.diff(regRange)
        xi = np.log(sol)
        rho = np.log(res)
        xi_prime = np.diff(xi)/dl
        rho_prime = np.diff(rho)/dl
        xi_prime_prime = np.diff(xi_prime)/dl[:-1]
        rho_prime_prime = np.diff(rho_prime)/dl[:-1]

        curv = 2*(rho_prime[:-1]*xi_prime_prime - rho_prime_prime*xi_prime[:-1])/np.power(np.power(rho_prime[:-1], 2) + np.power(xi_prime[:-1], 2), 3./2)

        line, = self.ax.plot(regRange[1:-1], curv)

        maxCurv = np.max(curv[~np.isnan(curv)])
        self.opt_l = regRange[1:-1][curv == maxCurv]
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

        
        return self.opt_l[0]

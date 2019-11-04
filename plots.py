import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from utils import *

DUNEblue = '#7FAED5'
DUNElightOrange = '#F19E54'
DUNEdarkOrange = '#F0652B'
DUNEgreen = '#8ACA6F'
DUNEgray = '#626466'
DUNEyellow = '#FBD03F'

DUNEcolors = [DUNEblue,
              DUNElightOrange,
              DUNEdarkOrange,
              DUNEgreen,
              DUNEgray,
              DUNEyellow]

matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
matplotlib.rc('text', usetex = True)
matplotlib.rc('axes', prop_cycle = matplotlib.cycler(color = DUNEcolors))

class fit_and_ratio_plot:
    def __init__(self, fitter, title = None, **kwargs):
        
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
        
        self.add(fitter, plotTarget = True, **kwargs)

        self.axUp.set_xlim(0, 10)
        self.axUp.set_xticklabels([])
        self.axUp.grid(True, which = 'both')
        self.axUp.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]', labelpad = 15)

        self.axLo.grid(True)
        self.axLo.set_xlim(0, 10)
        self.axLo.set_ylim(-0.12, 0.12)
        self.axLo.set_xlabel(r'$E_\nu$ [GeV]')
        self.axLo.set_ylabel(r'$\frac{ND - FD (osc.)}{FD (unosc.)}$', labelpad = 5)
        self.axLo.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
        self.axLo.set_yticklabels(["-10\%", "-5\%", "0\%", "5\%", "10\%"])
        self.axLo.grid(True, which = 'both')

        self.fig.tight_layout()

    def add(self, fitter, label = None, color = None, plotTarget = False):
        if plotTarget:
            targetNomLine, = self.axUp.plot(fitter.Ebins,
                                            fitter.FD_oscillated,
                                            color = 'black')
            self.legLineList.append(targetNomLine)
            self.legLabelList.append(r'Target Flux')

            
            self.axUp.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.axUp.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')
            self.axUp.arrow(fitter.Ebounds[0], 3.4e-15,
                            0.15, 0,
                            width = 2.e-17,
                            head_length = 0.05,
                            color = 'red')
            self.axUp.arrow(fitter.Ebounds[1], 1.5e-15,
                            -0.15, 0,
                            width = 2.e-17,
                            head_length = 0.05,
                            color = 'red')
        
            self.axLo.axvline(x = fitter.Ebounds[0],
                              ls = '--',
                              color = 'red')
            self.axLo.axvline(x = fitter.Ebounds[1],
                              ls = '--',
                              color = 'red')

        NDNomLine, = self.axUp.plot(fitter.Ebins,
                                    np.dot(fitter.ND_full, fitter.c),
                                    color = color)
        self.legLineList.append(NDNomLine)
        
        if not label:
            self.legLabelList.append(r'Fluxes up to ' + str(fitter.maxOA) + r'm')
        else:
            self.legLabelList.append(label)

        self.axLo.plot(fitter.Ebins,
                       (np.dot(fitter.ND_full, fitter.c) - fitter.FD_oscillated)/fitter.FD_unoscillated,
                       color = color)

    def show(self):
        self.axUp.legend(self.legLineList,
                         self.legLabelList,
                         **self.legArgs)
        self.fig.show()

    def save(self, fileName):
        self.fig.savefig(fileName)

class coeff_plot:
    def __init__(self, fitter, title = "Coefficients", HC = True, **kwargs):
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

        self.add(fitter, **kwargs)
            
        ylim = 2.5e-7

        self.OAax.set_ylim(-ylim, ylim)
        self.OAax.grid(True)
        self.OAax.set_xlabel(r'$D_{OA}$ [m]')
        self.OAax.set_ylabel(r'$c_i$')

        self.HCax.set_ylim(-ylim, ylim)
        self.HCax.grid(True)
        self.HCax.set_xlabel(r'Horn Current [kA]')
        self.HCax.set_yticklabels([])

        self.fig.suptitle(title)
        self.fig.tight_layout()

    def add(self, fitter, label = None, color = None, HCplot = False):
        self.OAax.plot(fitter.OAbins[fitter.OAbins <= fitter.maxOA],
                       fitter.cOA,
                       color = color,
                       label = label)
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

    def show(self):
        self.OAax.legend()
        self.fig.show()
        
    def save(self, fileName):
        self.fig.savefig(fileName)

class ND_flux_plot:
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

    def show(self):
        self.fig.show()
        
    def save(self, fileName):
        self.fig.savefig(fileName)

class ND_flux_slice_plot:
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
            print meanLines, labels
            self.ax.legend(meanLines,
                           labels,
                           frameon = False)
            
        self.ax.set_xlim(np.min(fitter.Ebins), np.max(fitter.Ebins))
        
        self.ax.set_title("ND Flux")
        self.fig.tight_layout()

    def show(self):
        self.fig.show()
        
    def save(self, fileName):
        self.fig.savefig(fileName)

class FD_flux_plot:
    def __init__(self, fitter, color = None):
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        if "FD_unosc_univs" in dir(fitter):
            oscillated = fitter.FD_unosc_univs*self.Posc
            lower, upper = np.quantile(oscillated, [0.16, 0.84], axis = 0)
            self.ax.fill_between(fitter.Ebins,
                                 lower,
                                 upper,
                                 alpha = 0.5,
                                 color = color)
            
        self.ax.plot(fitter.Ebins,
                     fitter.FD_oscillated,
                     color = '#ff7f0e')

        self.ax.set_xlabel(r'$E_\nu$ [GeV]')
        self.ax.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        self.ax.grid()
        self.ax.set_title("FD Oscillated Flux")
        self.fig.tight_layout()

    def show(self):
        self.fig.show()

    def save(self, fileName):
        self.fig.savefig(fileName)

class FD_flux_osc_and_unosc_plot:
    def __init__(self, fitter, title = "FD Flux", inset_text = None, logScale = False):
        self.fig = plt.figure()
        self.ax = plt.gca()

        self.legLabelList = [r'Oscillated',
                             r'Unoscillated']

        oscLine, = self.ax.plot(fitter.Ebins,
                                fitter.FD_oscillated,
                                color = DUNEblue)

        unoscLine, = self.ax.plot(fitter.Ebins,
                                  fitter.FD_unoscillated,
                                  color = DUNElightOrange)
        
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

        print self.legLineList, self.legLabelList
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

    def show(self):
        self.fig.show()
        
    def save(self, fileName):
        self.fig.savefig(fileName)

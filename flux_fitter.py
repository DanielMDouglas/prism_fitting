from utils import *
from fluxes import *
from oscProbs import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class flux_fitter:
    def __init__(self, beamMode, FDfromFlavor, FDtoFlavor, NDflavor, oscParam = None):
        self.FD_unoscillated = FD_nominal[beamMode][FDfromFlavor].load()
        self.ND_full = ND_nominal[beamMode][NDflavor].load()
        self.ND = self.ND_full
        self.Ebins = Ebins
        self.OAbins = OAbins

        self.Ebounds = (0, np.max(self.Ebins))
        self.OutOfRegionFactors = (0, 0)
        self.maxOA = np.max(self.OAbins)
        
        if not oscParam:
            self.Posc = oscProb(FDfromFlavor, FDtoFlavor).load(self.Ebins)
        else:
            self.Posc = oscParam.load(self.Ebins)
            
        self.FD_oscillated = self.FD_unoscillated*self.Posc
        # add target shaping later
        # maybe don't need it, since target is weight-able?
        self.target = self.FD_oscillated

    def set_fit_region(self, energies = None, peaks = None):
        if energies:
            if not len(energies) == 2:
                print "Energy bounds must be an iterable of length 2!"
                return
            self.Ebounds = energies
        elif peaks:
            if not len(peaks) == 2:
                print "Peak indices must be an iterable of length 2!"
                return
            else:
                # this should be an odd integer
                # otherwise I have to think harder about how to do this...
                window_size = 5
                fringe_size = window_size/2
                smoothed = np.convolve(self.FD_oscillated,
                                       (1/float(window_size))*np.ones(window_size))[fringe_size:-fringe_size]
                threshold = 0.1*np.max(smoothed)
                foundPeaks = self.Ebins[(np.diff(smoothed, n = 1, prepend = 0) > 0) &
                                        (np.diff(smoothed, n = 1, append = 0) < 0) &
                                        (smoothed > threshold)]

                # peaks should be specified from the right,
                # i.e (1, 3) fits between the first and third-highest energy peaks
                self.Ebounds = (foundPeaks[-peaks[1]], foundPeaks[-peaks[0]])                
        
    def set_OOR(self, weights):
        self.OutOfRegionFactors = weights

    def set_oscHypothesis(self, Posc):
        self.Posc = Posc
        self.FD_oscillated = self.FD_unoscillated*Posc
        self.target = self.FD_oscillated

    def set_maxOA(self, maxOA):
        self.maxOA = maxOA
        self.ND = self.ND_full.T[OAbins <= maxOA].T

    def calc_coeffs(self, reg):
        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        Gamma = reg*A

        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))

        # ND matrix
        NDmatr = np.matrix(self.ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T, Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), self.target)

        self.c = np.array(np.dot(RHS, LHS.I)).squeeze()

    def set_coeffs(self, other):
        self.coeffs = other.coeffs


    def plot_target_fit_and_ratio(self, outfileName = None, title = None, resid = False):
        matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
        matplotlib.rc('text', usetex = True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.close('all')
        fig = plt.figure(figsize = (8.5, 5))
        gs = GridSpec(2, 1,
                      figure=fig,
                      height_ratios = [0.7, 0.3],
                      hspace = 0)
        axUp = fig.add_subplot(gs[0, :])
        axLo = fig.add_subplot(gs[1, :], sharex=axUp)

        axUp.plot(self.Ebins,
                  self.FD_oscillated,
                  color = 'k')
        axUp.plot(self.Ebins,
                  np.dot(self.ND, self.c),
                  label = r'Fluxes up to '+str(self.maxOA)+r'm',
                  color = '#ff7f0e')
        axUp.axvline(x = self.Ebounds[0],
                     ls = '--',
                     color = '#1f77b4',
                     label = r'Fit region')
        axUp.axvline(x = self.Ebounds[1],
                     ls = '--',
                     color = '#1f77b4')
        axUp.arrow(self.Ebounds[0], 3.4e-15,
                   0.15, 0,
                   width = 2.e-17,
                   head_length = 0.05,
                   color = '#1f77b4')
        axUp.arrow(self.Ebounds[1], 1.5e-15,
                   -0.15, 0,
                   width = 2.e-17,
                   head_length = 0.05,
                   color = '#1f77b4')
        axUp.legend(title = title)
        axUp.grid()
        axUp.set_xlim(0, 5)
        # axUp.set_ylim(-3.e-16, 4e-15)
        axUp.set_ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]', labelpad = 15)

        axLo.plot(self.Ebins,
                  (np.dot(self.ND, self.c) - self.FD_oscillated)/self.FD_unoscillated,
                  color = '#ff7f0e')
        axLo.axvline(x = self.Ebounds[0],
                     ls = '--',
                     color = '#1f77b4')
        axLo.axvline(x = self.Ebounds[1],
                     ls = '--',
                     color = '#1f77b4')
        axLo.grid()
        axLo.set_ylim(-0.5, 0.5)
        axLo.set_xlabel(r'$E_\nu$ [GeV]')
        axLo.set_ylabel(r'$\frac{ND - FD (osc.)}{FD (unosc.)}$', labelpad = 5)

        if resid:
            resid = np.sum(np.power(np.dot(self.ND, self.c) - self.FD_oscillated, 2))/len(self.Ebins)
            frac = str(round(float(str(resid).split('e')[0]), 2)) 
            exp = str(int(str(resid).split('e')[1]))
            axLo.text(self.Ebounds[0]+0.3, 0.2, r'Sum of Squared Residuals per E bin: $' + frac + r'\times 10^{' + exp + r'}$')
        
        plt.tight_layout()

        if outfileName:
            plt.savefig(outfileName)
            plt.clf()
        else:
            plt.show()

    def plot_ND_flux(self, outfileName = None):
        matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
        matplotlib.rc('text', usetex = True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.close('all')
        fig = plt.figure()
        plt.pcolormesh(self.Ebins, self.OAbins, self.ND_full.T)

        plt.xlabel(r'$E_\nu$ [GeV]')
        plt.ylabel(r'$D_{OA}$ [m]')
        
        cb = plt.colorbar()
        cb.set_label(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        plt.title("ND Flux")
        plt.tight_layout()

        if outfileName:
            plt.savefig(outfileName)
            plt.clf()
        else:
            plt.show()

    def plot_FD_flux(self, outfileName = None):
        matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
        matplotlib.rc('text', usetex = True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.close('all')
        fig = plt.figure()
        plt.plot(self.Ebins, self.FD_oscillated)

        plt.xlabel(r'$E_\nu$ [GeV]')
        plt.ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        plt.grid()
        plt.title("FD Oscillated Flux")
        plt.tight_layout()

        if outfileName:
            plt.savefig(outfileName)
            plt.clf()
        else:
            plt.show()

    def plot_FD_flux_osc_and_unosc(self, outfileName = None, title = None, inset_text = None):
        matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
        matplotlib.rc('text', usetex = True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.close('all')
        fig = plt.figure()
        plt.plot(self.Ebins, self.FD_oscillated, color = '#7eadd4', label = r'Oscillated')
        plt.plot(self.Ebins, self.FD_unoscillated, color = '#ef6530', label = r'Unoscillated')

        plt.xlabel(r'$E_\nu$ [GeV]')
        plt.ylabel(r'$\Phi$ [cm$^{-2}$ per POT per GeV]')

        # plt.legend(frameon = False)
        
        plt.xlim(np.min(self.Ebins), np.max(self.Ebins))
        plt.ylim(0, 1.2*np.max(self.FD_unoscillated))
        
        # plt.grid()
        logScale = False
        if logScale:
            plt.semilogx()
            
        if title:
            plt.title(title)
        if inset_text:
            top_line_y_loc = 1.05*np.max(self.FD_unoscillated)
            if logScale:
                top_line_x_loc = 2.e-2
            else:
                top_line_x_loc = 6.
                
            if type(inset_text) == list:
                line_spacing = 0.15*top_line_y_loc
                for i, line in enumerate(inset_text):
                    plt.text(top_line_x_loc,
                             top_line_y_loc - i*line_spacing,
                             line)

            else:
                plt.text(top_line_x_loc,
                         top_line_y_loc,
                         inset_text)
        plt.tight_layout()

        if outfileName:
            plt.savefig(outfileName)
            plt.clf()
        else:
            plt.show()

    def plot_coeffs(self, outfileName = None, title = None):
        matplotlib.rc('font', family = 'FreeSerif', size = 16, weight = 'bold')
        matplotlib.rc('text', usetex = True)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        plt.close('all')
        fig = plt.figure()
        plt.plot(self.OAbins[self.OAbins <= self.maxOA],
                 self.c)

        plt.grid()
        plt.title("Coefficients")
        plt.xlabel(r'$D_{OA}$ [m]')
        plt.ylabel(r'$c_i$')
        plt.ylim(-1.e-7, 1.e-7)
        
        if title:
            plt.title(title)
        plt.tight_layout()

        if outfileName:
            plt.savefig(outfileName)
            plt.clf()
        else:
            plt.show()

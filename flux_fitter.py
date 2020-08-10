from utils import *
from fluxes import *
from xSec import *
from oscProbs import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class flux_fitter:
    def __init__(self, beamMode, FDfromFlavor, FDtoFlavor, NDflavor, oscParam = None, Erebin = 1, OArebin = 1, useHC = True):
        """
        flux_fitter object which holds fluxes, settings, and solutions to 
        DUNE-PRISM style flux matching problems
        """
        # initial binning.  This is determined by the binning of the flux sim
        self.initNEbins = Ebins.size
        self.initNOAbins = OAbins.size

        self.Ebins = Ebins
        self.OAbins = OAbins

        self.EbinEdges = EbinEdges
        self.OAbinEdges = OAbinEdges
        
        # if we're using alternative horn currents, establish another set of bins for those
        if useHC:
            # currents is defined in fluxes.py each value should have an accompanying flux object
            self.HCbins = currents
            self.maxHC = np.max(currents)
        else:
            self.HCbins = []
            self.maxHC = 293. # 293 kA is the nominal horn current setting

        # initial maximum off-axis distance.  Determined by flux sim, usually 40m
        self.maxOA = np.max(self.OAbins)

        # Are we rebinning in energy?
        # if so, redefine the energy bins
        self.Erebin = Erebin
        if type(self.Erebin) == np.ndarray:
            self.Ebins = average_by_bin_edge(self.Ebins, Ebins, self.Erebin)
            self.EbinEdges = self.Erebin
        elif type(self.Erebin) == int:
            self.Ebins = average(self.Ebins, self.Erebin)
            self.EbinEdges = self.EbinEdges[::self.Erebin]

        # likewise with off-axis position
        self.OArebin = OArebin
        if type(self.OArebin) == np.ndarray:
            self.OAbins = average_by_bin_edge(self.OAbins, OAbins, self.OArebin)
            self.OAbinEdges = self.OArebin
        elif type(self.OArebin) == int:
            self.OAbins = average(self.OAbins, self.OArebin)
            self.OAbinEdges = self.OAbinEdges[::self.OArebin]

        # set flavor information
        self.beamMode = beamMode
        self.FDfromFlavor = FDfromFlavor
        self.FDtoFlavor = FDtoFlavor
        self.NDflavor = NDflavor

        # if not specified, load the default oscillation profile
        if not oscParam:
            self.Posc = oscProb(FDfromFlavor, FDtoFlavor).load(self.Ebins)
        else:
            self.Posc = oscParam.load(self.Ebins)
            
        # load cross sections
        self.load_xSec()

        # load nominal fluxes
        self.load_nom()

        # set default bounds for the fit
        self.Ebounds = (0, np.max(self.Ebins))
        self.OutOfRegionFactors = (0, 0)

        # Don't load systematics by default
        self.ppfx_systs_loaded = False
        self.other_systs_loaded = False

    def load_xSec(self):
        """
        Load cross sections for the set flavors
        """
        self.NDxSec = xSec[self.NDflavor]["total"].load(binEdges = [self.EbinEdges])
        self.FDxSec = xSec[self.FDtoFlavor]["total"].load(binEdges = [self.EbinEdges])
        
    def load_FD_nom(self):
        """
        Load the nominal far detector flux and update the FD_oscillated and target 
        """
        flux = FD_nominal[self.beamMode][self.FDfromFlavor]
        self.FD_unoscillated = flux.load(binEdges = [self.EbinEdges])

        self.FD_oscillated = self.FD_unoscillated*self.Posc
        self.target = self.FD_oscillated

        self.FD_rate = FDscale*self.FD_oscillated*self.FDxSec
        self.FD_rate_statErr = np.sqrt(self.FD_rate)

    def load_ND_OA_nom(self):
        """
        Load the nominal near detector flux for off-axis positions at 293 kA
        """
        flux = ND_nominal[self.beamMode][self.NDflavor]
        self.ND_OA = flux.load(binEdges = [self.EbinEdges, self.OAbinEdges])

        self.ND_OA_rate = NDscale*(self.ND_OA.T*self.NDxSec).T
        self.ND_OA_rate_statErr = np.sqrt(self.ND_OA_rate)
        
        if "ND_HC" in dir(self):
            self.ND = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
            self.ND_rate = np.concatenate((self.ND_OA_rate, self.ND_HC_rate), axis = 1)
            self.ND_rate_statErr = np.concatenate((self.ND_OA_rate_statErr, self.ND_HC_rate_statErr), axis = 1)
            
    def load_ND_HC_nom(self):
        """
        Load the nominal near detector flux for all alternative horn currents
        """
        if list(self.HCbins):
            self.ND_HC = 12*np.array([ND_HC_shifts[self.beamMode][self.NDflavor][current].load(binEdges = [self.EbinEdges])
                                      for current in self.HCbins]).T
        else:
            self.ND_HC = np.ndarray((len(self.Ebins), 0))
            
        self.ND_HC_rate = NDscale*(self.ND_HC.T*self.NDxSec).T
        self.ND_HC_rate_statErr = np.sqrt(self.ND_HC_rate)

        if "ND_OA" in dir(self):
            self.ND = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
            self.ND_rate = np.concatenate((self.ND_OA_rate, self.ND_HC_rate), axis = 1)
            self.ND_rate_statErr = np.concatenate((self.ND_OA_rate_statErr, self.ND_HC_rate_statErr), axis = 1)
            
    def load_nom(self):
        """
        Load all nominal fluxes
        """
        self.load_FD_nom()
        self.load_ND_OA_nom()
        self.load_ND_HC_nom()
        
    def load_FD_ppfx_systs(self, nUniv = 100):
        """
        Load the hadron production throws for the far detector flux
        """
        flux = FD_ppfx_CV[self.beamMode][self.FDfromFlavor]
        FD_CV = flux.load(binEdges = [self.EbinEdges])
        FD = np.array([shift.load(binEdges = [self.EbinEdges]) for shift
                       in FD_ppfx_shifts[self.beamMode][self.FDfromFlavor][:nUniv]])
        
        FD /= FD_CV            
        FD *= self.FD_unoscillated

        self.FD_unosc_ppfx_univs = FD

    def load_ND_OA_ppfx_systs(self, nUniv = 100):
        """
        Load the hadron production throws for the near detector flux for off-axis positions
        """
        ND_OA_CV = ND_ppfx_CV[self.beamMode][self.FDfromFlavor].load(binEdges = [self.EbinEdges, self.OAbinEdges])
        ND_OA = np.array([shift.load(binEdges = [self.EbinEdges, self.OAbinEdges]) for shift
                          in ND_ppfx_shifts[self.beamMode][self.NDflavor][:nUniv]])

        ND_OA /= ND_OA_CV
        ND_OA *= self.ND_OA

        self.ND_OA_ppfx_univs = ND_OA
        if "ND_HC_ppfx_univs" in dir(self):
            self.ND_ppfx_univs = np.concatenate((self.ND_OA_ppfx_univs,
                                                 self.ND_HC_ppfx_univs),
                                                axis = 2)

    def load_ND_HC_ppfx_systs(self, nUniv = 100):
        """
        Load the hadron production throws for the near detector flux for alternative horn currents
        """
        ND_HC_CV = np.array([ND_HC_ppfx_CV[self.beamMode][self.FDfromFlavor][current].load(binEdges = [self.EbinEdges])
                             for current in self.HCbins]).T
        if np.any(self.HCbins):
            ND_HC = np.array([np.array([ND_HC_ppfx_shifts[self.beamMode][self.NDflavor][current][univ].load(binEdges = [self.EbinEdges])
                                        for univ in range(nUniv)]).T
                              for current in self.HCbins]).T
        else:
            ND_HC = np.ndarray((nUniv, self.initNEbins, len(self.HCbins)))

        ND_HC /= ND_HC_CV
        ND_HC *= 12*self.ND_HC
        
        self.ND_HC_ppfx_univs = ND_HC
        if "ND_OA_ppfx_univs" in dir(self):
            self.ND_ppfx_univs = np.concatenate((self.ND_OA_ppfx_univs,
                                                 self.ND_HC_ppfx_univs),
                                                axis = 2)

    def load_ppfx_systs(self, nUniv = 100):
        """
        Load all hadron production throws
        """
        self.load_FD_ppfx_systs(nUniv = nUniv)
        self.load_ND_OA_ppfx_systs(nUniv = nUniv)
        self.load_ND_HC_ppfx_systs(nUniv = nUniv)

        self.ppfx_systs_loaded = True

    def load_FD_other_systs(self):
        """
        Load systematic throw fluxes for specific beam parameters for far detector
        """
        FD = {key: shift.load(binEdges = [self.EbinEdges])
              for key, shift in FD_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.FD_unosc_other_univs = FD

    def load_ND_OA_other_systs(self):
        """
        Load systematic throw fluxes for specific beam parameters for near detector off-axis positions
        """
        ND_OA = {key: shift.load(binEdges = [self.EbinEdges])
                 for key, shift in ND_other_shifts[self.beamMode][self.NDflavor].items()}
        self.ND_OA_other_univs = ND_OA

        if "ND_HC_other_univs" in dir(self):
            self.ND_other_univs = {key: np.concatenate((self.ND_OA_other_univs[key],
                                                        self.ND_HC_other_univs[key]),
                                                       axis = 1)
                                   for key in systKeys}

    def load_ND_HC_other_systs(self):
        """
        Load systematic throw fluxes for specific beam parameters for near detector alternative horn currents
        """
        # Don't have HC beam systs yet, so fake them by assuming
        # the same fractional error as seen in 293 kA
        if np.any(self.HCbins):
            ND_HC = {}
            for key in HCsystKeys:
                currentList = []
                for current in self.HCbins:
                    flux = ND_HC_other_shifts[self.beamMode][self.NDflavor][current][key]
                    currentList.append(12*flux.load(binEdges=[self.EbinEdges]))
                arr = np.array(currentList).T
                ND_HC.update({key: arr})
        else:
            ND_HC = {key: np.ndarray((self.Ebins.size, len(self.HCbins)))
                     for key in HCsystKeys}
            
        self.ND_HC_other_univs = ND_HC
        if "ND_OA_other_univs" in dir(self):
            self.ND_other_univs = {key: np.concatenate((self.ND_OA_other_univs[key],
                                                        self.ND_HC_other_univs[HCkey]),
                                                       axis = 1)
                                   for key, HCkey in zip(systKeys, HCsystKeys)}

    def load_other_systs(self):
        """
        Load all systematic throw fluxes for specific beam parameters
        """
        self.load_FD_other_systs()
        self.load_ND_OA_other_systs()
        self.load_ND_HC_other_systs()

        self.other_systs_loaded = True
        
    def set_fit_region(self, energies = None, peaks = None):
        """
        Specify the energy region in which to fit the ND and target fluxes. The region can be specified by the keyword argument 'energies', or 'peaks', which are counted from high to low energy.  The arguments should be specified as an iterable of length 2
        """
        if energies:
            if not len(energies) == 2:
                print("Energy bounds must be an iterable of length 2!")
                return
            self.Ebounds = energies
        elif peaks:
            # WARNING: this is very touchy and tuned specifically for FD numu -> numu fluxes
            if not len(peaks) == 2:
                print("Peak indices must be an iterable of length 2!")
                return
            else:
                # this should be an odd integer
                # otherwise I have to think harder about how to do this...
                window_size = 7
                fringe_size = window_size/2
                smoothed = np.convolve(self.target,
                                       (1/float(window_size))*np.ones(window_size))[:-window_size+1]
                threshold = 0.05*np.max(smoothed)
                foundPeaks = self.Ebins[(np.diff(smoothed, n = 1, prepend = 0) > 0) &
                                        (np.diff(smoothed, n = 1, append = 0) < 0) &
                                        (smoothed > threshold)]

                # peaks should be specified from the right,
                # i.e (1, 3) fits between the first and third-highest energy peaks
                if peaks[0] == 0:
                    self.Ebounds = (foundPeaks[-peaks[1]], max(self.Ebins))
                else:
                    self.Ebounds = (foundPeaks[-peaks[1]], foundPeaks[-peaks[0]])
                print("found peaks at ", self.Ebounds)

    def set_OOR(self, weights):
        """
        Specify the relative weight of the regions outside of the "fit region" when contributing to the residual. The weights should be specified as an iterable of length 2
        """
        self.OutOfRegionFactors = weights

    def set_oscHypothesis(self, Posc):
        """
        
        """
        self.Posc = Posc.load(self.Ebins)
        self.FD_oscillated = self.FD_unoscillated*self.Posc
        self.target = self.FD_oscillated

    def set_maxOA(self, maxOA):
        self.maxOA = maxOA
        self.OAbins = self.OAbins[self.OAbins <= maxOA]
        self.OAbinEdges = np.concatenate([self.OAbinEdges[self.OAbinEdges < maxOA],
                                          [maxOA]])

        self.load_ND_OA_nom()
        
        if self.ppfx_systs_loaded:
            self.load_ND_OA_ppfx_systs()
        if self.other_systs_loaded:
            self.load_ND_OA_other_systs()
            
    def set_maxHC(self, maxHC):
        self.maxHC = maxHC
        self.HCbins = self.HCbins[self.HCbins <= maxHC]

        self.load_ND_HC_nom()
       
        if self.ppfx_systs_loaded:
            self.load_ND_HC_ppfx_systs()
        if self.other_systs_loaded:
            self.load_ND_HC_other_systs()

    def use_currents(self, theseCurrents):
        self.maxHC = max(list(theseCurrents) + [293])
        self.HCbins = np.array(theseCurrents)
        
        self.load_ND_HC_nom()
        
        if self.ppfx_systs_loaded:
            self.load_ND_HC_ppfx_systs()
        if self.other_systs_loaded:
            self.load_ND_HC_other_systs()
            
    def calc_coeffs(self, OAreg, HCreg, ND = [None], target = [None], fluxTimesE = False, NDtoFD = NDtoFD, **kwargs):
        if not np.any(ND):
            ND = self.ND
        if not np.any(target):
            target = self.target
        if fluxTimesE:
            ND = ND*self.Ebins.reshape(ND.shape[0], 1)
            target = target*self.Ebins
        
        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        nBinsHC = self.ND_HC.shape[1]
        
        # OA penalty term is a difference matrix
        OA_penalty = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        # # OA penalty term is the L1 norm
        # OA_penalty = np.eye(nBinsOA)
        HC_penalty = np.eye(nBinsHC)
        self.A = block_diag(OA_penalty, HC_penalty)
        Gamma = block_diag(OAreg*OA_penalty, HCreg*HC_penalty)
        self.Gamma = Gamma

        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(self.Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))
        self.P = P
        # ND matrix
        NDmatr = np.matrix(ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T, Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), target)
        
        self.c = np.array(np.dot(RHS, LHS.I)).squeeze()
        self.cOA = self.c[:nBinsOA]
        self.cHC = self.c[nBinsOA:]

        self.fluxPred = np.dot(self.ND, self.c)
        self.ratePred = NDtoFD*np.dot(self.ND_rate, self.c)
        self.ratePred_statErr = NDtoFD*np.sqrt(np.dot(self.ND_rate, np.power(self.c, 2)))
        
    def calc_coeffs_DFT(self, OAreg, HCreg, filt, ND = [None], target = [None], fluxTimesE = False):
        if not np.any(ND):
            ND = self.ND
        if not np.any(target):
            target = self.target
        if fluxTimesE:
            ND = ND*self.Ebins.reshape(ND.shape[0], 1)
            target = target*self.Ebins
        # if not np.any(filt):
        #     filt = np.ones_like(self.OAbins <= self.maxOA, dtype = float)
            
        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        nBinsHC = self.ND_HC.shape[1]
        
        # # OA penalty term is a difference matrix
        # OA_penalty = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        # # OA penalty term is the L1 norm
        # OA_penalty = np.eye(nBinsOA)
        # OA penalty is a frequency filter based on the DFT
        from scipy.linalg import dft
        # OA_penalty = np.diag(filt)*dft(nBinsOA) + OAreg*(np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1))
        OA_penalty = np.diag(filt)*dft(nBinsOA) + OAreg*(np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1))
        HC_penalty = np.eye(nBinsHC)
        self.A = block_diag(OA_penalty, HC_penalty)
        # Gamma = block_diag(OAreg*OA_penalty, HCreg*HC_penalty)
        Gamma = block_diag(OA_penalty, HCreg*HC_penalty)
        self.Gamma = Gamma

        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(self.Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))
        self.P = P

        # ND matrix
        NDmatr = np.matrix(ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T.conj(), Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), target)

        self.c = np.array(np.dot(RHS, LHS.I)).squeeze()
        self.cOA = self.c[:nBinsOA]
        self.cHC = self.c[nBinsOA:]

        self.fluxPred = np.dot(self.ND, self.c)
        self.ratePred = np.dot(self.ND_rate, self.c)
        self.ratePred_statErr = np.sqrt(np.dot(self.ND_rate, np.power(self.c, 2)))

    def calc_var_coeff_correction(self, reg, fluxTimesE = False):
        """
        Calculate a correction to the nominal coefficients to minimize the variance of the residual
        over different systematic universes
        """
        ND = np.concatenate((self.ND,) + tuple(ND for ND in self.ND_ppfx_univs))
        target = np.concatenate((self.target,) + tuple(FD*self.Posc for FD in self.FD_unosc_ppfx_univs))
        target -= np.dot(ND, self.c)
        # if fluxTimesE:
        #     ND = ND*self.Ebins.reshape(ND.shape[0], 1)
        #     target = target*self.Ebins
        
        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        nBinsHC = self.ND_HC.shape[1]
        
        # # OA penalty term is a difference matrix
        # OA_penalty = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        # OA penalty term is the L1 norm
        OA_penalty = np.eye(nBinsOA)
        HC_penalty = np.eye(nBinsHC)
        OAreg = reg
        HCreg = reg
        Gamma = block_diag(OAreg*OA_penalty, HCreg*HC_penalty)
    
        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(self.Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))
        P = block_diag(*(101*(P,)))
        
        # ND matrix
        NDmatr = np.matrix(ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T, Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), target)

        deltac = np.array(np.dot(RHS, LHS.I)).squeeze()
        # self.c = np.linalg.solve(LHS, RHS)
        self.c += deltac
        self.cOA = self.c[:nBinsOA]
        self.cHC = self.c[nBinsOA:]

        self.fluxPred = np.dot(self.ND, self.c)
        self.ratePred = np.dot(self.ND_rate, self.c)
        self.ratePred_statErr = np.sqrt(np.dot(self.ND_rate, np.power(self.c, 2)))

    def compressed_sensing(self, diag = [None], reg = 1.e-9, ND = None, target = None):
        if not ND:
            ND = self.ND
        if not target:
            target = self.target

        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        # A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        if any(diag):
            Gamma = np.diag(diag)
        else:
            A = np.eye(nBinsOA)
            Gamma = reg*A
        
        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))

        # ND matrix
        NDmatr = np.matrix(ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T, Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), target)

        self.c = np.array(np.dot(RHS, LHS.I)).squeeze()

        self.fluxPred = np.dot(self.ND, self.c)
        self.ratePred = np.dot(self.ND_rate, self.c)
        self.ratePred_statErr = np.sqrt(np.dot(self.ND_rate, np.power(self.c, 2)))

    def residual_norm(self, P = None, ND = None, target = None, **kwargs):
        if not np.any(P):
            P = self.P
        if not np.any(ND):
            ND = self.ND
        if not np.any(target):
            target = self.target

        return np.sqrt(np.sum(np.power(np.matmul(P, np.dot(ND, self.c) - target), 2)))
    
    def solution_norm(self, A = None, **kwargs):
        if not np.any(A):
            A = self.A
        return np.sqrt(np.sum(np.power(np.dot(A, self.c), 2)))

    def variance_norm(self):
        if self.ppfx_systs_loaded:
            return np.sqrt(np.sum(np.sum(np.power(np.matmul(self.P,
                                                            np.dot(self.ND_ppfx_univs,
                                                                   self.c) - self.target),
                                                  2))))

    def set_coeffs(self, other):
        self.coeffs = other.coeff

from utils import *
from fluxes import *
from oscProbs import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class flux_fitter:
    def __init__(self, beamMode, FDfromFlavor, FDtoFlavor, NDflavor, oscParam = None, ErebinF = 1, useHC = True):
        self.initNEbins = Ebins.size
        self.initNOAbins = OAbins.size

        self.Ebins = Ebins
        self.OAbins = OAbins

        if useHC:
            self.HCbins = currents
            self.maxHC = np.max(currents)
        else:
            self.HCbins = []
            self.maxHC = 293.

        self.maxOA = np.max(self.OAbins)

        self.rebinF = ErebinF
        self.Ebins = average(self.Ebins, self.rebinF)

        self.beamMode = beamMode
        self.FDfromFlavor = FDfromFlavor
        self.FDtoFlavor = FDtoFlavor
        self.NDflavor = NDflavor
        
        if not oscParam:
            self.Posc = oscProb(FDfromFlavor, FDtoFlavor).load(self.Ebins)
        else:
            self.Posc = oscParam.load(self.Ebins)
        
        self.load_nom()
            
        self.Ebounds = (0, np.max(self.Ebins))
        self.OutOfRegionFactors = (0, 0)

        self.ppfx_systs_loaded = False
        self.other_systs_loaded = False

    def rebin_nom(self, rebinF):
        self.ND_OA = average(self.ND_OA, rebinF)
        self.ND_HC = average(self.ND_HC, rebinF)
        self.ND_full = average(self.ND_full, rebinF)
        self.FD_unoscillated = average(self.FD_unoscillated, rebinF)
        self.FD_oscillated = average(self.FD_oscillated, rebinF)
        self.Posc = average(self.Posc, rebinF)
        self.target = average(self.target, rebinF)
        self.Ebins = average(self.Ebins, rebinF)
        
    def rebin_ppfx_systs(self, rebinF):
        self.FD_unosc_ppfx_univs = average(self.FD_unosc_ppfx_univs, rebinF)
        self.ND_OA_unosc_ppfx_univs = average(self.ND_OA_unosc_ppfx_univs, rebinF)
        self.ND_HC_unosc_ppfx_univs = average(self.ND_HC_unosc_ppfx_univs, rebinF)
        
    def rebin_other_systs(self, rebinF):
        pass
    
    def rebin(self, rebinF):
        self.rebinF = rebinF
        
        self.rebin_nom(rebinF)

        if self.ppfx_systs_loaded:
            self.rebin_ppfx_systs(rebinF)
        if self.other_systs_loaded:
            self.rebin_other_systs(rebinF)

    def load_FD_nom(self):
        self.FD_unoscillated = FD_nominal[self.beamMode][self.FDfromFlavor].load()
        self.FD_unoscillated = average(self.FD_unoscillated, self.rebinF)
        self.FD_oscillated = self.FD_unoscillated*self.Posc
        self.target = self.FD_oscillated

    def load_ND_OA_nom(self):
        self.ND_OA = ND_nominal[self.beamMode][self.NDflavor].load()[:,OAbins <= self.maxOA]
        self.ND_OA = average(self.ND_OA, self.rebinF)

        if "ND_HC" in dir(self):
            self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
            self.ND = self.ND_full
        
    def load_ND_HC_nom(self):
        if list(self.HCbins):
            self.ND_HC = np.array([ND_HC_shifts[self.beamMode][self.NDflavor][current].load()
                                   for current in self.HCbins]).T
        else:
            self.ND_HC = np.ndarray((self.initNEbins, 0))
        self.ND_HC = average(self.ND_HC, self.rebinF)
            
        if "ND_OA" in dir(self):
            self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
            self.ND = self.ND_full

    def load_nom(self):
        self.load_FD_nom()
        self.load_ND_OA_nom()
        self.load_ND_HC_nom()
        
    def load_FD_ppfx_systs(self, nUniv = 100):
        FD_CV = FD_ppfx_CV[self.beamMode][self.FDfromFlavor].load()
        FD = np.array([shift.load() for shift
                       in FD_ppfx_shifts[self.beamMode][self.FDfromFlavor][:nUniv]])
        
        FD /= FD_CV
        FD = average(FD, self.rebinF, axis = 1)
        FD *= self.FD_unoscillated

        self.FD_unosc_ppfx_univs = FD

    def load_ND_OA_ppfx_systs(self, nUniv = 100):
        ND_OA_CV = ND_ppfx_CV[self.beamMode][self.FDfromFlavor].load()[:, OAbins <= self.maxOA]
        ND_OA = np.array([shift.load() for shift
                          in ND_ppfx_shifts[self.beamMode][self.FDfromFlavor][:nUniv]])[:, :, OAbins <= self.maxOA]

        ND_OA /= ND_OA_CV
        ND_OA = average(ND_OA, self.rebinF, axis = 1)
        ND_OA *= self.ND_OA

        self.ND_OA_ppfx_univs = ND_OA
        if "ND_HC_ppfx_univs" in dir(self):
            self.ND_ppfx_univs = np.concatenate((self.ND_OA_ppfx_univs,
                                                 self.ND_HC_ppfx_univs),
                                                axis = 2)

    def load_ND_HC_ppfx_systs(self, nUniv = 100):
        ND_HC_CV = np.array([ND_HC_ppfx_CV[self.beamMode][self.FDfromFlavor][current].load()
                             for current in self.HCbins]).T
        if np.any(self.HCbins):
            ND_HC = np.array([np.array([ND_HC_ppfx_shifts[self.beamMode][self.FDfromFlavor][current][univ].load()
                                        for univ in range(nUniv)]).T
                              for current in self.HCbins]).T
        else:
            ND_HC = np.ndarray((nUniv, self.initNEbins, len(self.HCbins)))

        ND_HC /= ND_HC_CV
        ND_HC = average(ND_HC, self.rebinF, axis = 1)
        ND_HC *= self.ND_HC
        
        self.ND_HC_ppfx_univs = ND_HC
        if "ND_OA_ppfx_univs" in dir(self):
            self.ND_ppfx_univs = np.concatenate((self.ND_OA_ppfx_univs,
                                                 self.ND_HC_ppfx_univs),
                                                axis = 2)

    def load_ppfx_systs(self, nUniv = 100):
        self.load_FD_ppfx_systs(nUniv = nUniv)
        self.load_ND_OA_ppfx_systs(nUniv = nUniv)
        self.load_ND_HC_ppfx_systs(nUniv = nUniv)

        self.ppfx_systs_loaded = True

    def load_FD_other_systs(self):
        FD = {key: average(shift.load(), self.rebinF) for key, shift
              in FD_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.FD_unosc_other_univs = FD

    def load_ND_OA_other_systs(self):
        ND_OA = {key: average(shift.load(), self.rebinF) for key, shift
                 in ND_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.ND_OA_other_univs = ND_OA
        if "ND_HC_other_univs" in dir(self):
            self.ND_other_univs = np.concatenate((self.ND_OA_other_univs,
                                                  self.ND_HC_other_univs),
                                                 axis = 2)

    def load_ND_HC_other_systs(self):
        ND = {key: shift.load() for key, shift
              in ND_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.ND_HC_other_univs = ND_OA
        if "ND_OA_other_univs" in dir(self):
            self.ND_other_univs = np.concatenate((self.ND_OA_other_univs,
                                                  self.ND_HC_other_univs),
                                                 axis = 2)
        
    def load_other_systs(self):
        
        self.FD_unosc_other_univs = FD
        self.ND_other_univs = ND

        self.load_FD_other_systs()
        self.load_ND_OA_other_systs()
        self.load_ND_HC_other_systs()

        self.other_systs_loaded = True
        
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
                print "found peaks at ", self.Ebounds
    def set_OOR(self, weights):
        self.OutOfRegionFactors = weights

    def set_oscHypothesis(self, Posc):
        self.Posc = Posc
        self.FD_oscillated = self.FD_unoscillated*Posc
        self.target = self.FD_oscillated

    def set_maxOA(self, maxOA):
        self.maxOA = maxOA
        self.OAbins = OAbins[OAbins <= maxOA]
        self.load_ND_OA_nom()

        if self.ppfx_systs_loaded:
            self.load_ND_OA_ppfx_systs()

    def set_maxHC(self, maxHC):
        self.maxHC = maxHC
        # self.ND_HC = np.array([ND_HC_shifts[self.beamMode][self.NDflavor][current].load()
        #                        for current in currents])[self.HCbins <= maxHC]
        self.HCbins = self.HCbins[self.HCbins <= maxHC]

        self.load_ND_HC_nom()
        # self.ND_HC = average(self.ND_HC, rebinF)
        # self.ND_full = np.concatenate((self.ND_OA, self.ND_HC.T), axis = 1)

        if self.ppfx_systs_loaded:
            self.load_ND_HC_ppfx_systs()

    def use_currents(self, theseCurrents):
        self.maxHC = max(list(theseCurrents) + [293])
        if np.any(theseCurrents):
            self.ND_HC = np.array([8*ND_HC_shifts[self.beamMode][self.NDflavor][thisCurrent].load()
                                   for thisCurrent in theseCurrents]).T
            self.ND_HC = average(self.ND_HC, self.rebinF)
        else:
            self.ND_HC = np.ndarray((self.ND_HC.shape[0], 0))
        self.HCbins = np.array(theseCurrents)
        
        self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
        self.ND = self.ND_full

        if self.ppfx_systs_loaded:
            self.load_ND_HC_ppfx_systs()

    def calc_coeffs(self, OAreg, HCreg, ND = [None], target = [None], fluxTimesE = False):
        if not np.any(ND):
            ND = self.ND_full
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

    def calc_var_coeff_correction(self, reg, fluxTimesE = False):
        """
        Calculate a correction to the nominal coefficients to minimize the variance of the residual
        over different systematic universes
        """
        ND = np.concatenate((self.ND_full,) + tuple(ND for ND in self.ND_ppfx_univs))
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

    def residual_norm(self, P = None, ND = None, target = None, **kwargs):
        if not np.any(P):
            P = self.P
        if not np.any(ND):
            ND = self.ND_full
        if not np.any(target):
            target = self.target

        return np.sqrt(np.sum(np.power(np.matmul(P, np.dot(ND, self.c) - target), 2)))
    
    def solution_norm(self, A = None, **kwargs):
        if not np.any(A):
            A = self.A
        return np.sqrt(np.sum(np.power(np.dot(self.A, self.c), 2)))

    def variance_norm(self):
        if self.ppfx_systs_loaded:
            return np.sqrt(np.sum(np.sum(np.power(np.matmul(self.P,
                                                            np.dot(self.ND_ppfx_univs,
                                                                   self.c) - self.target),
                                                  2))))

    def set_coeffs(self, other):
        self.coeffs = other.coeff

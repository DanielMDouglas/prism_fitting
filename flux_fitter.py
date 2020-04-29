from utils import *
from fluxes import *
from oscProbs import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class flux_fitter:
    def __init__(self, beamMode, FDfromFlavor, FDtoFlavor, NDflavor, oscParam = None, rebin = 1, useHC = True):
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

        self.rebin = rebin
        if type(self.rebin) == np.ndarray:
            self.Ebins = average_by_bin_edge(self.Ebins, self.Ebins, self.rebin)
        elif type(self.rebin) == int:
            self.Ebins = average(self.Ebins, self.rebin)

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

    def rebin_nom(self):
        if type(self.rebin) == np.ndarray:
            self.ND_OA = average_by_bin_edge(self.ND_OA,
                                             self.Ebins,
                                             self.rebin)
            self.ND_HC = average_by_bin_edge(self.ND_HC,
                                             self.Ebins,
                                             self.rebin)
            self.ND_full = average_by_bin_edge(self.ND_full,
                                               self.Ebins,
                                               self.rebin)
            self.FD_unoscillated = average_by_bin_edge(self.FD_unoscillated,
                                                       self.Ebins,
                                                       self.rebin)
            self.FD_oscillated = average_by_bin_edge(self.FD_oscillated,
                                                     self.Ebins,
                                                     self.rebin)
            self.Posc = average_by_bin_edge(self.Posc,
                                            self.Ebins,
                                            self.rebin)
            self.target = average_by_bin_edge(self.target,
                                              self.Ebins,
                                              self.rebin)
            self.Ebins = average_by_bin_edge(self.Ebins,
                                             self.Ebins,
                                             self.rebin)
        elif type(self.rebin) == int:
            self.ND_OA = average(self.ND_OA,
                                 self.rebin)
            self.ND_HC = average(self.ND_HC,
                                 self.rebin)
            self.ND_full = average(self.ND_full,
                                   self.rebin)
            self.FD_unoscillated = average(self.FD_unoscillated,
                                           self.rebin)
            self.FD_oscillated = average(self.FD_oscillated,
                                         self.rebin)
            self.Posc = average(self.Posc,
                                self.rebin)
            self.target = average(self.target,
                                  self.rebin)
            self.Ebins = average(self.Ebins,
                                 self.rebin)
        
    def rebin_ppfx_systs(self):
        if type(self.rebin) == np.ndarray:
            self.FD_unosc_ppfx_univs = average_by_bin_edge(self.FD_unosc_ppfx_univs,
                                                           self.Ebins,
                                                           self.rebin)
            self.ND_OA_unosc_ppfx_univs = average_by_bin_edge(self.ND_OA_unosc_ppfx_univs,
                                                              self.Ebins,
                                                              self.rebin)
            self.ND_HC_unosc_ppfx_univs = average_by_bin_edge(self.ND_HC_unosc_ppfx_univs,
                                                              self.Ebins,
                                                              self.rebin)
        elif type(self.rebin) == int:
            self.FD_unosc_ppfx_univs = average(self.FD_unosc_ppfx_univs,
                                               self.rebin)
            self.ND_OA_unosc_ppfx_univs = average(self.ND_OA_unosc_ppfx_univs,
                                                  self.rebin)
            self.ND_HC_unosc_ppfx_univs = average(self.ND_HC_unosc_ppfx_univs,
                                                  self.rebin)
        
    def rebin_other_systs(self):
        pass
    
    def set_rebin(self, rebin):
        self.rebin = rebin

        if self.ppfx_systs_loaded:
            self.rebin_ppfx_systs()
        if self.other_systs_loaded:
            self.rebin_other_systs()
        
        self.rebin_nom()

    def load_FD_nom(self):
        self.FD_unoscillated = FD_nominal[self.beamMode][self.FDfromFlavor].load()
        if type(self.rebin) == np.ndarray:
            self.FD_unoscillated = average_by_bin_edge(self.FD_unoscillated,
                                                       Ebins,
                                                       self.rebin)
        elif type(self.rebin) == int:
            self.FD_unoscillated = average(self.FD_unoscillated,
                                           self.rebin)
        self.FD_oscillated = self.FD_unoscillated*self.Posc
        self.target = self.FD_oscillated

    def load_ND_OA_nom(self):
        self.ND_OA = ND_nominal[self.beamMode][self.NDflavor].load()[:,OAbins <= self.maxOA]
        if type(self.rebin) == np.ndarray:
            print self.ND_OA.shape
            self.ND_OA = average_by_bin_edge(self.ND_OA,
                                             Ebins,
                                             self.rebin)
            print self.ND_OA.shape
        elif type(self.rebin) == int:
            self.ND_OA = average(self.ND_OA, self.rebin)

        if "ND_HC" in dir(self):
            self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
            self.ND = self.ND_full
        
    def load_ND_HC_nom(self):
        if list(self.HCbins):
            self.ND_HC = np.array([ND_HC_shifts[self.beamMode][self.NDflavor][current].load()
                                   for current in self.HCbins]).T
        else:
            self.ND_HC = np.ndarray((self.initNEbins, 0))

        if type(self.rebin) == np.ndarray:
            self.ND_HC = 12*average_by_bin_edge(self.ND_HC,
                                                Ebins,
                                                self.rebin)
        elif type(self.rebin) == int:
            self.ND_HC = 12*average(self.ND_HC, self.rebin)
            
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
        if type(self.rebin) == np.ndarray:
            print FD.shape
            FD = average_by_bin_edge(FD, Ebins, self.rebin, axis = 1)
            print FD.shape
        elif type(self.rebin) == int:
            FD = average(FD, self.rebin, axis = 1)
            
        FD *= self.FD_unoscillated

        self.FD_unosc_ppfx_univs = FD

    def load_ND_OA_ppfx_systs(self, nUniv = 100):
        ND_OA_CV = ND_ppfx_CV[self.beamMode][self.FDfromFlavor].load()[:, OAbins <= self.maxOA]
        ND_OA = np.array([shift.load() for shift
                          in ND_ppfx_shifts[self.beamMode][self.NDflavor][:nUniv]])[:, :, OAbins <= self.maxOA]

        ND_OA /= ND_OA_CV
        if type(self.rebin) == np.ndarray:
            print ND_OA.shape
            ND_OA = average_by_bin_edge(ND_OA, Ebins, self.rebin, axis = 1)
            print ND_OA.shape
        elif type(self.rebin) == int:
            ND_OA = average(ND_OA, self.rebin, axis = 1)

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
            ND_HC = np.array([np.array([ND_HC_ppfx_shifts[self.beamMode][self.NDflavor][current][univ].load()
                                        for univ in range(nUniv)]).T
                              for current in self.HCbins]).T
        else:
            ND_HC = np.ndarray((nUniv, self.initNEbins, len(self.HCbins)))

        ND_HC /= ND_HC_CV

        if type(self.rebin) == np.ndarray:
            ND_HC = average_by_bin_edge(ND_HC, Ebins, self.rebin, axis = 1)
        elif type(self.rebin) == int:
            ND_HC = average(ND_HC, self.rebin, axis = 1)

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
        if type(self.rebin) == np.ndarray:
            FD = {key: average_by_bin_edge(shift.load(),
                                           Ebins,
                                           self.rebin) for key, shift
                  in FD_other_shifts[self.beamMode][self.FDfromFlavor].items()}
        elif type(self.rebin) == int:
            FD = {key: average(shift.load(), self.rebin) for key, shift
                  in FD_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.FD_unosc_other_univs = FD

    def load_ND_OA_other_systs(self):
        if type(self.rebin) == np.ndarray:
            ND_OA = {key: average_by_bin_edge(shift.load(),
                                              Ebins,
                                              self.rebin)[:, OAbins <= self.maxOA] for key, shift
                     in ND_other_shifts[self.beamMode][self.NDflavor].items()}
        elif type(self.rebin) == int:
            ND_OA = {key: average(shift.load(), self.rebin)[:, OAbins <= self.maxOA] for key, shift
                     in ND_other_shifts[self.beamMode][self.NDflavor].items()}
        
        self.ND_OA_other_univs = ND_OA
        if "ND_HC_other_univs" in dir(self):
            self.ND_other_univs = {key: np.concatenate((self.ND_OA_other_univs[key],
                                                        self.ND_HC_other_univs[key]),
                                                       axis = 1)
                                   for key in systKeys}

    def load_ND_HC_other_systs(self):
        # Don't have HC beam systs yet, so fake them by assuming
        # the same fractional error as seen in 293 kA
        if np.any(self.HCbins):
            ND_HC = {}
            for key in HCsystKeys:
                currentList = []
                for current in self.HCbins:
                    thisFlux = ND_HC_other_shifts[self.beamMode][self.NDflavor][current][key].load()
                    if type(self.rebin) == np.ndarray:
                        currentList.append(12*average_by_bin_edge(thisFlux, Ebins, self.rebin))
                    elif type(self.rebin) == int:
                        currentList.append(12*average(thisFlux, self.rebin))
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
            
    def calc_coeffs(self, OAreg, HCreg, ND = [None], target = [None], fluxTimesE = False, **kwargs):
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

    def calc_coeffs_DFT(self, OAreg, HCreg, filt, ND = [None], target = [None], fluxTimesE = False):
        if not np.any(ND):
            ND = self.ND_full
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
        # print np.diag(filt)*dft(nBinsOA)
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

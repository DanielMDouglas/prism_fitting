from utils import *
from fluxes import *
from oscProbs import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

class flux_fitter:
    def __init__(self, beamMode, FDfromFlavor, FDtoFlavor, NDflavor, oscParam = None, ErebinF = 1, useHC = True):
        self.beamMode = beamMode
        self.FDfromFlavor = FDfromFlavor
        self.FDtoFlavor = FDtoFlavor
        self.NDflavor = NDflavor

        self.FD_unoscillated = FD_nominal[beamMode][FDfromFlavor].load()
        self.ND_OA = ND_nominal[beamMode][NDflavor].load()
        if useHC:
            self.ND_HC = np.array([ND_HC_shifts[beamMode][NDflavor][current].load()
                                   for current in currents]).T
        else:
            print self.ND_OA.shape
            self.ND_HC = np.ndarray((self.ND_OA.shape[0], 0))
            print self.ND_HC.shape
        self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
        self.ND = self.ND_full

        self.Ebins = Ebins
        self.OAbins = OAbins
        if useHC:
            self.HCbins = currents
        else:
            self.HCbins = []
            
        self.Ebounds = (0, np.max(self.Ebins))
        self.OutOfRegionFactors = (0, 0)
        self.maxOA = np.max(self.OAbins)
        if useHC:
            self.maxHC = np.max(currents)
        else:
            self.maxHC = 293.
            
        if not oscParam:
            self.Posc = oscProb(FDfromFlavor, FDtoFlavor).load(self.Ebins)
        else:
            self.Posc = oscParam.load(self.Ebins)
            
        self.FD_oscillated = self.FD_unoscillated*self.Posc
        # add target shaping later
        # maybe don't need it, since target is weight-able?
        self.target = self.FD_oscillated

        self.rebin(ErebinF)

    def rebin(self, rebinF):
        self.rebinF = rebinF
        
        self.ND_OA = average(self.ND_OA, rebinF)
        self.ND_HC = average(self.ND_HC, rebinF)
        self.ND_full = average(self.ND_full, rebinF)
        self.FD_unoscillated = average(self.FD_unoscillated, rebinF)
        self.FD_oscillated = average(self.FD_oscillated, rebinF)

        self.target = average(self.target, rebinF)
        self.Ebins = average(self.Ebins, rebinF)
    def load_ppfx_systs(self, nUniv = 100):
        FD_CV = FD_ppfx_CV[self.beamMode][self.FDfromFlavor].load()
        FD = np.array([shift.load((self.Ebins,)) for shift
                       in FD_ppfx_shifts[self.beamMode][self.FDfromFlavor][:nUniv]])
        FD *= self.FD_unoscillated/FD_CV
        
        ND_OA_CV = ND_ppfx_CV[self.beamMode][self.FDfromFlavor].load()
        ND_OA = np.array([shift.load((self.Ebins, self.OAbins)) for shift
                          in ND_ppfx_shifts[self.beamMode][self.FDfromFlavor][:nUniv]])
        ND_OA *= self.ND_OA/ND_OA_CV

        ND_HC_CV = np.array([2*np.sum(ND_HC_ppfx_CV[self.beamMode][self.FDfromFlavor][current].load()[:,:4], axis = 1)
                             for current in currents]).T
        ND_HC = np.array([[2*np.sum(ND_HC_ppfx_shifts[self.beamMode][self.FDfromFlavor][current][univ].load((self.Ebins, self.OAbins))[:,:4], axis = 1)
                           for univ in range(nUniv)]
                          for current in currents]).reshape(nUniv, len(self.Ebins), len(self.HCbins))
        ND_HC = ND_HC*(self.ND_HC.T/ND_HC_CV)
        
        self.FD_unosc_ppfx_univs = FD
        self.ND_OA_ppfx_univs = ND_OA
        self.ND_HC_ppfx_univs = ND_HC
        self.ND_ppfx_univs = np.concatenate((ND_OA,
                                             ND_HC),
                                            axis = 2)

    def load_other_systs(self):
        FD = {key: shift.load((self.Ebins,)) for key, shift
              in FD_other_shifts[self.beamMode][self.FDfromFlavor].items()}
        ND = {key: shift.load((self.Ebins, self.OAbins)) for key, shift
              in ND_other_shifts[self.beamMode][self.FDfromFlavor].items()}

        self.FD_unosc_other_univs = FD
        self.ND_other_univs = ND
        
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
                smoothed = np.convolve(self.target,
                                       (1/float(window_size))*np.ones(window_size))[fringe_size:-fringe_size]
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
        self.ND_OA = ND_nominal[self.beamMode][self.NDflavor].load()[:,OAbins <= maxOA]

        self.ND_OA = average(self.ND_OA, self.rebinF)

        self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
        if 'ND_ppfx_univs' in dir(self):
            self.ND_ppfx_univs = self.ND_ppfx_univs[:, :, np.concatenate((OAbins <= self.maxOA, currents <= self.maxHC))]
            
    def set_maxHC(self, maxHC, rebinF = 1):
        self.maxHC = maxHC
        # self.ND_HC = np.array([2*np.sum(thisHC.load()[:,:4], axis = 1)
        #                        for thisHC in ND_HC_shifts[self.beamMode][self.NDflavor].values()])[currents <= maxHC]
        self.ND_HC = np.array([ND_HC_shifts[self.beamMode][self.NDflavor][current].load()
                               for current in currents])[self.HCbins <= maxHC]
        self.HCbins = self.HCbins[self.HCbins <= maxHC]
        
        self.ND_HC = np.sum(self.ND_HC[:, i::rebinF] for i in range(rebinF))/float(rebinF)
        self.ND_full = np.concatenate((self.ND_OA, self.ND_HC.T), axis = 1)

        if 'ND_ppfx_univs' in dir(self):
            self.ND_ppfx_univs = self.ND_ppfx_univs[:, :, np.concatenate((OAbins <= self.maxOA, currents <= self.maxHC))]

    def use_currents(self, theseCurrents):
        self.maxHC = max(list(theseCurrents) + [293])
        if np.any(theseCurrents):
            self.ND_HC = np.array([8*ND_HC_shifts[self.beamMode][self.NDflavor][thisCurrent].load()
                                   for thisCurrent in theseCurrents]).T
            self.ND_HC = average(self.ND_HC, self.rebinF)
        else:
            self.ND_HC = np.ndarray((self.ND_OA.shape[0], 0))
        self.HCbins = np.array(theseCurrents)
        
        self.ND_full = np.concatenate((self.ND_OA, self.ND_HC), axis = 1)
        self.ND = self.ND_full

        if 'ND_ppfx_univs' in dir(self):
            self.ND_ppfx_univs = self.ND_ppfx_univs[:, :, np.concatenate((OAbins <= self.maxOA, currents == theseCurrents))]

    def calc_coeffs(self, OAreg, HCreg, ND = [None], target = [None], fluxTimesE = False):
        if not np.any(ND):
            ND = self.ND_full
        if not np.any(target):
            target = self.target
        if fluxTimesE:
            ND = ND*self.Ebins.reshape(200, 1)
            target = target*self.Ebins
        
        # penalty matrix for coefficients
        nBinsOA = np.sum(self.OAbins <= self.maxOA)
        nBinsHC = self.ND_HC.shape[1]
        
        # # OA penalty term is a difference matrix
        # OA_penalty = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)
        # OA penalty term is the L1 norm
        OA_penalty = np.eye(nBinsOA)
        HC_penalty = np.eye(nBinsHC)
        Gamma = block_diag(OAreg*OA_penalty, HCreg*HC_penalty)

        # weighting matrix for target flux
        P = np.diag(np.where(self.Ebins > self.Ebounds[1],
                             self.OutOfRegionFactors[1],
                             np.where(self.Ebins < self.Ebounds[0],
                                      self.OutOfRegionFactors[0],
                                      1)))

        # ND matrix
        NDmatr = np.matrix(ND)
        LHS = np.matmul(NDmatr.T, np.matmul(P, NDmatr)) + np.matmul(Gamma.T, Gamma)
        RHS = np.dot(np.matmul(NDmatr.T, P), target)

        self.c = np.array(np.dot(RHS, LHS.I)).squeeze()
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
        
    def set_coeffs(self, other):
        self.coeffs = other.coeffs

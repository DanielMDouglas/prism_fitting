from utils import *

import os

class flux:
    def __init__(self, infileName, branchName):
        self.infileName = infileName
        self.branchName = branchName
    def load(self, bins = []):
        content, stat = root_to_array(self.infileName, self.branchName, bins = bins)
        return content

# Everything from here on depends on the tree structure
# of the flux file you're using.
# Please be aware of this and update this accordingly
# if you switch to using different fluxes

# Currently, this program expects a single large file
# which contains the nominal ND and FD fluxes
# as well as a large number of systematic variations

# Right now, only the nominal fluxes are loaded,
# but please add the systematic shift definitions here
# as they are needed

if 'DP_FLUX_FILE' in os.environ:
    nomFileName = os.environ['DP_FLUX_FILE']
else:
    # nomFileName = "../flux/syst/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_fine.root"
    # nomFileName = "../flux/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_finebin.root"
    # nomFileName = "../flux/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_finebin_HHCOnAxis.root"
    nomFileName = "../flux/All_HC.root"
    print "[WARNING] Environment variable DP_FLUX_FILE is unset! Using default location:"
    print nomFileName

if 'DP_SYST_FILE' in os.environ:
    systFileName = os.environ['DP_SYST_FILE']
else:
    # systFileName = "../flux/syst/FluxErrors_40mOffAxis_Total_BothBeamModes_AllSpecies.root"
    systFileName = nomFileName
    print "[WARNING] Environment variable DP_SYST_FILE is unset! Using default location:"
    print systFileName

# HCFileName = "../flux/syst/HigherHCFluxes.Fine.root"
HCFileName = nomFileName

flavors = ["nue", "numu", "nuebar", "numubar"]
# modes = ["nu", "nubar"]
modes = ["nu"]
currentNames = ["200",
                "250",
                "280",
                "295.5",
                "298",
                "300.5",
                "303",
                "305.5",
                "308",
                "310.5",
                "313",
                "323",
                "333",
                "343"]
currents = np.array([float(i) for i in currentNames])
systKeys = ["HC_p1",
            "HC_m1",
            "Horn1_XShift",
            "Horn1_YShift",
            "Horn1_XNegShift",
            "Horn1_X3mmShift",
            "Horn1_XNeg3mmShift",
            "Horn2_XShift",
            "Horn2_YShift",
            "Horn2_XNegShift",
            "WL_p1",
            "DPR_p1",
            "TargetDensity_p1",
            "TargetDensity_m1",
            "BeamOffsetX_p1",
            "BeamOffsetX_m1",
            "BeamSigma_p1",
            "BeamSigma_m1",
            "BeamTheta_p1",
            "BeamThetaPhi_p1"]

# Using central value (CV) as nominal for now...
# ND_nominal = {mode: {flavor: flux(nomFileName, "ND_"+mode+"_ppfx/LBNF_"+flavor+"_flux_Nom")
#                      for flavor in ["nue", "numu", "nuebar", "numubar"]}
#               for mode in ["nu", "nubar"]}
ND_nominal = {mode: {flavor: flux(nomFileName, "ND_"+mode+"_ppfx/LBNF_"+flavor+"_flux_Nom")
                     for flavor in flavors}
              for mode in modes}

ND_HC_shifts = {mode: {flavor: {current: flux(HCFileName, "ND_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_Nom")
                                for current, currentName in zip(currents, currentNames)}
                       for flavor in flavors}
                for mode in modes}
ND_HC_ppfx_shifts = {mode: {flavor: {current: [flux(HCFileName, "ND_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_univ_"+str(i))
                                               for i in range(100)]
                                     for current, currentName in zip(currents, currentNames)}
                            for flavor in flavors}
                     for mode in modes}
ND_HC_ppfx_CV = {mode: {flavor: {current: flux(HCFileName, "ND_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_CV")
                                 for current, currentName in zip(currents, currentNames)}
                        for flavor in flavors}
                 for mode in modes}

# ND_shifts = {mode: {flavor: [flux(systFileName,
#                                   "EffectiveFluxParameters/param_"+str(i)+"/ND_"+mode+"_"+flavor)
#                              for i in range(50)]
#                     for flavor in ["nue", "numu", "nuebar", "numubar"]}
#              for mode in ["nu", "nubar"]}
ND_ppfx_shifts = {mode: {flavor: [flux(systFileName,
                                       "ND_"+mode+"_ppfx/LBNF_"+flavor+"_flux_univ_"+str(i))
                                  for i in range(100)]
                         for flavor in flavors}
                  for mode in modes}
ND_ppfx_CV = {mode: {flavor: flux(systFileName,
                                  "ND_"+mode+"_ppfx/LBNF_"+flavor+"_flux_CV")
                     for flavor in flavors}
              for mode in modes}
                 
ND_other_shifts = {mode: {flavor: {key: flux(systFileName,
                                             "ND_"+mode+"_"+key+"/LBNF_"+flavor+"_flux")
                                   for key in systKeys}
                          for flavor in flavors}
                   for mode in modes}
                                             
FD_nominal = {mode: {flavor: flux(nomFileName, "FD_"+mode+"_ppfx/LBNF_"+flavor+"_flux_Nom")
                     for flavor in flavors}
              for mode in modes}

FD_HC_shifts = {mode: {flavor: {current: flux(HCFileName, "FD_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_Nom")
                                for current, currentName in zip(currents, currentNames)}
                       for flavor in flavors}
                for mode in modes}
FD_HC_ppfx_shifts = {mode: {flavor: {current: [flux(HCFileName, "FD_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_univ_"+str(i))
                                               for i in range(100)]
                                     for current, currentName in zip(currents, currentNames)}
                            for flavor in flavors}
                     for mode in modes}
FD_HC_ppfx_CV = {mode: {flavor: {current: flux(HCFileName, "FD_"+mode+"_HC_"+str(currentName)+"/LBNF_"+flavor+"_flux_CV")
                                 for current, currentName in zip(currents, currentNames)}
                        for flavor in flavors}
                 for mode in modes}

FD_ppfx_shifts = {mode: {flavor: [flux(systFileName,
                                       "FD_"+mode+"_ppfx/LBNF_"+flavor+"_flux_univ_"+str(i))
                                  for i in range(100)]
                         for flavor in flavors}
                  for mode in modes}
FD_ppfx_CV = {mode: {flavor: flux(systFileName,
                                  "FD_"+mode+"_ppfx/LBNF_"+flavor+"_flux_CV")
                     for flavor in flavors}
              for mode in modes}
FD_other_shifts = {mode: {flavor: {key: flux(systFileName,
                                             "FD_"+mode+"_"+key+"/LBNF_"+flavor+"_flux")
                                   for key in systKeys}
                          for flavor in flavors}
                   for mode in modes}

# The binning within this file is consistent,
# so get them from any old histogram
# but this is not always the case!

Ebins, OAbins, zbins = root_to_axes(nomFileName, "ND_nu_ppfx/LBNF_numu_flux_Nom")

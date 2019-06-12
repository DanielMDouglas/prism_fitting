from utils import *

import os

class flux:
    def __init__(self, infileName, branchName):
        self.infileName = infileName
        self.branchName = branchName
    def load(self):
        content, stat = root_to_array(self.infileName, self.branchName)
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
    systFileName = os.environ['DP_FLUX_FILE']
else:
    systFileName = "../flux/syst/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_fine.root"
    print "[WARNING] Environment variable DP_FLUX_FILE is unset! Using default location:"
    print systFileName

ND_nominal = {mode: {flavor: flux(systFileName, "ND_"+mode+"_ppfx/LBNF_"+flavor+"_flux_Nom")
                     for flavor in ["nue", "numu", "nuebar", "numubar"]}
              for mode in ["nu", "nubar"]}

FD_nominal = {mode: {flavor: flux(systFileName, "FD_"+mode+"_ppfx/LBNF_"+flavor+"_flux_Nom")
                     for flavor in ["nue", "numu", "nuebar", "numubar"]}
              for mode in ["nu", "nubar"]}

# The binning within this file is consistent,
# so get them from any old histogram
# but this is not always the case!

Ebins, OAbins, zbins = root_to_axes(systFileName, "ND_nu_ppfx/LBNF_numu_flux_Nom")

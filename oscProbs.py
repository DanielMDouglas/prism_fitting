from utils import *
import os

# previously computed probabilities are stored here
# to avoid recalculating them every time
if 'DP_OSC_PROB_DIR' in os.environ:
    oscProbBaseDir = os.environ['DP_OSC_PROB_DIR']
else:
    oscProbBaseDir = "../oscProb/"
    print "[WARNING] Environment variable DP_OSC_PROB_DIR is unset! Using default location:"
    print oscProbBaseDir

# Make sure that DUNEPrismTools is installed
# And its setup.sh has been sourced
if 'DUNEPRISMTOOLSROOT' in os.environ:
    dunePrismToolsRoot = os.environ['DUNEPRISMTOOLSROOT']
else:
    raise EnvironmentError("Environment variable DUNEPRISMTOOLSROOT is unset!\nPlease make sure that DUNEPrismTools is installed and setup.")

# Oscillation helper needs a flux file
# get this in the same way as in fluxes.py
if 'DP_FLUX_FILE' in os.environ:
    systFileName = os.environ['DP_FLUX_FILE']
else:
    systFileName = "../flux/syst/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_fine.root"
    print "[WARNING] Environment variable DP_FLUX_FILE is unset! Using default location:"
    print systFileName

pdgCodes = {"nue": "12",
            "numu": "14",
            "nutau": "16",
            "nuebar": "-12",
            "numubar": "-14",
            "nutaubar": "-16"}

class oscProb:
    def __init__(self, fromFlavor, toFlavor, s23 = 0.5, dm32 = 2.6e-3, dcp = np.pi):
        
        s23 = str(s23)
        dm32 = str(dm32)
        dcp = str(dcp)

        self.branchName = "POsc"

        subdir = s23+"_"+dm32+"_"+dcp
        filename = fromFlavor+"_to_"+toFlavor+".root"
        self.infileName = oscProbBaseDir + subdir+"/"+filename

        if not subdir in os.listdir(oscProbBaseDir):
            os.mkdir(oscProbBaseDir+subdir)
        if not filename in os.listdir(oscProbBaseDir+subdir):
            exe = dunePrismToolsRoot+"/bin/dp_OscillateFlux"
            # arguments for dp_OscillatedFlux:
            #    -d dip angle (degrees)
            #    -p oscillation parameters (comma separated):
            #        sin^2(theta_{12})
            #        sin^2(theta_{13})
            #        sin^2(theta_{23})
            #        delta m^2_{12}
            #        delta m^2_{32}
            #        delta cp
            #    -i input flux (not used, but required by dp_OscillateFlux)
            #    -o output root file
            #    -n PDG codes for oscillation mode (e.g., numu -> nue)
            args = ["-d 5.8",
                    "-p 0.297,0.0215,"+s23+",7.37E-5,"+dm32+","+dcp,
                    "-i "+systFileName+",FD_nu_ppfx/LBNF_numu_flux_Nom",
                    "-o "+oscProbBaseDir+subdir+"/"+fromFlavor+"_to_"+toFlavor+".root,tmp",
                    "-n "+pdgCodes[fromFlavor]+","+pdgCodes[toFlavor]
            ]
            cmdString = " ".join([exe] + args)
            
            os.system(cmdString)
        
    def load(self, Ebins):
        return tgraph_to_array(self.infileName, self.branchName, Ebins)

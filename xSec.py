from __future__ import print_function

from utils import *
from fluxes import *

class crossSection:
    def __init__(self, infileName, branchName):
        self.infileName = infileName
        self.branchName = branchName
    def load(self, binEdges = [EbinEdges]):
        binCenters = 0.5*(binEdges[0][1:] + binEdges[0][:-1])
        content, stat = root_to_array(self.infileName,
                                      self.branchName,
                                      method = "interpolate",
                                      binEdges = binEdges)
        return binCenters*content

# Everything from here on depends on the tree structure
# of the cross section file you're using.
# Please be aware of this and update this accordingly
# if you switch to using different cross sections

# Currently, this program expects a single large file
# which contains the cross sections for nue, numu, nuebar, and numubar

if 'DP_XSEC_FILE' in os.environ:
    xSecFileName = os.environ['DP_XSEC_FILE']
else:
    xSecFileName = "../Ar40_xSec.root"
    print ("[WARNING] Environment variable DP_FLUX_FILE is unset! Using default location:")
    print (nomFileName)

xSec = {"nue": {"CCInc": crossSection(xSecFileName, "Ar40.nu12.E_0_10GeV.sigmaenu.root/CCInc"),
                "NCInc": crossSection(xSecFileName, "Ar40.nu12.E_0_10GeV.sigmaenu.root/NCInc"),
                "total": crossSection(xSecFileName, "Ar40.nu12.E_0_10GeV.sigmaenu.root/TotalXSec")},
        "numu": {"CCInc": crossSection(xSecFileName, "Ar40.nu14.E_0_10GeV.sigmaenu.root/CCInc"),
                 "NCInc": crossSection(xSecFileName, "Ar40.nu14.E_0_10GeV.sigmaenu.root/NCInc"),
                 "total": crossSection(xSecFileName, "Ar40.nu14.E_0_10GeV.sigmaenu.root/TotalXSec")},
        "nuebar": {"CCInc": crossSection(xSecFileName, "Ar40.nu-12.E_0_10GeV.sigmaenu.root/CCInc"),
                   "NCInc": crossSection(xSecFileName, "Ar40.nu-12.E_0_10GeV.sigmaenu.root/NCInc"),
                   "total": crossSection(xSecFileName, "Ar40.nu-12.E_0_10GeV.sigmaenu.root/TotalXSec")},
        "numubar": {"CCInc": crossSection(xSecFileName, "Ar40.nu-14.E_0_10GeV.sigmaenu.root/CCInc"),
                    "NCInc": crossSection(xSecFileName, "Ar40.nu-14.E_0_10GeV.sigmaenu.root/NCInc"),
                    "total": crossSection(xSecFileName, "Ar40.nu-14.E_0_10GeV.sigmaenu.root/TotalXSec")}}

NDvol = 4*2*3*(100**3) # cm^3
# FDvol = 6.2 * 11 * 11.96 * (100**3) # cm^3
FDvol = 12.0 * 14.0 * 58.2 * (100**3) # cm^3

LArDens = 1.3954 # g cm^-3
NDmass = NDvol*LArDens # g
FDmass = FDvol*LArDens # g

ArMass = 39.948*1.66054e-24 # g atom^-1
ArNnuc = 40. # nucleon atom^-1
MassPerNuc = ArMass/ArNnuc # g nucleon^-1
NDnNuc = NDmass/MassPerNuc # nucleon
FDnNuc = FDmass/MassPerNuc # nucleon

TarRate = 1.1e21 # POT/year

NDOAruntimeSplit = np.array(4*[2*14*7/365.25] + 77*[12/365.25])
NDHCruntimeSplit = np.array([7/365.25])

T = 3.5 # year
# T = 10
POT = TarRate*T # POT
    
NDscale = NDnNuc*POT/4 # nucleon POT
FDscale = FDnNuc*POT/2 # nucleon POT
NDtoFD = FDscale/NDscale

# # Vol ND
# LAr_dens = 1.3954 * (1E-3 / 1E-6) # g/cm3 -> kg/m^3
# vol_ND = 0.5 * 2 * 3              # x * y * z FV m^3
# # mass_ND = vol_ND * LAr_dens
# mass_ND = 67.e3
# massNucKg = 39.9623831238 * (1.998467052E-26 / 12.0)
# # Vol FD
# vol_FD = 6.2 * 11 * 11.96 # x * y * z FV m^3
# # mass_FD = 4* vol_FD * LAr_dens
# mass_FD = 40.e6
# # bw = 0.05
# # bw = FD_flux->GetXaxis()->GetBinWidth(i + 1)
# # E = FD_flux->GetXaxis()->GetBinCenter(i + 1)
# # Flux = FD_flux->GetBinContent(i + 1) * bw
# # FDEvRate_POT(i) = Flux * xsec_ccinc->Interpolate(E) * calc->P(14, 14, E) * mass_FD * 40.0 / MassNucKg

# T = 3.5 # years
# POT = 1.1e21 # POT/year

# NDscale = mass_ND*POT*T/massNucKg
# FDscale = mass_FD*POT*T/massNucKg
# NDtoFD = mass_FD/mass_ND



from utils import *
import sys
import os

# Make sure that Prob3 is installed
if 'PROB3ROOT' in os.environ:
    Prob3Root = os.environ['PROB3ROOT']
else:
    Prob3Root = os.path.join(os.getcwd(), "Prob3")
    os.environ['PROB3ROOT'] = Prob3Root
    print "[WARNING] Environment variable PROB3ROOT is unset! Using default location:"
    print Prob3Root

sys.path.append(os.path.dirname(Prob3Root))    
from Prob3 import BargerPropagator

prob3Codes = {"nue": 1,
              "numu": 2,
              "nutau": 3,
              "nuebar": -1,
              "numubar": -2,
              "nutaubar": -3}

class oscProb:
    def __init__(self, fromFlavor, toFlavor,
                 s12 = 0.297, s13 = 0.0215, s23 = 0.53,
                 dm21 = 7.37e-5, dm32 = 2.46e-3, dcp = -np.pi/2):

        DipAngle_degrees = 5.8
        LengthParam = np.cos(np.radians(90.0 + DipAngle_degrees))

        self.s12 = s12
        self.s13 = s13
        self.s23 = s23
        self.dm21 = dm21
        self.dm32 = dm32
        self.dcp = dcp
        
        self.bp = BargerPropagator()
        self.bp.DefinePath(LengthParam, 0)

        self.fromFlavor = prob3Codes[fromFlavor]
        self.toFlavor = prob3Codes[toFlavor]
        
        
    def load(self, Ebins):
        P = []
        for E in Ebins:
            self.bp.SetMNS(self.s12, self.s13, self.s23,
                           self.dm21, self.dm32, self.dcp,
                           E, True, self.fromFlavor)
            self.bp.propagate(self.toFlavor)
            prob = self.bp.GetProb(self.fromFlavor, self.toFlavor)
            P.append(prob)
        return np.array(P)

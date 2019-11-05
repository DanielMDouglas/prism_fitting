from flux_fitter import *
from plots import *

if __name__ == '__main__':

    # make a new fitter object
    # use FHC
    # use numu -> numu at FD
    # use numu at ND
    fitter = flux_fitter("nu", "numu", "numu", "numu", useHC = False)

    # use the following oscillation parameters
    dcp = -np.pi/2
    s23 = 0.5
    dm32 = 3.0e-3
    # a nicely formatted string for plotting
    dm32Pretty = r'$2.4 \times 10^{-3}$'

    # get the numu -> numu probability for these parameters
    osc_hyp = oscProb("numu", "numu", s23 = s23, dm32 = dm32, dcp = dcp)
    # load them into the fitter
    fitter.set_oscHypothesis(osc_hyp.load(Ebins))

    # we want to fit between the first and fourth highest energy peaks
    fitter.set_fit_region(peaks = [1, 4])
    # with a 5% effect on the low-E side and no weight on the high-E side
    fitter.set_OOR([0.05, 0])
    # find the coefficients with a regularization factor of 8.e-9
    fitter.calc_coeffs(8.e-9, 8.e-9)

    # a nicely formatted title 
    title = r'$\sin^2 \theta_{23} = $'+str(s23)+r'$, \Delta m^2_{32} = $'+dm32Pretty
    # plot the fit and the (fit - target)/target_unosc
    FD_flux_plot(fitter)
    ND_flux_slice_plot(fitter)
    fit_and_ratio_plot(fitter, title = title)
    coeff_plot(fitter, HC = False)

    plt.show()    
    

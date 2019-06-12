from flux_fitter import *

if __name__ == '__main__':

    fitter = flux_fitter("nu", "numu", "numu", "numu")
    fitter.set_OOR([0.05, 0])
    
    dcp = -np.pi/2
    s23 = 0.5
    dm32 = 2.5e-3
    dm32Pretty = r'$2.4 \times 10^{-3}$'

    osc_hyp = oscProb("numu", "numu", s23 = s23, dm32 = dm32, dcp = dcp)
    fitter.set_oscHypothesis(osc_hyp.load(Ebins))
    fitter.set_fit_region(peaks = [1, 4])
    fitter.calc_coeffs(8.e-9)

    title = r'$\sin^2 \theta_{23} = $'+str(s23)+r'$, \Delta m^2_{32} = $'+dm32Pretty

    fitter.plot_target_fit_and_ratio(title = title)
    fitter.plot_coeffs(title = title)

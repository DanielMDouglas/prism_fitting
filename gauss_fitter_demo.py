import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ROOT import TH1D, TH2D, TFile

def root_to_array(infileName, branchName):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    arr = np.ndarray(shape)
    for i in xrange(shape[0]):
        for j in xrange(shape[1]):
            for k in xrange(shape[2]):
                arr[i,j,k] = TH.GetBinContent(i+1, j+1, k+1)
    infile.Close()
    return arr.squeeze()

def root_to_axes(infileName, branchName):
    infile = TFile(infileName)
    TH = infile.Get(branchName)
    shape = [TH.GetNbinsX(),
             TH.GetNbinsY(),
             TH.GetNbinsZ()]
    mins = [TH.GetXaxis().GetXmin(),
            TH.GetYaxis().GetXmin(),
            TH.GetZaxis().GetXmin()]
    maxs = [TH.GetXaxis().GetXmax(),
            TH.GetYaxis().GetXmax(),
            TH.GetZaxis().GetXmax()]
    axes = [np.linspace(min_i, max_i, nBins_i)
            for min_i, max_i, nBins_i
            in zip(mins, maxs, shape)]
    return tuple(axes)

Ebins, OAbins, zbins = root_to_axes("../flux/nominal_wider.50cm.FHC.uniform.root", "LBNF_numu_flux")
nBinsE = len(Ebins)
nBinsOA = len(OAbins)

gauss_loc = 1
target_vect = st.norm.pdf(Ebins, loc = gauss_loc, scale = 0.1*gauss_loc)

plt.clf()
plt.plot(Ebins, target_vect)
plt.xlim(Ebins[0], Ebins[-1])
plt.xlabel(r'$E$ (GeV)')
plt.ylabel(r'$\Phi$ (arb.)$')
plt.grid()
plt.title("Gaussian Target")
plt.show()

ND_matrix = np.matrix(root_to_array("../flux/nominal_wider.50cm.FHC.uniform.root", "LBNF_numu_flux"))

plt.clf()
plt.imshow(ND_matrix.T, origin='lower')
ax = plt.gca()
ax.set_xticks([0, len(Ebins)/2, len(Ebins)])
ax.set_xticklabels([0, 5, 10])
ax.set_yticks([0, len(OAbins)])
ax.set_yticklabels([0, 33])
plt.colorbar()
plt.xlabel(r'$E$ (GeV)')
plt.ylabel(r'$D_{OA}$ (m)')
plt.title("ND off-axis spectrum")
plt.show()

def plot_fit(coeffs):
    # fig = plt.figure()
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    
    # ax1 = fig.add_subplot(221)
    ax1.plot(Ebins, target_vect, label = 'Target')
    ax1.plot(Ebins, np.dot(ND_matrix, c), label = 'Fit')
    ax1.set_xlim(Ebins[0], Ebins[-1])
    # ax1.set_xlabel(r'$E$ (GeV)')
    ax1.set_ylabel(r'$\Phi$ (arb.)')
    ax1.grid()
    ax1.legend()
    ax1.set_title("Flux fits")

    ax5 = fig.add_subplot(122)
    ax5.plot(OAbins, coeffs)
    ax5.set_xlabel(r'$D_{OA}$ (m)')
    ax5.set_xlim(OAbins[0], OAbins[-1])
    ax5.grid()
    ax5.set_title('Coefficients')

    # ax3 = fig.add_subplot(223)
    ax3.plot(Ebins, np.zeros_like(Ebins), linestyle='--')
    ax3.plot(Ebins, np.array(np.dot(ND_matrix, c)).flatten() - target_vect)
    ax3.set_xlim(Ebins[0], Ebins[-1])
    ax3.set_xlabel(r'$E$ (GeV)')
    ax3.set_ylabel(r'Residual (arb.)')
    ax3.grid()

    ax2.axis('off')
    ax4.axis('off')

    fig.subplots_adjust(hspace=0)
    
    plt.show()
    

c = np.dot(ND_matrix.I, target_vect.T).T

plot_fit(c)

D = np.diag(nBinsOA*[2]) - np.diag((nBinsOA - 1)*[1], k = 1) - np.diag((nBinsOA - 1)*[1], k = -1)
A = np.diag(nBinsOA*[1]) - np.diag((nBinsOA - 1)*[1], k = 1)

l = 1.e-12

def coefficients(flux_matrix, target, l):
    Gamma = l*A
    return np.dot(np.dot((np.dot(flux_matrix.T,flux_matrix) + np.dot(Gamma.T, Gamma)).I, flux_matrix.T), target).T


c = coefficients(ND_matrix, target_vect, l)

plot_fit(c)

def solution_norm(lamb):
    c = coefficients(ND_matrix, target_vect, lamb)
    return 0.5*np.log(np.sum(np.power(np.dot(A, c),2)))

def residual_norm(lamb):
    c = coefficients(ND_matrix, target_vect, lamb)
    return 0.5*np.log(np.sum(np.power(np.dot(ND_matrix, c).T - target_vect, 2)))

l_space = np.logspace(-15, -5, 1000)
sol_norm = np.array([solution_norm(l) for l in l_space])
res_norm = np.array([residual_norm(l) for l in l_space])
    
plt.clf()
plt.plot(res_norm, sol_norm)
plt.xlabel(r'Residual norm $|\Phi^{ND} \vec{c} - \Phi^{FD}|$')
plt.ylabel(r'Solution norm $|\Gamma \vec{c}|$')
plt.grid()
plt.title("L-curve")
plt.show()

dl = np.diff(l_space)
xi = sol_norm
rho = res_norm

xi_prime = np.diff(xi)/dl
rho_prime = np.diff(rho)/dl

xi_prime_prime = np.diff(xi_prime)/dl[:-1]
rho_prime_prime = np.diff(rho_prime)/dl[:-1]

curv = 2*(rho_prime[:-1]*xi_prime_prime - rho_prime_prime*xi_prime[:-1])/np.power(np.power(rho_prime[:-1], 2) + np.power(xi_prime[:-1], 2), 3./2)

opt_lambda = l_space[1:-1][curv==np.max(curv)][0]

print r'optimum l =', opt_lambda

plt.clf()
plt.plot(l_space[1:-1], curv)
plt.axvline(x=opt_lambda, ls = '--', c = 'r')
plt.xlabel(r'$\lambda$')
plt.grid()
plt.title(r'Curvature')
plt.semilogx()
plt.show()

Gamma = opt_lambda*A
c = np.dot(np.dot((np.dot(ND_matrix.T,ND_matrix) + np.dot(Gamma.T, Gamma)).I, ND_matrix.T), target_vect).T

opt_sol_norm = 0.5*np.log(np.sum(np.power(np.dot(Gamma, c),2)))
opt_res_norm = 0.5*np.log(np.sum(np.power(np.dot(ND_matrix, c) - target_vect, 2)))

plt.clf()
plt.plot(res_norm, sol_norm)
# plt.scatter(opt_res_norm, opt_sol_norm, label=str(opt_lambda))
plt.scatter(residual_norm(opt_lambda), solution_norm(opt_lambda), label = str(opt_lambda))
plt.scatter(residual_norm(1.e-12), solution_norm(1.e-12), label = '1.e-12')
plt.xlabel(r'Residual norm $|\Phi^{ND} \vec{c} - \Phi^{FD}|$')
plt.ylabel(r'Solution norm $|\Gamma \vec{c}|$')
plt.grid()
plt.title("L-curve")
plt.legend()
plt.show()
    
plot_fit(c)

__author__ = 'kwebb'

import numpy as np
import os
from scipy.optimize import curve_fit
import glob
from astropy.io import fits
import ppxf_util as util
from matplotlib import pyplot as plt


def remove_absorp_lines(bin_sci, ppxf_bestfit, em_bin_sci, plot=False):
    """
    """

    for j in range(len(glob.glob(bin_sci.format('*')))):
        if not os.path.exists(em_bin_sci.format(j)):
            abs_fit, ems_fit = fit_spectra(bin_sci.format(j), ppxf_bestfit.format(j), plot)

            ems_hdu = fits.PrimaryHDU()
            ems_hdu.data = ems_fit
            ems_hdu.header = fits.getheader(bin_sci.format(j), 0)
            ems_hdu.writeto(em_bin_sci.format(j))


def remove_emission_lines(bin_sci, abs_bin_sci, ppxf_bestfit, plot=False):
    """
    """

    for j in range(len(glob.glob(bin_sci.format('*')))):
        if not os.path.exists(abs_bin_sci.format(j)):

            if plot:
                print('>>>>> Removing emission lines from spectra {}'.format(j))

            abs_fit, ems_fit = fit_spectra(bin_sci.format(j), ppxf_bestfit.format(j), plot)

            abs_hdu = fits.PrimaryHDU()
            abs_hdu.data = abs_fit
            abs_hdu.header = fits.getheader(bin_sci.format(j), 0)
            abs_hdu.writeto(abs_bin_sci.format(j))


def fit_spectra(bin_sci_j, ppxf_bestfit_j, plot=False):
    """
    """

    with fits.open(bin_sci_j) as hdu:
        odata = hdu[0].data
        ohdr = hdu[0].header

    bestfit = fits.getdata(ppxf_bestfit_j)

    lamRange = ohdr['CRVAL1'] + np.array([1. - ohdr['CRPIX1'], ohdr['NAXIS1'] - ohdr['CRPIX1']]) * ohdr['CD1_1']

    # Log bin the spectra to match the best fit absorption template
    galaxy, logLam1, velscale = util.log_rebin(lamRange, odata)
    log_bins = np.exp(logLam1)
    emlns, lnames, lwave = util.emission_lines(logLam1, lamRange, 2.3)

    # Hgamma 4340.47, Hbeta 4861.33, OIII [4958.92, 5006.84]
    # lwave = [4340.47, 4861.33, 4958.92, 5006.84]

    # find the index of the emission lines
    iHg = (np.abs(log_bins - lwave[0])).argmin()
    iHb = (np.abs(log_bins - lwave[1])).argmin()
    iOIII = (np.abs(log_bins - lwave[2])).argmin()
    iOIIIb = (np.abs(log_bins - (lwave[2] - 47.92))).argmin()

    # There are BOTH absorption and emission features about the wavelength of Hgamma and Hbeta, so we need
    # to use a specialized fitting function (convolved Guassian and Lorentzian -> pVoight) to remove the
    # emission lines
    popt_Hg, pcov_Hg = fit_gaussian_pvoightcont(galaxy, iHg, bestfit, log_bins)
    popt_Hb, pcov_Hb = fit_gaussian_pvoightcont(galaxy, iHb, bestfit, log_bins)

    # There are only emission features about the OIII doublet so we only fit the emission line with a Gaussian
    popt_OIII, pcov_OIII = fit_gaussian_lorentz_lincont(galaxy, iOIII, bestfit, log_bins, 5042, [5039, 5045])
    popt_OIIIb, pcov_OIIIb = fit_gaussian_lorentz_lincont(galaxy, iOIIIb, bestfit, log_bins, 5042, [5039, 5045])

    x = np.linspace(lamRange[0], lamRange[1], ohdr['NAXIS1'])

    em_fit = gaussian(x, popt_Hg[0], popt_Hg[1], popt_Hg[2], popt_Hg[3]) + \
             gaussian(x, popt_Hb[0], popt_Hb[1], popt_Hb[2], popt_Hb[3]) + \
             gaussian(x, popt_OIII[0], popt_OIII[1], popt_OIII[2], 0.) + \
             lorentz(x, popt_OIII[3], popt_OIII[4], popt_OIII[5], 0.) + \
             gaussian(x, popt_OIIIb[0], popt_OIIIb[1], popt_OIIIb[2], 0.) + \
             lorentz(x, popt_OIIIb[3], popt_OIIIb[4], popt_OIIIb[5], 0.)

    abs_fit = odata - em_fit

    if plot:
        plt.plot(x, odata, '-k', label="spectra")
        plt.plot(x, bestfit, '--r', label="bestfit absorption line")
        plt.plot(x, abs_fit, '-b', label="absorption spectra - gauss")
        plt.legend()
        plt.show()

    return abs_fit, em_fit


def fit_ems_pvoightcont(galaxy, iline, bestfit, log_bins):
    w=100
    cutout = galaxy[iline - w / 2:iline + w / 2]

    # Find the peak within this cutout, emission line may be shifted from where it is expected to be
    iline2 = np.where(galaxy == np.max(cutout))[0][0]  # This is the index of the real emission line
    x = log_bins[iline2 - w / 2:iline2 + w / 2]
    cutout2 = galaxy[iline2 - w / 2:iline2 + w / 2]

    #plt.plot(log_bins[iline - w / 2:iline + w / 2], cutout, '--b', label="expected location of emission line")
    #plt.plot(x, cutout2, '-r', label="location of emission line")
    #plt.legend()
    #plt.show()

    # To fit the spectra to a pVoight (absorption) and a Guassian (emission) we first want to determine what the
    # optimal parameters for the absorption is by using the best fit absorption template
    bfcutout = bestfit[iline2 - w / 2:iline2 + w / 2]

    b_init = np.mean([np.mean(bfcutout[0:int(0.25*len(bfcutout))]), np.mean(bfcutout[int(0.75*len(bfcutout)):-1])])
    a_init = b_init - np.min(bfcutout)
    x0_init = x[np.argmin(bfcutout)]

    # Fit the absorption spectra to a pVoight function (weighted sum of a Gaussian and Lorentzian function)
    poptbf, pcovbf = curve_fit(pvoight, x, bfcutout, p0=[a_init, 5., x0_init, b_init, 0.8])
    bf_voi = pvoight(x, poptbf[0], poptbf[1], poptbf[2], poptbf[3], poptbf[4])

    #plt.plot(x, bfcutout, '--b', label="bestfit absorption line")
    #plt.plot(x, cutout2, '-r', label="spectra")
    #plt.plot(x, bf_voi, '-g', label="fit of bestfit spectra")
    #plt.legend()
    #plt.show()

    # Now fit the cutout region to a Gaussian function and and the bf_voi pVoight function
    b_em = poptbf[3]
    a_em = np.max(cutout2) + np.abs(poptbf[0]) - b_em
    popt, pcov = curve_fit(fit_emline_over_absline(bf_voi), x, cutout2, p0=[a_em, 2., poptbf[2], a_em, 2., poptbf[2]])

    em_fit = gauss_lorentz(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
    abs_fit = np.subtract(cutout2, em_fit)

    #plt.plot(x, bfcutout, '--b', label="bestfit absorption line")
    #plt.plot(x, cutout2, '-r', label="spectra")
    #plt.plot(x, em_fit, '-k', label="emission spectra")
    #plt.plot(x, abs_fit, '-g', label="absorption spectra")
    #plt.legend()
    #plt.show()

    return popt, pcov


def fit_ems_lincont(galaxy, iline, bestfit, log_bins, bad_pt, mask_region):
    w = 90
    cutout = galaxy[iline - w / 2:iline + w / 2]

    # Find the peak within this cutout, emission line may be shifted from where it is expected to be
    iline2 = np.where(galaxy == np.max(cutout))[0][0]  # This is the index of the real emission line
    wline2 = log_bins[iline2]

    x = log_bins[iline2 - w / 2:iline2 + w / 2]
    cutout2 = galaxy[iline2 - w / 2:iline2 + w / 2]

    b_init = np.mean([np.mean(cutout2[0:int(len(cutout2) / 4)]), np.mean(cutout2[int(3 * len(cutout2) / 4):-1])])
    a_init = np.max(cutout2) - b_init

    bfcutout = bestfit[iline2 - w / 2:iline2 + w / 2]
    poptbf, pcovbf = curve_fit(linear, x, bfcutout, p0=[-0.1, b_init])
    bf_lin = linear(x, poptbf[0], poptbf[1])

    if bad_pt in range(int(x[0]), int(x[-1])):  # apply mask to region with weird bump that bugs fitting method
        x_ma = np.ma.masked_inside(x, mask_region[0], mask_region[1])
        # Now get data only for points that are not masked
        x_ma_data = x[~x_ma.mask]
        cutout2_ma_data = cutout2[~x_ma.mask]
        bf_lin_ma_data = bf_lin[~x_ma.mask]
        popt, pcov = curve_fit(fit_emline_over_cont(bf_lin_ma_data), x_ma_data, cutout2_ma_data,
                               p0=[a_init, 2., wline2, a_init / 2., 2., wline2])

    else:
        popt, pcov = curve_fit(fit_emline_over_cont(bf_lin), x, cutout2,
                               p0=[a_init, 2., wline2, a_init / 2., 2., wline2])

    spec_fit = gauss_lorentz(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) + bf_lin
    em_fit = gauss_lorentz(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
    abs_fit = np.subtract(cutout2, em_fit)

    # plt.plot(x, cutout2, '-k', label="spectra")
    #plt.plot(x, abs_fit, '-r', label="absorption spectra")
    #plt.show()

    return popt, pcov


def pvoight(x, a, w, x0, b, frac):
    return (1 - frac) * gaussian(x, a, w, x0, b) + frac * lorentz(x, a, w, x0, b)


def lorentz(x, a, w, x0, b):
    return a / (1 + ((x - x0) / (w / 2)) ** 2) + b


def gaussian(x, a, w, x0, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * w ** 2)) + b


def linear(x, m, b):
    return m * x + b


def fit_emline_over_cont(bf_lin):
    def emline_over_cont(x, a, w, x0, a2, w2, x02):
        return gauss_lorentz(x, a, w, x0, a2, w2, x02) + bf_lin
    return emline_over_cont


def fit_emline_over_absline(bf_voi):
    def emline_over_absline(x, a, w, x0, a2, w2, x02):
        return gauss_lorentz(x, a, w, x0, a2, w2, x02) + bf_voi
    return emline_over_absline

def gauss_lorentz(x, a, w, x0, a2, w2, x02):
    return gaussian(x, a, w, x0, 0.) + lorentz(x, a2, w2, x02, 0.)


def clean_spec(bin_spec, feat_spec, bad_region):
    """
    Replace noisy region with continuum
    """

    # use splot to determine bad region, the '$' will cahnge scale from wavelength to pixels

    # Use splot to remove the continuum from the spectra - 't' then '-' then 'q' then 'i' and choose output file name
    # To get the continuum (which I will use to fill in the bad region) subtract the continuum removed spectra
    # from the original spectra. i.e. (spectra) - (spectra - continuum) = continuum

    bin_spec_orig = bin_spec.split('bin')[0] + 'orig_bin' + bin_spec.split('bin')[1]  # save original under new name

    assert os.path.exists(
        feat_spec), 'Features fits file has not yet been created, refer to instructions in clean_spec_30'

    if os.path.isfile(bin_spec_orig):
        print('Spectra {} has already been cleaned'.format(bin_spec))
        return

    with fits.open(bin_spec) as spec_hdu:
        spec_data = spec_hdu[0].data
        spec_hdr = spec_hdu[0].header

    with fits.open(bin_spec_features) as feat_hdu:
        feat_data = feat_hdu[0].data

    continuum_data = np.subtract(spec_data, feat_data)

    # now replace bad region in spectra with that of the values of the continuum
    clean_data = spec_data
    clean_data[bad_region[0]:bad_region[1]] = continuum_data[bad_region[0]:bad_region[1]]

    # write into a new fits image the clean data
    clean_hdu = fits.PrimaryHDU()
    clean_hdu.data = clean_data
    clean_hdu.header = spec_hdr
    clean_hdu.writeto(bin_spec, clobber=True)

    orig_hdu = fits.PrimaryHDU()
    orig_hdu.data = spec_data
    orig_hdu.header = spec_hdr
    orig_hdu.writeto(bin_spec_orig, clobber=True)
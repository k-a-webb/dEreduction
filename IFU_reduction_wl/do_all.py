# written by Kristi Webb for data of IC225 taken by Bryan Miller
# similar to an earlier method developed by Gwen Rudie in IDL

from __future__ import print_function
from astropy.io import fits
import numpy as np
from time import clock
from voronoi_2d_binning import voronoi_2d_binning
import glob
import pandas as pd
import os
from scipy import ndimage
from ppxf import ppxf
import ppxf_util as util
import glob
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from sauron_colormap import sauron
from cap_plot_velfield import plot_velfield

"""
This script covers the necessary steps to preform the pPXF kinematic analysis of a given 3D image cube.
In order the steps are:
    - make_table: compile coordinates of signal and noise measurements for the binning
    - voronoi_binning: preform voronoi 2d binning to acheive desired signal to noise
    # - bin_spectra: organise pixel values within their respective bins
    # no longer necessary as not using imcombine anymore
    - combine_spectra: combine spectra of the same bin into fits files
    - ppxf_kinematics: preform pPXF on the binned spectra
    - make_overplot_data: copy the desired spectral range and read mean flux values to be used for plotting
        - continuum: 4370-4870 A
        - OIII: 4980-5035 A (4959-5700)
        - H_beta: 4883-4887 A
        - H_gamma: 4360-4362 A
    - plot_velfield_setup: access arrays to input into cap_plot_velfield()

The next steps would be to make the necessary overplot data with make_overplot_data.cl and plot the results
"""

# These are the files to change based on desired use
DIR_PATH = '/Users/kwebb/IFU_reduction_wl'
PROC_PATH = 'ppxf_proc'  # all output of this script will be in this folder
IMAGE_CUBE = os.path.join(DIR_PATH, 'dcsteqpxbprgN20051205S0006_add.fits')
SCI_EXT = os.path.join(DIR_PATH, 'IC225_2D_james_sci.fits')
VAR_EXT = os.path.join(DIR_PATH, 'IC225_2D_james_var.fits')
TARGET_SN = 30
SCROP_WAVE_RANGE = [4370, 4870]
CONT_FLUX_SCI = 'cont_flux_{}.fits'  # naming convention of continuum wavelength range combined sci spectra from flux*
FLUX_FILE = os.path.join(DIR_PATH, PROC_PATH, 'binned_cont_flux_{}.txt'.format(TARGET_SN))  # make overplot output

# Voronoi binning parameters
VEL_INIT = 1500  # initial guess for velocity
SIG_INIT = 100.  # inital guess of sigma distribution
LAM_RANGE = [4186, 5369]  # wavelength range for logarithmic rebinning
TEMPLATE_PATH = '/Users/kwebb/idl/cappellari/ppxf/spectra/'  # location of template files
# template from http://archives.pd.astro.it/2500-10500/
# adsabs paper: http://adsabs.harvard.edu/abs/2005A%26A...442.1127M
TEMPLATE_FITS = 'Mun1.30z*.fits'  # ****** Make sure to change sigma in pPXF if change template ******
TEMPLATE_RESOLUTION = 2.  # FWHM of the template spectra

# These variables define the chosen file structure I use to organise the output
XYSN_FILE = os.path.join(DIR_PATH, PROC_PATH, 'x_y_signal_noise.txt')  # output of make_table
V2B_FILE = os.path.join(DIR_PATH, PROC_PATH, 'v2b_output_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning
V2B_XY_FILE = os.path.join(DIR_PATH, PROC_PATH, 'v2b_output_xy_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning

DIR_SCI_COMB = os.path.join(DIR_PATH, PROC_PATH, 'comb_fits_sci_{}'.format(TARGET_SN))  # folder with combined spectra
DIR_VAR_COMB = os.path.join(DIR_PATH, PROC_PATH, 'comb_fits_var_{}'.format(TARGET_SN))  # folder with combined var spec
BIN_SCI = 'bin_sci_{}.fits'  # naming convention of combined (sum) sci spectra
IN_FILE_PREFIX = 'bin_sci_*'  # glob prefix for input into pPXF
BIN_VAR = 'bin_var_{}.fits'  # naming convention of combined (sum) var spectra
FLUX_SCI = 'flux_sci_{}.fits'  # naming convention of combined (average) sci spectra
FLUX_VAR = 'flux_var_{}.fits'  # naming convention of combined (average) var spectra

PPXF_PROC_PATH = 'ppxf_output_{}'.format(TARGET_SN)
PPXF_FILE = os.path.join(DIR_PATH, PROC_PATH, 'ppxf_output_sn{}.txt'.format(TARGET_SN))
PPXF_BESTFIT = os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH, 'ppxf_bestfit_{}.fits')

# Various header names
XAXIS_HEADER = 'NAXIS1'  # number of pixels in the x dimension
YAXIS_HEADER = 'NAXIS2'  # number of pixels in the y dimension
CUBE_MATRIX_ELEM_HEADER = 'CD3_3'  # wcs matrix element in dimension 3 3 (ie wavelength dimension)
IMCOMB_MATRIX_ELEM_HEADER = 'CD1_1'  # wcs matrix element in dimension 1 1 (because new array is 1D)
CUBE_REF_PIX_HEADER = 'CRPIX3'  # reference pixel of wavelength dimension
IMCOMB_REF_PIX_HEADER = 'CRPIX1'  # header given to wavelength reference pixel in combined 1D array
CUBE_RA_REF_PIX_HEADER = 'CRVAL3'  # RA at reference pixel of wavelength dimension
IMCOMB_RA_REF_PIX_HEADER = 'CRVAL1'  # header given to RA of wavelength reference pixel in combined 1D array


def make_table():
    """
    Read in the pixel values of the science and varience planes of the image cube and make a table of the
    coordinates with the respective signal and noise measurements
    INPUT: SCI_EXT (2D_IC225_cal_integer_crop_sci.fits), VAR_EXT (2D_IC225_cal_integer_crop_var.fits)
    OUTPUT: XYSN_FILE (x_y_signal_noise.txt)
    """

    with fits.open(SCI_EXT) as sci_hdu:
        sci_data = sci_hdu[0].data
        sci_xaxis = sci_hdu[0].header[XAXIS_HEADER]
        sci_yaxis = sci_hdu[0].header[YAXIS_HEADER]

    with fits.open(VAR_EXT) as var_hdu:
        var_data = var_hdu[0].data
        var_xaxis = var_hdu[0].header[XAXIS_HEADER]
        var_yaxis = var_hdu[0].header[YAXIS_HEADER]

    assert sci_xaxis == var_xaxis, 'Sci and var planes have diffent x dimensions'
    assert sci_yaxis == var_yaxis, 'Sci and var planes have diffent y dimensions'

    with open(XYSN_FILE, 'w') as outfile:
        for i in range(sci_yaxis - 1):
            for j in range(sci_xaxis - 1):
                noise = np.sqrt(var_data[i, j])
                outfile.write('   {}   {}   {}   {}\n'.format(i, j, sci_data[i, j], noise))


def voronoi_binning():
    """
    Follows example script provided from Michele Cappellari for voronoi 2d binning
    INPUT: XYSN_FILE (x_y_signal_noise.txt)
    OUTPUT: V2B_FILE (v2b_output_sn30.txt), V2B_XY_FILE (v2b_output_xy_sn30.txt)

    Output variables of voronoi_2d_binning (copied straight from the description in the script):
        BINNUMBER: Vector (same size as X) containing the bin number assigned to each input pixel. The index goes from
            zero to Nbins-1. This vector alone is enough to make *any* subsequent computation on the binned data.
            Everything else is optional!
        XBIN: Vector (size Nbins) of the X coordinates of the bin generators. These generators uniquely define the
            Voronoi tessellation.
        YBIN: Vector (size Nbins) of Y coordinates of the bin generators.
        XBAR: Vector (size Nbins) of X coordinates of the bins luminosity weighted centroids. Useful for plotting
            interpolated data.
        YBAR: Vector (size Nbins) of Y coordinates of the bins luminosity weighted centroids.
        SN: Vector (size Nbins) with the final SN of each bin.
        NPIXELS: Vector (size Nbins) with the number of pixels of each bin.
        SCALE: Vector (size Nbins) with the scale length of the Weighted Voronoi Tessellation, when the /WVT keyword is
            set. In that case SCALE is *needed* together with the coordinates XBIN and YBIN of the generators, to
            compute the tessellation (but one can also simply use the BINNUMBER vector).
    """

    x, y, signal, noise = np.loadtxt(XYSN_FILE, unpack=True)  # , skiprows=3)

    ''' CHECK VALIDITY OF THIS '''
    noise[noise == 0] = 0.000001

    # Perform the actual computation. The vectors (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
    # are all generated in *output*

    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(x, y, signal, noise, TARGET_SN,
                                                                              plot=1, quiet=0)

    # Save to a text file the initial coordinates of each pixel together with the corresponding bin number computed
    # by this procedure. binNum uniquely specifies the bins and for this reason it is the only
    # number required for any subsequent calculation on the bins.

    np.savetxt(V2B_FILE, np.column_stack([x, y, binNum]), header='x  y  binNum', fmt=b'%10.6f %10.6f %8i')
    np.savetxt(V2B_XY_FILE, np.column_stack([xBar, yBar, xNode, yNode]), header='xBar  yBar  xNode   yNode',
               fmt=b'%10.6f %10.6f %10.6f %10.6f')


def bin_spectra():
    """
    THIS STEP IS NO LONGER NECESSARY AS I NO LONGER USE IMCOMBINE

    Organise the spectra according to the voronoi binning output into list of pixels with the same bin
    INPUT: V2B_FILE (v2b_output_sn30.txt)
    OUTPUT: DIR_SCI (comb_lists_sci_{S/N}/imcomb_sci_{S/N}.lis), DIR_VAR (comb_lists_var_{S/N}/imcomb_var_{S/N}.lis)
    """

    if not os.path.exists(DIR_SCI):
        os.mkdir(DIR_SCI)
    if not os.path.exists(DIR_VAR):
        os.mkdir(DIR_VAR)

    v2b_output = pd.read_table(V2B_FILE, sep=r"\s*", engine='python', skiprows=1, names=["x", "y", "bin"])
    print('Highest bin number: ', np.amax(v2b_output['bin'].values))

    for i in range(np.amax(v2b_output['bin'].values) + 1):
        line = v2b_output.query('bin == {}'.format(i))

        file_sci_i = FILE_SCI.format(i)
        file_var_i = FILE_VAR.format(i)

        with open('{}/{}'.format(DIR_SCI, file_sci_i), 'w') as outfile:
            for j in range(len(line)):
                # if int(line["x"].values[j - 1]) != 0:
                outfile.write('{}[sci][{},{},*]\n'.format(IMAGE_CUBE, int(line["x"].values[j - 1]) + 1,
                                                          int(line["y"].values[j - 1]) + 1))

        with open('{}/{}'.format(DIR_VAR, file_var_i), 'w') as outfile:
            for j in range(len(line)):
                # if int(line["x"].values[j - 1]) != 0:
                outfile.write('{}[var][{},{},*]\n'.format(IMAGE_CUBE, int(line["x"].values[j - 1] + 1),
                                                          int(line["y"].values[j - 1]) + 1))


def combine_spectra():
    """
    Combine each pixel of the same bin (according to the lists output by bin_spectra) into a single fits file
    INPUT: V2B_FILE (v2b_output_sn30.txt)
    OUTPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), DIR_VAR_COMB (comb_fis_var_{S/N}/bin_var_{S/N}.fits)
    """

    # make output folders if they don't exist already
    if not os.path.exists(DIR_SCI_COMB):
        os.mkdir(DIR_SCI_COMB)
    if not os.path.exists(DIR_VAR_COMB):
        os.mkdir(DIR_VAR_COMB)

    # read in the output of the voronoi binning
    v2b_output = pd.read_table(V2B_FILE, sep=r"\s*", engine='python', skiprows=1, names=["x", "y", "bin"])

    # for each bin, combine the spectra of the same bin
    for i in range(np.amax(v2b_output['bin'].values) + 1):
        bin_sci_i = os.path.join(DIR_SCI_COMB, BIN_SCI.format(i))
        bin_var_i = os.path.join(DIR_VAR_COMB, BIN_VAR.format(i))
        flux_sci_i = os.path.join(DIR_SCI_COMB, FLUX_SCI.format(i))
        flux_var_i = os.path.join(DIR_VAR_COMB, FLUX_VAR.format(i))

        # Check to see how many pixels are in this bin
        spaxel_list = v2b_output.query('bin == {}'.format(i))
        print('Number of pixels in bin {}: {}'.format(i, len(spaxel_list)))

        '''
        file_sci_i = os.path.join(DIR_SCI, FILE_SCI.format(i))
        file_var_i = os.path.join(DIR_VAR, FILE_VAR.format(i))
        assert os.path.exists(file_sci_i), "File {} does not exist".format(file_sci_i)
        assert os.path.exists(file_var_i), "File {} does not exist".format(file_var_i)

        import pyraf.iraf as iraf
        iraf.imcombine("@{}".format(file_sci_i), "{}".format(bin_sci_i), combine="sum")
        iraf.imcombine("@{}".format(file_var_i), "{}".format(bin_var_i), combine="sum")
        iraf.imcombine("@{}".format(file_sci_i), "{}".format(flux_sci_i), combine="average")
        iraf.imcombine("@{}".format(file_var_i), "{}".format(flux_var_i), combine="average")
        '''

        imcombine(bin_sci_i, bin_var_i, spaxel_list, 'sum')
        imcombine(flux_sci_i, flux_var_i, spaxel_list, 'average')


def imcombine(outsci, outvar, spaxel_list, combine_method):
    """
    In order to avoid a 'floating point error' in imexamine I'm going to try and combine the fits sections as
    numpy arrays instead and manually copy the header values. Only works for combine='average'|'sum' at the moment
    """

    assert combine_method is 'sum' or 'average', 'Combine is not "sum" or "average"'

    # IMAGE_CUBE[int(line["x"].values[j - 1]) + 1, int(line["y"].values[j - 1]) + 1, *]

    with fits.open(IMAGE_CUBE) as image_cube_hdu:
        # Filename: dcsteqpxbprgN20051205S0006_add.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU     214   ()
        # 1    SCI         ImageHDU        68   (76, 49, 1300)   float32
        # 2    VAR         ImageHDU        68   (76, 49, 1300)   float32
        image_cube_sci = image_cube_hdu[1].data
        image_cube_var = image_cube_hdu[2].data
        cube_header = image_cube_hdu['SCI'].header

    imcomb_sci = []
    imcomb_var = []

    # image_cube_sci.shape = (1300, 49, 76) ie (z,y,x)
    for k in range(image_cube_sci.shape[0]):
        pix_sci = []
        pix_var = []
        for j in range(len(spaxel_list)):
            pix_sci.append(image_cube_sci[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1])
            pix_var.append(image_cube_var[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1])

        if combine_method == 'sum':
            imcomb_sci.append(np.array(pix_sci).sum())
            imcomb_var.append(np.array(pix_var).sum())
        elif combine_method == 'average':
            imcomb_sci.append(np.array(pix_sci).mean())
            imcomb_var.append(np.array(pix_var).mean())

    write_imcomb_fits(imcomb_sci, outsci, cube_header)
    write_imcomb_fits(imcomb_var, outvar, cube_header)


def write_imcomb_fits(outdata, outfile, cube_header):
    """
    Write the combined spectra to a fits file with headers from the wavelength dimension of the 3D cube
    """
    outfile_hdu = fits.PrimaryHDU()
    outfile_hdu.header['CD1_1'] = cube_header['CD3_3']
    outfile_hdu.header['CRPIX1'] = cube_header['CRPIX3']
    outfile_hdu.header['CRVAL1'] = cube_header['CRVAL3']
    outfile_hdu.data = outdata
    outfile_hdu.writeto(outfile, clobber=True)


def ppxf_kinematics():
    """
    Follow the pPXF usage example by Michile Cappellari
    INPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), TEMPLATE_* (spectra/Mun1.30z*.fits)
    OUTPUT: PPXF_FILE (ppxf_output_sn30.txt), OUT_BESTFIT_FITS (ppxf_fits/ppxf_bestfit_sn30.fits)
    """

    # Read a galaxy spectrum and define the wavelength range
    in_file = glob.glob(os.path.join(DIR_SCI_COMB, IN_FILE_PREFIX))

    hdu = fits.open(in_file[0])
    gal_lin = hdu[0].data
    h1 = hdu[0].header

    # lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
    # FWHM_gal = 4.2 # SAURON has an instrumental resolution FWHM of 4.2A.

    FWHM_gal = 2.3  # GMOS IFU has an instrumental resolution FWHM of 2.3 A

    # If the galaxy is at a significant redshift (z > 0.03), one would need to apply a large velocity shift in PPXF to
    # match the template to the galaxy spectrum.This would require a large initial value for the velocity (V > 1e4 km/s)
    # in the input parameter START = [V,sig]. This can cause PPXF to stop! The solution consists of bringing the galaxy
    # spectrum roughly to the rest-frame wavelength, before calling PPXF. In practice there is no need to modify the
    # spectrum before the usual LOG_REBIN, given that a red shift corresponds to a linear shift of the log-rebinned
    # spectrum. One just needs to compute the wavelength range in the rest-frame and adjust the instrumental resolution
    # of the galaxy observations. This is done with the following three commented lines:
    #
    # z = 1.23 # Initial estimate of the galaxy redshift
    # lamRange1 = lamRange1/(1+z) # Compute approximate restframe wavelength range
    # FWHM_gal = FWHM_gal/(1+z)   # Adjust resolution in Angstrom

    # galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_lin)
    # galaxy = galaxy/np.median(galaxy) # Normalize spectrum to avoid numerical issues
    # noise = galaxy*0 + 0.0049           # Assume constant noise per pixel here

    galaxy, logLam1, velscale = util.log_rebin(LAM_RANGE, gal_lin)
    galaxy = galaxy / np.median(galaxy)  # Normalize spectrum to avoid numerical issues

    # Read the list of filenames from the Single Stellar Population library by Vazdekis (1999, ApJ, 513, 224). A subset
    # of the library is included for this example with permission. See http://purl.org/cappellari/software
    # for suggestions of more up-to-date stellar libraries.

    # vazdekis = glob.glob(dir + 'Rbi1.30z*.fits')
    # vazdekis.sort()
    # FWHM_tem = 1.8 # Vazdekis spectra have a resolution FWHM of 1.8A.

    Mun = glob.glob(TEMPLATE_PATH + TEMPLATE_FITS)
    FWHM_tem = TEMPLATE_RESOLUTION

    # Extract the wavelength range and logarithmically rebin one spectrum to the same velocity scale of the SAURON
    # galaxy spectrum, to determine the size needed for the array which will contain the template spectra.

    # hdu = fits.open(vazdekis[0])
    # ssp = hdu[0].data
    # h2 = hdu[0].header
    # lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    # sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    # templates = np.empty((sspNew.size,len(vazdekis)))

    hdu = fits.open(Mun[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0., h2['CDELT1'] * (h2['NAXIS1'] - 1)])
    sspNew, logLam2, velscale2 = util.log_rebin(lamRange2, ssp)
    templates = np.empty((sspNew.size, len(Mun)))

    # Convolve the whole Vazdekis library of spectral templates with the quadratic difference between the SAURON and the
    # Vazdekis instrumental resolution. Logarithmically rebin and store each template as a column in the array TEMPLATES.

    # Quadratic sigma difference in pixels Vazdekis --> SAURON
    # The formula below is rigorously valid if the shapes of the instrumental spectral profiles are well approximated
    # by Gaussians.

    # FROM GWENS IDL SCRIPT ----------------------------------------------------------
    # Example: BC03 spectra have a resolution of 3A at 5100A, this corresponds to sigma ~ 3.0/5100*3e5/2.35 = 75 km/s.
    # The GMOS IFU has an instrumental resolution of 2.3 A, this corresponds to sigma ~ 2.3/5100*3e5/2.35 = 56.8 km/s.
    # The quadratic difference is sigma = sqrt(56.8^2 - 75^2) = undefined
    # (the above reasoning can be applied if the shape of the instrumental spectral profiles can be well approximated
    # by a Gaussian).
    # For the lower resolution models, we must degrade the DATA to fit the models.
    # Thus: The quadratic difference is sigma = sqrt(75^2 - 56.8^2) = 49.0 km/s
    # ---------------------------------------------------------------------------------

    FWHM_dif = np.sqrt(FWHM_gal ** 2 - FWHM_tem ** 2)
    sigma = FWHM_dif / 2.355 / h2['CDELT1']  # Sigma difference in pixels

    # Logarithmically rebin the whole Mun library of spectra, and store each template as a column in the array TEMPLATES,

    # for j in range(len(vazdekis)):
    # hdu = fits.open(vazdekis[j])
    # ssp = hdu[0].data
    # ssp = ndimage.gaussian_filter1d(ssp, sigma)
    # sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    #        templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

    for j in range(len(Mun)):
        hdu = fits.open(Mun[j])
        ssp = hdu[0].data
        ssp = ndimage.gaussian_filter1d(ssp, sigma)
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale2)
        templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

    # The galaxy and the template spectra do not have the same starting wavelength. For this reason an extra velocity
    # shift DV has to be applied to the template to fit the galaxy spectrum. We remove this artificial shift by using
    # the keyword VSYST in the call to PPXF below, so that all velocities are measured with respect to DV. This assume
    # the redshift is negligible.In the case of a high-redshift galaxy one should de-redshift its wavelength to the
    # rest frame before using the line below (see above).

    vel_list = []
    sig_list = []
    h3_list = []
    h4_list = []
    dV_list = []
    dsigma_list = []
    dh3_list = []
    dh4_list = []

    for j in range(len(in_file)):
        print('>>>>> {}  {}'.format(j, in_file[j]))

        hdu = fits.open(in_file[j])
        b_gal = hdu[0].data
        b_gal = ndimage.gaussian_filter1d(b_gal, sigma)

        galaxy, logLam1, velscale = util.log_rebin(LAM_RANGE, b_gal, velscale=velscale)
        #   error = galaxy*0 + 1 # Assume constant error

        ''' Replace 0 values with interpretad values, check validity of this '''
        mask = np.isnan(galaxy)
        galaxy[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), galaxy[~mask])

        noise = galaxy * 0 + 0.0049  # Assume constant noise per pixel here

        c = 299792.458
        dv = (logLam2[0] - logLam1[0]) * c  # km/s

        #   vel = 1500.  # Initial estimate of the galaxy velocity in km/s
        z = np.exp(VEL_INIT / c) - 1  # Relation between velocity and redshift in pPXF

        goodPixels = util.determine_goodpixels(logLam1, lamRange2, z)
        # Gwen uses goodpixel range: [range(2,195),range(235,790),range(830,890),range(930,940),range(980,1300)]

        # Here the actual fit starts. The best fit is plotted on the screen.
        # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.

        start = [VEL_INIT, SIG_INIT]  # (km/s), starting guess for [V,sigma]
        t = clock()

        pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodPixels, plot=True, moments=4,
                  degree=4, vsyst=dv)

        print("Formal errors:")
        print("     dV    dsigma   dh3      dh4")
        print("".join("%8.2g" % f for f in pp.error * np.sqrt(pp.chi2)))

        print('Elapsed time in PPXF: %.2f s' % (clock() - t))

        # If the galaxy is at significant redshift z and the wavelength has been de-redshifted with the three lines
        # "z = 1.23..." near the beginning of this procedure, the best-fitting redshift is now given by the following
        # commented line (equation 2 of Cappellari et al. 2009, ApJ, 704, L34;
        # http://adsabs.harvard.edu/abs/2009ApJ...704L..34C)

        #   print, 'Best-fitting redshift z:', (z + 1)*(1 + sol[0]/c) - 1

        # Gwen obtains the velocity and sigma information from the SOL parameter
        vel_list.append(pp.sol[0])
        sig_list.append(pp.sol[1])
        h3_list.append(pp.sol[2])
        h4_list.append(pp.sol[3])
        dV_list.append((pp.error * np.sqrt(pp.chi2))[0])
        dsigma_list.append((pp.error * np.sqrt(pp.chi2))[1])
        dh3_list.append((pp.error * np.sqrt(pp.chi2))[2])
        dh4_list.append((pp.error * np.sqrt(pp.chi2))[3])

        hdu_best = fits.PrimaryHDU()
        hdu_best.data = pp.bestfit
        hdu_best.writeto(PPXF_BESTFIT.format(j), clobber=True)

    np.savetxt(PPXF_FILE,
               np.column_stack([vel_list, sig_list, h3_list, h4_list, dV_list, dsigma_list, dh3_list, dh4_list]),
               fmt=b'%10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f',
               header='velocity  sigma  h3  h4  dV  dsigma  dh3  dh4')


def make_overplot_data():
    """
    NOT YET WORKING, REFER TO make_overplot_data.cl FOR A WORKING COPY
    """

    from pyraf import iraf

    iraf.noao()
    iraf.onedspec()

    files_in_dir = glob.glob(os.path.join(DIR_SCI_COMB, FLUX_SCI.format('*')))
    assert len(files_in_dir) > 0, 'No files match {}'.format(os.path.join(DIR_SCI_COMB, FLUX_SCI.format('*')))

    mean_flux = []

    for i in range(len(files_in_dir)):

        cont_flux_sci_i = os.path.join(DIR_SCI_COMB, CONT_FLUX_SCI.format(i))

        if not os.path.exists(cont_flux_sci_i):
            iraf.scopy(files_in_dir[i], cont_flux_sci_i, w1=SCROP_WAVE_RANGE[0], w2=SCROP_WAVE_RANGE[1])
        cont_flux_data = fits.getdata(cont_flux_sci_i, 0)
        assert cont_flux_data.ndim != 0, "Scrop'd array is empty"
        mean_flux.append(cont_flux_data.mean())

    np.array(mean_flux).tofile(FLUX_FILE, sep='\n')


def plot_velfield_setup():
    """
    Plot velfield as done in Gwen's plot_ppxf_30.pro
    """

    xbar, ybar, xnode, ynode = np.loadtxt(V2B_XY_FILE, unpack=True, skiprows=1)
    vel, sig, h3, h4, dvel, dsig, dh3, dh4 = np.loadtxt(PPXF_FILE, unpack=True)
    flux = np.loadtxt(FLUX_FILE, unpack=True)

    assert len(xbar) == len(ybar), 'Xbar is not the same length as Ybar'
    assert len(xbar) == len(vel), 'Xbar is not the same length as vel'
    assert len(xbar) == len(flux), 'Xbar is not the same length as flux'

    plt.clf()
    plt.title('Velocity')
    plot_velfield(xbar, ybar, vel, flux=flux, colorbar=True, label='km/s')
    plt.show()


if __name__ == '__main__':
    """
    Follow the steps described above
    """

    if not os.path.exists(os.path.join(DIR_PATH, PROC_PATH)):
        os.mkdir(os.path.join(DIR_PATH, PROC_PATH))

    if not os.path.exists(XYSN_FILE):
        print('>>>>> Making table')
        make_table()

    if not os.path.exists(V2B_FILE):
        print('>>>>> Voronoi binning')
        t = clock()
        voronoi_binning()
        print('Elapsed time: %.2f seconds' % (clock() - t))

    # if not os.path.exists(os.path.join(DIR_SCI, FILE_VAR.format(1))):
    # print('>>>>> Binning spectra')
    #     bin_spectra()

    if not os.path.exists(os.path.join(DIR_VAR_COMB, BIN_VAR.format('1'))):
        print('>>>>> Combining spectra')
        combine_spectra()

    if not os.path.exists(PPXF_FILE):
        print('>>>>> pPXF')
        if not os.path.exists(os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH)):
            os.mkdir(os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH))
        ppxf_kinematics()

    if not os.path.exists(FLUX_FILE):
        print('>>>>> overplot data')
        make_overplot_data()

    print('>>>>> Plotting vel field')
    plot_velfield_setup()

    print('>>>>> Done')

""" At S/N 30 with James' reduction and cubing method
41  initial bins.
Reassign bad bins...
26  good bins.
Modified Lloyd algorithm...
Iter:    1  Diff: 139.2
Iter:    2  Diff: 13.32
Iter:    3  Diff: 6.005
Iter:    4  Diff: 3.506
Iter:    5  Diff: 1.943
Iter:    6  Diff: 1.511
Iter:    7  Diff: 0.9909
Iter:    8  Diff: 0.4755
Iter:    9  Diff: 0.5045
Iter:   10  Diff: 0.2212
Iter:   11  Diff: 0.1989
Iter:   12  Diff: 0.06859
Iter:   13  Diff: 0.1068
Iter:   14  Diff: 0.1636
Iter:   15  Diff: 0.08666
Iter:   16  Diff: 0.07709
Iter:   17  Diff: 0.1417
Iter:   18  Diff: 0.06984
Iter:   19  Diff: 0.0317
Iter:   20  Diff: 0.009915
Iter:   21  Diff: 0
20  iterations.
Unbinned pixels:  0  /  3600
Fractional S/N scatter (%): 13.3416729619
"""
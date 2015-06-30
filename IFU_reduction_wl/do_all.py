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

"""
This script covers the necessary steps to preform the pPXF kinematic analysis of a given 3D image cube.
In order the steps are:
    - make_table: compile coordinates of signal and noise measurements for the binning
    - voronoi_binning: preform voronoi 2d binning to acheive desired signal to noise
    - bin_spectra: organise pixel values within their respective bins
    - combine_spectra: combine spectra of the same bin into fits files
    - ppxf_kinematics: preform pPXF on the binned spectra

The next steps would be to make the necessary overplot data with make_overplot_data.cl and plot the results
"""

# These are the files to change based on desired use
DIR_PATH = '/Users/kwebb/IFU_reduction_wl'
IMAGE_CUBE = os.path.join(DIR_PATH, 'dcsteqpxbprgN20051205S0006_add.fits')
SCI_EXT = os.path.join(DIR_PATH, 'IC225_2D_james_sci.fits')
VAR_EXT = os.path.join(DIR_PATH, 'IC225_2D_james_var.fits')
TARGET_SN = 30
VEL_INIT = 1500
SIG_INIT = 100.
XAXIS_HEADER = 'NAXIS1'
YAXIS_HEADER = 'NAXIS2'
XDIM_HEADER = 'CD1_1'
YDIM_HEADER = 'CD2_2'
LAM_RANGE = [4186, 5369]
TEMPLATE_PATH = '/Users/kwebb/idl/cappellari/ppxf/spectra/'
TEMPLATE_FITS = 'Mun1.30z*.fits'  # Make sure to change sigma in pPXF if change template

# These variables define the chosen file structure I use to organise the output
PROC_PATH = 'ppxf_proc'
if not os.path.exists(os.path.join(DIR_PATH, PROC_PATH)):
    os.mkdir(os.path.join(DIR_PATH, PROC_PATH))
XYSN_FILE = os.path.join(DIR_PATH, PROC_PATH, 'x_y_signal_noise.txt')
V2B_FILE = os.path.join(DIR_PATH, PROC_PATH, 'v2b_output_sn30.txt')
V2B_XY_FILE = os.path.join(DIR_PATH, PROC_PATH, 'v2b_output_sn30_xy.txt')
DIR_SCI = os.path.join(DIR_PATH, PROC_PATH, 'comb_lists_sci_{}'.format(TARGET_SN))
DIR_VAR = os.path.join(DIR_PATH, PROC_PATH, 'comb_lists_var_{}'.format(TARGET_SN))
DIR_SCI_COMB = os.path.join(DIR_PATH, PROC_PATH, 'comb_fits_sci_{}'.format(TARGET_SN))
DIR_VAR_COMB = os.path.join(DIR_PATH, PROC_PATH, 'comb_fits_var_{}'.format(TARGET_SN))
FILE_SCI = 'imcomb_sci_{}.lis'
FILE_VAR = 'imcomb_var_{}.lis'
BIN_SCI = 'bin_sci_{}.fits'
BIN_VAR = 'bin_var_{}.fits'
IN_FILE_PREFIX = 'bin_sci_*'
PPXF_PROC_PATH = 'ppxf_output'
OUT_VEL_SIGMA = os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH, 'vel_sigma_output_sn30.txt')
OUT_VEL = os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH, 'vel_output_sn30.txt')
OUT_BESTFIT_FITS = os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH, 'ppxf_bestfit_sn30.fits')
OUT_ERROR_FITS = os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH, 'ppxf_error_sn30.fits')
DIR_SPEC = os.path.join(DIR_PATH, PROC_PATH, 'bin_spec_sn30')

MIN_BINS = 23
NUM_SMALL_FITS = 6


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
        sci_xdim = sci_hdu[0].header[XDIM_HEADER]
        sci_ydim = sci_hdu[0].header[XDIM_HEADER]

    with fits.open(VAR_EXT) as var_hdu:
        var_data = var_hdu[0].data
        var_xaxis = var_hdu[0].header[XAXIS_HEADER]
        var_yaxis = var_hdu[0].header[YAXIS_HEADER]
        var_xdim = var_hdu[0].header[XDIM_HEADER]
        var_ydim = var_hdu[0].header[XDIM_HEADER]

    assert sci_xaxis == var_xaxis, 'Sci and var planes have diffent x dimensions'
    assert sci_yaxis == var_yaxis, 'Sci and var planes have diffent y dimensions'

    with open(XYSN_FILE, 'w') as outfile:
        for i in range(sci_yaxis - 1):
            for j in range(sci_xaxis - 1):
                # noise = np.sqrt(var_data[i, j] * var_xdim * var_ydim)
                # outfile.write('   {}   {}   {}   {}\n'.format(i, j, sci_data[i, j] * sci_xdim * sci_ydim, noise))
                noise = np.sqrt(var_data[i, j])
                outfile.write('   {}   {}   {}   {}\n'.format(i, j, sci_data[i, j], noise))


def voronoi_binning():
    """
    Follows example script provided from Michele Cappellari for voronoi 2d binning
    INPUT: XYSN_FILE (x_y_signal_noise.txt)
    OUTPUT: V2B_FILE (v2b_output_sn30.txt), V2B_XY_FILE (v2b_output_sn30_xy.txt)
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

    # assert np.amax(np.transpose(binNum)) > MIN_BINS, 'Too few bins ({}) of expected at least {}'.format(np.amax(
    # np.transpose(binNum)), MIN_BINS)

    np.savetxt(V2B_FILE, np.column_stack([x, y, binNum]), fmt=b'%10.6f %10.6f %8i')
    np.savetxt(V2B_XY_FILE, np.column_stack([xBar, yBar, xNode, yNode]), fmt=b'%10.6f %10.6f %10.6f %10.6f')


def bin_spectra():
    """
    Organise the spectra according to the voronoi binning output into list of pixels with the same bin
    INPUT: V2B_FILE (v2b_output_sn30.txt)
    OUTPUT: DIR_SCI (comb_lists_sci_{S/N}/imcomb_sci_{S/N}.lis), DIR_VAR (comb_lists_var_{S/N}/imcomb_var_{S/N}.lis)
    """

    if not os.path.exists(DIR_SCI):
        os.mkdir(DIR_SCI)
    if not os.path.exists(DIR_VAR):
        os.mkdir(DIR_VAR)

    v2b_output = pd.read_table(V2B_FILE, sep=r"\s*", engine='python', names=["x", "y", "bin"])  # , skiprows=1)
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
    INPUT: DIR_SCI (comb_lists_sci_{S/N}/imcomb_sci_{S/N}.lis), DIR_VAR (comb_lists_var_{S/N}/imcomb_var_{S/N}.lis)
    OUTPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), DIR_VAR_COMB (comb_fis_var_{S/N}/bin_var_{S/N}.fits)
    """

    if not os.path.exists(DIR_SCI_COMB):
        os.mkdir(DIR_SCI_COMB)
    if not os.path.exists(DIR_VAR_COMB):
        os.mkdir(DIR_VAR_COMB)

    import pyraf.iraf as iraf

    v2b_output = pd.read_table(V2B_FILE, sep=r"\s*", engine='python', names=["x", "y", "bin"])  # , skiprows=1)

    for i in range(12, 13):  # , np.amax(v2b_output['bin'].values) + 1):
        file_sci_i = FILE_SCI.format(i)
        file_var_i = FILE_VAR.format(i)
        bin_sci_i = BIN_SCI.format(i)
        bin_var_i = BIN_VAR.format(i)

        assert os.path.exists(os.path.join(DIR_SCI, file_sci_i)), "File {}/{} does not exist".format(DIR_SCI,
                                                                                                     file_sci_i)
        assert os.path.exists(os.path.join(DIR_VAR, file_var_i)), "File {}/{} does not exist".format(DIR_VAR,
                                                                                                     file_var_i)
        # Check to see how many pixels are in this bin
        spaxel_list = v2b_output.query('bin == {}'.format(i))
        print('\nNumber of pixels in bin {}: {}'.format(i, len(spaxel_list)))

        try:
            iraf.imcombine("@{}/{}".format(DIR_SCI, file_sci_i), "{}/{}".format(DIR_SCI_COMB, bin_sci_i), combine="sum")
            iraf.imcombine("@{}/{}".format(DIR_VAR, file_var_i), "{}/{}".format(DIR_VAR_COMB, bin_var_i), combine="sum")
        except Exception, e:
            print (e)
            print('>>>>> WARNING: DOING MATH MYSELF <<<<<')
            avoid_imexam_err(i, spaxel_list)


def avoid_imexam_err(i, spaxel_list):
    """
    In order to avoid a 'floating point error' in imexamine I'm going to try and combine the fits sections as
    numpy arrays
    """

    bin_sci_i = BIN_SCI.format(i)
    bin_var_i = BIN_VAR.format(i)

    # IMAGE_CUBE[int(line["x"].values[j - 1]) + 1, int(line["y"].values[j - 1]) + 1, *]

    with fits.open(IMAGE_CUBE) as image_cube_hdu:
        # Filename: dcsteqpxbprgN20051205S0006_add.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU     214   ()
        # 1    SCI         ImageHDU        68   (76, 49, 1300)   float32
        # 2    VAR         ImageHDU        68   (76, 49, 1300)   float32
        image_cube_sci = image_cube_hdu[1].data
        image_cube_var = image_cube_hdu[2].data

    # image_cube_sci_zero = np.zeros(image_cube_sci.shape)
    # image_cube_var_zero = np.zeros(image_cube_var.shape)

    imcomb_sci = []
    imcomb_var = []

    print(image_cube_sci.shape)
    print(len(spaxel_list))

    # image_cube_sci.shape = (1300, 49, 76) ie (z,y,x)
    for k in range(image_cube_sci.shape[0]):
        sum_sci = 0
        sum_var = 0
        for j in range(len(spaxel_list)):
            sum_sci += image_cube_sci[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1]
            sum_var += image_cube_var[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1]
        imcomb_sci.append(sum_sci)
        imcomb_var.append(sum_var)

    image_cube_sci_hdu = fits.PrimaryHDU()
    image_cube_sci_hdu.data = imcomb_sci
    image_cube_sci_hdu.writeto("{}/{}".format(DIR_SCI_COMB, bin_sci_i), clobber=True)

    image_cube_var_hdu = fits.PrimaryHDU()
    image_cube_var_hdu.data = imcomb_var
    image_cube_var_hdu.writeto("{}/{}".format(DIR_VAR_COMB, bin_var_i), clobber=True)


def avoid_imexam_err2(i, spaxel_list):
    """
    In order to avoid a 'floating point error' in imexamine if you attempt to add too many files at once I split the
    group into four parts, then add the four parts at the end.
    """

    import pyraf.iraf as iraf

    file_sci_i = FILE_SCI.format(i)
    file_var_i = FILE_VAR.format(i)
    bin_sci_i = BIN_SCI.format(i)
    bin_var_i = BIN_VAR.format(i)

    spaxels = np.genfromtxt("{}/{}".format(DIR_SCI, file_sci_i), dtype='str')
    spaxels_var = np.genfromtxt("{}/{}".format(DIR_VAR, file_var_i), dtype='str')

    spaxels_files = []
    spaxels_var_files = []

    for j in range(0, NUM_SMALL_FITS):
        (spaxels[j:(j + 1) * len(spaxel_list) / 4]).tofile(
            "{}/tmp_{}".format(DIR_SCI, FILE_SCI.format('{}_{}'.format(i, j + 1))), sep='\n')
        (spaxels_var[j:(j + 1) * len(spaxel_list) / 4]).tofile(
            "{}/tmp_{}".format(DIR_VAR, FILE_VAR.format('{}_{}'.format(i, j + 1))), sep='\n')

        spaxels_files.append("{}/tmp_{}[0]".format(DIR_SCI_COMB, BIN_SCI.format('{}_{}'.format(i, j + 1))))
        spaxels_var_files.append("{}/tmp_{}[0]".format(DIR_VAR_COMB, BIN_VAR.format('{}_{}'.format(i, j + 1))))

    (np.array(spaxels_files)).tofile("{}/tmp_{}".format(DIR_SCI, FILE_SCI.format('{}_list'.format(i))), sep='\n')
    (np.array(spaxels_var_files)).tofile("{}/tmp_{}".format(DIR_VAR, FILE_VAR.format('{}_list'.format(i))), sep='\n')

    for j in range(0, NUM_SMALL_FITS):
        iraf.imcombine("@{}/tmp_{}".format(DIR_SCI, FILE_SCI.format('{}_{}'.format(i, j + 1))),
                       "{}/tmp_{}".format(DIR_SCI_COMB, BIN_SCI.format('{}_{}'.format(i, j + 1))), combine="sum")
        iraf.imcombine("@{}/tmp_{}".format(DIR_VAR, FILE_VAR.format('{}_{}'.format(i, j + 1))),
                       "{}/tmp_{}".format(DIR_VAR_COMB, BIN_VAR.format('{}_{}'.format(i, j + 1))), combine="sum")

    iraf.imcombine("@{}/tmp_{}".format(DIR_SCI, FILE_SCI.format('{}_list'.format(i))),
                   "{}/{}".format(DIR_SCI_COMB, bin_sci_i), combine="sum")
    iraf.imcombine("@{}/tmp_{}".format(DIR_VAR, FILE_VAR.format('{}_list'.format(i))),
                   "{}/{}".format(DIR_VAR_COMB, bin_var_i), combine="sum")

    for tmpfile in glob.glob('{}/tmp*'.format(DIR_SCI)):
        os.remove(tmpfile)
    for tmpfile in glob.glob('{}/tmp*'.format(DIR_VAR)):
        os.remove(tmpfile)
    for tmpfile in glob.glob('{}/tmp*'.format(DIR_SCI_COMB)):
        os.remove(tmpfile)
    for tmpfile in glob.glob('{}/tmp*'.format(DIR_VAR_COMB)):
        os.remove(tmpfile)


def ppxf_kinematics():
    """
    Follow the pPXF usage example by Michile Cappellari
    INPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), TEMPLATE_* (spectra/Mun1.30z*.fits)
    OUTPUT: OUT_VEL_SIGMA (vel_sigma_output_sn30.txt), OUT_VEL (vel_output_sn30.txt),
            OUT_BESTFIT_FITS (ppxf_bestfit_sn30.fits), OUT_ERROR_FITS(ppxf_error_sn30.fits)
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

    # FWHM_dif = np.sqrt(FWHM_gal ** 2 - FWHM_tem ** 2)
    # sigma = FWHM_dif / 2.355 / h2['CDELT1']  # Sigma difference in pixels
    sigma = 49.

    # Logarithmically rebin the whole Mun library of spectra, and store each template as a column in the array TEMPLATES,

    # for j in range(len(vazdekis)):
    # hdu = fits.open(vazdekis[j])
    #        ssp = hdu[0].data
    #        ssp = ndimage.gaussian_filter1d(ssp, sigma)
    #        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
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
        # Because MOMENTS=2 the SOL parameter contains [velocity, sigma]
        vel_list.append(pp.sol[0])
        sig_list.append(pp.sol[1])

    if not os.path.exists(os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH)):
        os.mkdir(os.path.join(DIR_PATH, PROC_PATH, PPXF_PROC_PATH))

    np.savetxt(OUT_VEL_SIGMA, np.column_stack([range(len(vel_list)), vel_list, sig_list]), fmt=b'%8i   %10.6f   %10.6f',
               header='   index   velocity   sigma')
    np.savetxt(OUT_VEL, np.column_stack([vel_list]), fmt=b'%10.6f')

    hdu_best = fits.PrimaryHDU()
    hdu_best.data = pp.bestfit
    hdu_best.writeto(OUT_BESTFIT_FITS)

    hdu_error = fits.PrimaryHDU()
    hdu_error.data = pp.error
    hdu_error.writeto(OUT_ERROR_FITS)
    # ---------------------------------------------------------------------------------------------


def make_overplot_data():
    """
    NOT YET WORKING, REFER TO make_overplot_data.cl FOR A WORKING COPY
    """

    if not os.path.exists(DIR_SPEC):
        os.mkdir(DIR_SPEC)

    import pyraf.iraf as iraf

    files_in_dir = glob.glob(os.path.join(DIR_SCI, '*.lis'))
    num_files = len(files_in_dir)

    for i in range(len(files_in_dir)):
        file_sci_i = FILE_SCI.format(i)
        flux_sci_i = 'flux_{}.lis'.format(i)
        cont_flux_sci_i = 'cont_flux_{}.fits'.format(i)

        iraf.imcombine("@{}/{}".format(DIR_SCI, file_sci_i), "{}/{}".format(DIR_SPEC, flux_sci_i), combine="average")
        iraf.scopy("{}/{}".format(DIR_SPEC, flux_sci_i), "{}/{}".format(DIR_SPEC, cont_flux_sci_i), w1=4370, w2=4870)
    '''
    iraf.imstat(
        idl_proc / bin_spec_sn30 / cont_flux_ *.fits, fields = "mean", format = no > idl_proc / bin_spec_sn30 / binned_cont_flux.lis
    '''


if __name__ == '__main__':
    """
    Follow the steps necessary to prepare input to pPXF, then run pPXF
    """

    if not os.path.exists(XYSN_FILE):
        print('>>>>> Making table')
        make_table()
    if not os.path.exists(V2B_FILE):
        print('>>>>> Voronoi binning')
        t = clock()
        voronoi_binning()
        print('Elapsed time: %.2f seconds' % (clock() - t))
    if not os.path.exists(os.path.join(DIR_SCI, FILE_VAR.format(1))):
        print('>>>>> Binning spectra')
        bin_spectra()
    if not os.path.exists(os.path.join(DIR_SCI, BIN_VAR.format(1))):
        print('>>>>> Combining spectra')
        combine_spectra()
    '''
    if not os.path.exists(OUT_VEL_SIGMA):
        print('>>>>> pPXF')
        ppxf_kinematics()
    '''
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
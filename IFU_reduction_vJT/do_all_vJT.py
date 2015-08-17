# written by Kristi Webb for data of IC225 taken by Bryan Miller
# similar to an earlier method developed by Gwen Rudie in IDL/cl/awk

# from __future__ import print_function
from time import clock
import glob
import os
import numpy as np
from voronoi_2d_binning import voronoi_2d_binning
from scipy import ndimage, signal
from ppxf import ppxf
import ppxf_util as util
from matplotlib import pyplot as plt
from astropy.io import fits
import pandas as pd
from cap_plot_velfield import plot_velfield
import remove_lines
import bin_by_object
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import specutils as su
#import ppxf_IC225_vJT


"""
This script covers all the kinematic analysis of a reduced 3D image cube by the method in reduce_and_cubel.cl for the
  example case of IC225

To visualize the 2D flattened cube for a given spectral range the steps are as follows:
    - Crop the 3D cube for a given spectral range
            scrop_cube(IMAGE_CUBE, SCROP_RANGE, CUBE_SCROPD)
    - Flatten the cropped 3D cubed
            imcombine_flatten(CUBE_SCROPD, SCI_EXT_SCROPD, VAR_EXT_SCROPD)
*** These flattened images are not used in the kinematic analysis, but are useful for visual confirmation of
    relative flux levels between the nucleui and surrounding gas

    The relavent nebular emission lines are:
        - OIII: 4940-5065 A - rest 4958.92, 5006.84
        - H_beta: 4865-4915 A - rest 4861.33
        - H_gamma: 4330-4380 A - rest 4340.47
    The continuum spectrum is:
        - continuum: 4383-4850 A

To run pPXF:
    - Compile coordinates of signal and noise measurements for the binning
            make_table(IMAGE_CUBE, SCI_EXT, VAR_EXT, XYSN_FILE)
    - Preform voronoi 2d binning to acheive desired signal to noise
            voronoi_binning(XYSN_FILE, V2B_FILE, V2B_XY_FILE)
    - Combine spectra of the same bin into fits files
            combine_spectra(V2B_FILE, IMAGE_CUBE, BIN_SCI, FLUX_SCI, BIN_VAR, FLUX_VAR)
    - Determine the optimal BIAS, first run ppxf with BIAS=0 to find the h3 and h4 parameters, put these into
      ppxf_simulation and determine the optimal BIAS. See the information in ppxf.py or in the ppxf readme.txt for
      more information about this step. I have selected BIAS=0.6 as a good value.
            ppxf_simulation(PPXF_BESTFIT.strip('.fits') + '_bias0.fits', LAM_RANGE, TARGET_SN, bias=0.5, spaxel=0)
    - Preform pPXF on the binned spectra
           ppxf_kinematics(BIN_SCI, PPXF_FILE, PPXF_BESTFIT, TEMPLATE_FITS, TEMPLATE_RESOLUTION, LAM_RANGE, VEL_INIT, SIG_INIT, bias=0.6)
To plot pPXF output:
    - Copy the desired spectral range and read mean flux values to be used for plotting
            scopy_flux(FLUX_SCI, FLUX_SCOPY_FITS, FLUX_SCOPY_RANGE, FLUX_SCOPY_FILE)
    - Access arrays to input into cap_plot_velfield(), set vel as velocity information obtained from ppxf
            plot_velfield_setup(vel, V2B_XY_FILE, FLUX_SCOPY_FILE)

As the spectra contain both emission and absorption features, it is important to separate the two before using fxcor
or rvsao as (unlike ppxf) they do not remove the emission lines
    - Create a spectra of the emission lines to study the gas kinematics
            remove_lines.remove_absorp_lines(BIN_SCI, PPXF_BESTFIT, EM_BIN_SCI)
    - Create a spectra of the emission lines to study the stellar kinematics
            remove_lines.remove_emission_lines(BIN_SCI, ABS_BIN_SCI, PPXF_BESTFIT, plot=False)
    *** xcsao is 'supposed' to have the ability to remove emission lines with 'em_chop+' but this kept returning very
        strange results, so I instead use the fits files with the emission lines removed manually.

To run fxcor to preform cross correlation to determine relative velocity:
    - Create input file from list of spectral bins and select the bin to use as the template
    - Absorption lines
            fxcor(EM_BIN_SCI, TEMPLATE_SPECTRA, EM_FXCOR_BIN_LIST, EM_FXCOR_FILE)
    - Emission lines
            fxcor(ABS_BIN_SCI, TEMPLATE_SPECTRA, ABS_FXCOR_BIN_LIST, ABS_FXCOR_FILE)
To plot the fxcor output:
    - Access arrays to input into cap_plot_velfield(), set vel as velocity information obtained from fxcor
            plot_velfield_setup(vel, V2B_XY_FILE, FLUX_SCOPY_FILE)

To run rvsao (with xcsao or emsao) to preform cross correlation to determine relative velocity:
    - Use the absorption lines spectra to study the relvative velocity of the stars
            rvsao(ABS_BIN_SCI, 'xcsao', TEMPLATE_SPECTRA, XCSAO_FILE, XCSAO_BIN_LIST)
    - Use the emission lines spectra to study the relvative velocity of the gas
            rvsao(EM_BIN_SCI, 'emsao', TEMPLATE_SPECTRA, EMSAO_FILE, EMSAO_BIN_LIST)
To plot, same method as above, with vel as the velcity information obtained from rvsao
            plot_velfield_setup(vel, V2B_XY_FILE, FLUX_SCOPY_FILE)


Template spectra for fxcor/rvsoa (bin with high SN)
    SN 15 - spectra 41 has SN=22.3
    SN 30 - spectra 3 has SN=45.8
    SN 20 - spectra 22 has SN=30.4537 (using spectra 3 here gives a weird values at bin 30)

"""

# These are the files to change based on desired use
# --------------------------------------------------

DIR_PATH = '/Users/kwebb/IFU_reduction_vJT'  # Working directory (where the 3D cube is)
IMAGE_CUBE = os.path.join(DIR_PATH, 'dchsteqpxbprgN20051205S0006_scropd_add.fits')
TARGET_SN = 40  # REMEMBER TO CHANGE TEMPLATE SPECTRA FOR EACH S/N
TEMPLATE_SPECTRA = 8  # bin 1,6 at SN 20 is about in the middle of the frame

# To create 2D flattened science and variance images of specific spectral range
# This is only necessary for visual comparison of the intensity of regions for a given spectral range
SCROP_RANGE = [4330, 4380]  # wavelength range to scrop the image cube to which will be flattened
SCROP_PATH = os.path.join(DIR_PATH, 'scrop_proc')  # contains the output of all the scrop-type methods
CUBE_SCROPD = os.path.join(SCROP_PATH, 'IC225_3D_scropd_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))
SCI_EXT_SCROPD = os.path.join(SCROP_PATH, 'IC225_2D_sci_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))
VAR_EXT_SCROPD = os.path.join(SCROP_PATH, 'IC225_2D_var_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))

# These variables define the chosen file structure I use to organise the output
# -----------------------------------------------------------------------------

# Flattened 3D cubes images
SCI_EXT = os.path.join(DIR_PATH, 'IC225_2D_sci.fits')
VAR_EXT = os.path.join(DIR_PATH, 'IC225_2D_var.fits')

PROC_PATH = os.path.join(DIR_PATH, 'doall_proc_{}/voronoi_binning/'.format(TARGET_SN))  # , 'region_ocn')
# PROC_PATH = os.path.join(DIR_PATH, 'doall_proc_{}/spec_binning/'.format(TARGET_SN))  # , 'region_ocn')
# PROC_PATH = os.path.join(DIR_PATH, 'doall_proc_{}/'.format(TARGET_SN))  # , 'region_ocn')

# Organise output of Voronoi binning
XYSN_FILE = os.path.join(PROC_PATH, 'y_x_signal_noise.txt')  # output of make_table
V2B_FILE = os.path.join(PROC_PATH, 'v2b_output_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning
V2B_XY_FILE = os.path.join(PROC_PATH, 'v2b_output_xy_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning

# Files with equivalent output as v2b method but with self defined bins
SRCBIN_FILE = os.path.join(PROC_PATH, 'srcbin_output_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning
SRCBIN_XY_FILE = os.path.join(PROC_PATH, 'srcbin_output_xy_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning

# Organise combined spectra folder and names of binned (summed) and flux'd (averaged) spectra
DIR_SCI_COMB = os.path.join(PROC_PATH, 'comb_fits_sci_{}'.format(TARGET_SN))  # folder with combined spectra
DIR_VAR_COMB = os.path.join(PROC_PATH, 'comb_fits_var_{}'.format(TARGET_SN))  # folder with combined var spec
BIN_SCI = os.path.join(DIR_SCI_COMB, 'bin_sci_{}.fits')  # naming convention of combined (sum) sci spectra
BIN_VAR = os.path.join(DIR_VAR_COMB, 'bin_var_{}.fits')  # naming convention of combined (sum) var spectra
FLUX_SCI = os.path.join(DIR_SCI_COMB, 'flux_sci_{}.fits')  # naming convention of combined (average) sci spectra
FLUX_VAR = os.path.join(DIR_VAR_COMB, 'flux_var_{}.fits')  # naming convention of combined (average) var spectra

EM_BIN_SCI = os.path.join(DIR_SCI_COMB, 'em_bin_sci_{}.fits')
ABS_BIN_SCI = os.path.join(DIR_SCI_COMB, 'abs_bin_sci_{}.fits')

# Organise output of pPXF
PPXF_PATH = os.path.join(PROC_PATH, 'ppxf_proc')
PPXF_FILE = os.path.join(PPXF_PATH, 'ppxf_output_sn{}.txt'.format(TARGET_SN))
PPXF_BESTFIT = os.path.join(PPXF_PATH, 'bestfit_{}.fits')

# Organise output of fxcor
FXCOR_PATH = os.path.join(PROC_PATH, 'fxcor_proc')
EM_FXCOR_FILE = os.path.join(FXCOR_PATH, 'fxcor_em_bin_sci_sn{}_tmpl{}'.format(TARGET_SN, TEMPLATE_SPECTRA))
EM_FXCOR_BIN_LIST = os.path.join(FXCOR_PATH, 'em_bin_sci_sn{}_list.lis'.format(TARGET_SN))
ABS_FXCOR_FILE = os.path.join(FXCOR_PATH, 'fxcor_abs_bin_sci_sn{}_tmpl{}'.format(TARGET_SN, TEMPLATE_SPECTRA))
ABS_FXCOR_BIN_LIST = os.path.join(FXCOR_PATH, 'abs_bin_sci_sn{}_list.lis'.format(TARGET_SN))

# Organize output of rvsao
RVSAO_PATH = os.path.join(PROC_PATH, 'rvsao_proc')

BIN_LIST = os.path.join(RVSAO_PATH, 'bin_sci_sn{}_list.lis'.format(TARGET_SN))

XCSAO_BIN_LIST = os.path.join(RVSAO_PATH, 'abs_bin_sci_sn{}_list.lis'.format(TARGET_SN))
XCSAO_TEMPLATE = os.path.join(RVSAO_PATH, BIN_SCI.format(TEMPLATE_SPECTRA))
XCSAO_FILE = os.path.join(RVSAO_PATH, 'xcsao_bin_sci_sn{}_tmpl{}.txt'.format(TARGET_SN, TEMPLATE_SPECTRA))
EMSAO_FILE = os.path.join(RVSAO_PATH, 'emsao_bin_sci_sn{}_tmpl{}.txt'.format(TARGET_SN, TEMPLATE_SPECTRA))
EMSAO_BIN_LIST = os.path.join(RVSAO_PATH, 'em_bin_sci_sn{}_list.lis'.format(TARGET_SN))

# pPXF parameters
VEL_INIT = 1500.  # initial guess for velocity
SIG_INIT = 150.  # inital guess of sigma distribution
LAM_RANGE = [4173.89, 5404.51]  # wavelength range for logarithmic rebinning (full range)
# template from MILES Library spanning 3540-7410 A, with resolution spectral resolution of 2.54 A (FWHM),
# sigma~64 km/s, R~2000. From Sanchez-Blazquez, et al. (2006)
# (http://www.iac.es/proyecto/miles/pages/stellar-libraries/miles-library.php)
# TEMPLATE_FITS = '/Users/kwebb/idl/cappellari/ppxf/spectra/Mun1.30z*.fits'  # Used by Gwen
TEMPLATE_FITS = '/Users/kwebb/python/ppxf/MILES_BaSTI_un_1.30_fits/Mun1.30Z*.fits'
TEMPLATE_RES = 2.54  # FWHM of the template spectra *** IN ANGSTROMS ***

# To create combined (averaged) spectra to determine mean flux of a specific wavelength range to plot
# YOU DO NOT NEED TO CROP THE SPECTRA HERE, it is just a good idea to check that you measure the same velocity for
# any given spectral range
FLUX_SCOPY_RANGE = [4173.89, 5404.51]
FLUX_SCOPY_FITS_SUFFIX = 'flux_scopy_{}.fits'  # wavelength range combined sci spectra from flux*
FLUX_SCOPY_FILE = os.path.join(PROC_PATH, 'binned_flux_{}.txt'.format(TARGET_SN))
FLUX_SCOPY_FITS = os.path.join(DIR_SCI_COMB, FLUX_SCOPY_FITS_SUFFIX)

XYSN_FILE_CN = os.path.join(PROC_PATH, 'y_x_signal_noise_cn.txt')  # output of make_table
XYSN_FILE_OCN = os.path.join(PROC_PATH, 'y_x_signal_noise_ocn.txt')  # output of make_table


def flatten_cube(image_cube, sci_ext, var_ext):
    """
    equivalent of last step of preprocessing
    """

    from pyraf import iraf

    iraf.images()

    if not os.path.exists(sci_ext):
        iraf.imcombine(image_cube + '[sci]', sci_ext, project="yes")
    if not os.path.exists(var_ext):
        iraf.imcombine(image_cube + '[var]', var_ext, project="yes")


def scrop_cube(image_cube, scrop_range, cube_scropd):
    """
    As the output cube of the preprocessing is NOT an MEF but instread a simple fits file with a SCI and VAR extension
    Gwen's scrop task cannot be used. This should replace the functionality.
    This takes a 3D cube and returns a 3D cube cropped to the wavelength specified in scrop_range

    NOTE: Still need to implement changes to header values

    INPUT: IMAGE_CUBE (dcsteqpxbprgN20051205S0006_add.fits), SCROP_RANGE ([4360, 4362])
    OUTPUT: CUBE_SCOPYD
    """

    if os.path.exists(cube_scropd):
        print('File {} already exists'.format(cube_scropd))
        return

    print ('Scropping cube {} to wavelength range {}'.format(image_cube, scrop_range))

    with fits.open(image_cube) as cube_hdu:
        # cube_hdu.info()
        # Filename: dcsteqpxbprgN20051205S0006_add.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU     214   ()
        # 1    SCI         ImageHDU        68   (76, 49, 1300)   float32
        # 2    VAR         ImageHDU        68   (76, 49, 1300)   float32
        cube_data0 = cube_hdu[0].data
        cube_data1 = cube_hdu[1].data
        cube_data2 = cube_hdu[2].data
        cube_hdr0 = cube_hdu[0].header
        cube_hdr1 = cube_hdu[1].header
        cube_hdr2 = cube_hdu[2].header

    crval3 = cube_hdr1['CRVAL3']
    crpix3 = cube_hdr1['CRPIX3']
    cd33 = cube_hdr1['CD3_3']
    npix = cube_hdr1['NAXIS3']

    wmin = crval3 + (1. - crpix3) * cd33
    wmax = crval3 + (npix - crpix3) * cd33
    # dwav = cd33

    assert scrop_range[0] >= wmin, 'Cannot crop spectra outside of wavelength range [{},{}]'.format(wmin, wmax)
    assert scrop_range[1] <= wmax, 'Cannot crop spectra outside of wavelength range [{},{}]'.format(wmin, wmax)

    x1 = (scrop_range[0] - crval3) / cd33 + crpix3
    x2 = (scrop_range[1] - crval3) / cd33 + crpix3

    scrop_cube1 = cube_data1[x1:x2, :, :]
    scrop_cube2 = cube_data2[x1:x2, :, :]

    cube_hdr1['NAXIS3'] = scrop_cube1.shape[0]
    cube_hdr2['NAXIS3'] = scrop_cube2.shape[0]

    hdu_out = fits.HDUList()
    hdu_out.append(fits.PrimaryHDU(data=cube_data0, header=cube_hdr0))
    hdu_out.append(fits.ImageHDU(data=scrop_cube1, name='SCI', header=cube_hdr1))
    hdu_out.append(fits.ImageHDU(data=scrop_cube2, name='VAR', header=cube_hdr2))
    hdu_out.header = cube_hdr0
    hdu_out.writeto(cube_scropd)


def imcombine_flatten(cube, sci_ext, var_ext):
    """
    Imcombine both variance and science extensions to flatten the cube to 2D
    INPUT: CUBE (scrop_proc/IC225_3D_{}_{}.fits)
    OUTPUT: SCI_EXT (scrop_proc/IC225_2D_sci_{}_{}.fits), VAR_EXT (scrop_proc/IC225_2D_var_{}_{}.fits)
    """

    if os.path.exists(sci_ext):
        print('File {} already exists'.format(sci_ext))
        return
    if os.path.exists(var_ext):
        print('File {} already exists'.format(var_ext))
        return

    from pyraf import iraf

    iraf.imcombine('{}[sci]'.format(cube), sci_ext, project="yes")
    iraf.imcombine('{}[var]'.format(cube), var_ext, project="yes")


def make_table(image_cube, sci_ext, var_ext, xysn_file):
    """
    Read in the pixel values of the science and varience planes of the image cube and make a table of the
    coordinates with the respective signal and noise measurements
    INPUT: SCI_EXT (IC225_2D_sci.fits), VAR_EXT (IC225_2D_var.fits)
    OUTPUT: XYSN_FILE (x_y_signal_noise.txt)
    """

    if os.path.exists(xysn_file):
        print('File {} already exists'.format(xysn_file))
        return

    if not os.path.exists(sci_ext):
        imcombine_flatten(image_cube, sci_ext, var_ext)

    assert os.path.exists(sci_ext), 'Image {} does not exist'.format(sci_ext)
    assert os.path.exists(var_ext), 'Image {} does not exist'.format(var_ext)

    with fits.open(sci_ext) as sci_hdu:
        sci_data = sci_hdu[0].data
        sci_xaxis = sci_hdu[0].header['NAXIS1']
        sci_yaxis = sci_hdu[0].header['NAXIS2']

    with fits.open(var_ext) as var_hdu:
        var_data = var_hdu[0].data
        var_xaxis = var_hdu[0].header['NAXIS1']
        var_yaxis = var_hdu[0].header['NAXIS2']

    assert sci_xaxis == var_xaxis, 'Sci and var planes have diffent x dimensions'
    assert sci_yaxis == var_yaxis, 'Sci and var planes have diffent y dimensions'

    with open(xysn_file, 'w') as outfile:
        for i in range(sci_yaxis - 1):
            for j in range(sci_xaxis - 1):
                noise = np.sqrt(var_data[i, j])
                outfile.write('   {}   {}   {}   {}\n'.format(i, j, sci_data[i, j], noise))


def voronoi_binning(xysn_file, v2b_file, v2b_xy_file):
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

    if os.path.exists(v2b_file):
        print('File {} already exists'.format(v2b_file))
        return

    y, x, signal, noise = np.loadtxt(xysn_file, unpack=True)  # , skiprows=3)

    # Only select pixels where the signal and noise is nonzero

    ''' CHECK VALIDITY OF THIS '''
    noise = noise[noise != 0]
    signal = signal[noise != 0]
    x = x[noise != 0]
    y = y[noise != 0]

    # Perform the actual computation. The vectors (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
    # are all generated in *output*

    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(x, y, signal, noise, TARGET_SN,
                                                                              plot=1, quiet=0)

    # Save to a text file the initial coordinates of each pixel together with the corresponding bin number computed
    # by this procedure. binNum uniquely specifies the bins and for this reason it is the only
    # number required for any subsequent calculation on the bins.

    np.savetxt(v2b_file, np.column_stack([x, y, binNum]), header='x  y  binNum', fmt=b'%10.6f %10.6f %8i')
    np.savetxt(v2b_xy_file, np.column_stack([xBar, yBar, xNode, yNode]), header='xBar  yBar  xNode   yNode',
               fmt=b'%10.6f %10.6f %10.6f %10.6f')


def combine_spectra(v2b_file, image_cube, bin_sci, flux_sci, bin_var, flux_var):
    """
    Combine each pixel of the same bin (according to the lists output by bin_spectra) into a single fits file
    INPUT: V2B_FILE (v2b_output_sn30.txt)
    OUTPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), DIR_VAR_COMB (comb_fis_var_{S/N}/bin_var_{S/N}.fits)
    """

    # read in the output of the voronoi binning
    v2b_output = pd.read_table(v2b_file, sep=r"\s*", engine='python', skiprows=1, names=["x", "y", "bin"])

    # for each bin, combine the spectra of the same bin
    for i in range(np.amax(v2b_output['bin'].values) + 1):
        bin_sci_i = bin_sci.format(i)
        bin_var_i = bin_var.format(i)
        flux_sci_i = flux_sci.format(i)
        flux_var_i = flux_var.format(i)

        if not os.path.exists(bin_sci_i):
            # Check to see how many pixels are in this bin
            spaxel_list = v2b_output.query('bin == {}'.format(i))
            print('Number of pixels in bin {}: {}'.format(i, len(spaxel_list)))

            imcombine(image_cube, bin_sci_i, bin_var_i, spaxel_list, 'sum')
            imcombine(image_cube, flux_sci_i, flux_var_i, spaxel_list, 'average')


def imcombine(image_cube, outsci, outvar, spaxel_list, combine_method):
    """
    In order to avoid a 'floating point error' in imexamine I'm going to try and combine the fits sections as
    numpy arrays instead and manually copy the header values. Only works for combine='average'|'sum' at the moment
    """

    assert combine_method is 'sum' or 'average', 'Combine is not "sum" or "average"'

    # IMAGE_CUBE[int(line["x"].values[j - 1]) + 1, int(line["y"].values[j - 1]) + 1, *]

    with fits.open(image_cube) as image_cube_hdu:
        # Filename: dcsteqpxbprgN20051205S0006_add.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU     214   ()
        # 1    SCI         ImageHDU        68   (76, 49, 1300)   float32
        # 2    VAR         ImageHDU        68   (76, 49, 1300)   float32
        image_cube_sci = image_cube_hdu[1].data
        image_cube_var = image_cube_hdu[2].data
        cube_header = image_cube_hdu['SCI'].header
        cube_header_0 = image_cube_hdu[0].header

    imcomb_sci = []
    imcomb_var = []

    # image_cube_sci.shape = (1300, 49, 76) ie (z,y,x)
    for k in range(image_cube_sci.shape[0]):
        pix_sci = []
        pix_var = []
        for j in range(len(spaxel_list)):
            pix_sci.append(image_cube_sci[k, spaxel_list["y"].values[j - 1] + 1, spaxel_list["x"].values[j - 1] + 1])
            pix_var.append(image_cube_var[k, spaxel_list["y"].values[j - 1] + 1, spaxel_list["x"].values[j - 1] + 1])

        if combine_method == 'sum':
            imcomb_sci.append(np.array(pix_sci).sum())
            imcomb_var.append(np.array(pix_var).sum())
        elif combine_method == 'average':
            imcomb_sci.append(np.array(pix_sci).mean())
            imcomb_var.append(np.array(pix_var).mean())

    write_imcomb_fits(imcomb_sci, outsci, cube_header, cube_header_0)
    write_imcomb_fits(imcomb_var, outvar, cube_header, cube_header_0)


def write_imcomb_fits(outdata, outfile, cube_header, cube_header_0):
    """
    Write the combined spectra to a fits file with headers from the wavelength dimension of the 3D cube
    """
    outfile_hdu = fits.PrimaryHDU()
    outfile_hdu.header['OBSERVAT'] = cube_header_0['OBSERVAT']
    outfile_hdu.header['CD1_1'] = cube_header['CD3_3']
    outfile_hdu.header['CRPIX1'] = cube_header['CRPIX3']
    outfile_hdu.header['CRVAL1'] = cube_header['CRVAL3']

    outfile_hdu.header['DATE'] = cube_header_0['DATE']
    outfile_hdu.header['INSTRUME'] = cube_header_0['INSTRUME']
    outfile_hdu.header['OBJECT'] = cube_header_0['OBJECT']
    outfile_hdu.header['GEMPRGID'] = cube_header_0['GEMPRGID']
    outfile_hdu.header['OBSID'] = cube_header_0['OBSID']
    outfile_hdu.header['OBSERVAT'] = cube_header_0['OBSERVAT']
    outfile_hdu.header['TELESCOP'] = cube_header_0['TELESCOP']
    outfile_hdu.header['EQUINOX'] = cube_header_0['EQUINOX']
    outfile_hdu.header['EPOCH'] = cube_header_0['EPOCH']
    outfile_hdu.header['RA'] = cube_header_0['RA']
    outfile_hdu.header['DEC'] = cube_header_0['DEC']
    outfile_hdu.header['DATE-OBS'] = cube_header_0['DATE-OBS']
    outfile_hdu.header['TIME-OBS'] = cube_header_0['TIME-OBS']
    outfile_hdu.header['UTSTART'] = cube_header_0['UTSTART']
    outfile_hdu.header['UTEND'] = cube_header_0['UTEND']
    outfile_hdu.header['EXPTIME'] = cube_header_0['EXPTIME']

    outfile_hdu.header['DC-FLAG'] = 0  # ie linear binning (1 = linear-log binning)
    # Header values readable by rvsao
    outfile_hdu.header['DEC--TAN'] = 'LAMBDA'
    outfile_hdu.data = outdata
    outfile_hdu.writeto(outfile, clobber=True)


def ppxf_kinematics(bin_sci, ppxf_file, ppxf_bestfit, template_fits, template_resolution, lam_range=[4173., 5404.],
                    vel_init=1500., sig_init=100., bias=0., plot=False, quiet=True, badPixels=[], clean=False):
    """
    Follow the pPXF usage example by Michile Cappellari
    INPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), TEMPLATE_* (spectra/Mun1.30z*.fits)
    OUTPUT: PPXF_FILE (ppxf_output_sn30.txt), OUT_BESTFIT_FITS (ppxf_fits/ppxf_bestfit_sn30.fits)
    """

    if os.path.exists(ppxf_file):
        print('File {} already exists'.format(ppxf_file))
        return

    # Read a galaxy spectrum and define the wavelength range
    assert len(glob.glob(bin_sci.format('*'))) > 0, 'Binned spectra {} not found'.format(glob.glob(bin_sci.format('*')))

    # hdu = fits.open(in_file[0])
    # gal_lin = hdu[0].data
    # h1 = hdu[0].header

    with fits.open(bin_sci.format(0)) as gal_hdu:
        gal_lin = gal_hdu[0].data
        gal_hdr = gal_hdu[0].header

    # lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
    # FWHM_gal = 4.2 # SAURON has an instrumental resolution FWHM of 4.2A.

    # lamRange1 is now variable lam_range (because I was lazy and didnt put headers into the binned spectra)
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

    galaxy, logLam1, velscale = util.log_rebin(lam_range, gal_lin)
    galaxy = galaxy / np.median(galaxy)  # Normalize spectrum to avoid numerical issues

    # Read the list of filenames from the Single Stellar Population library by Vazdekis (1999, ApJ, 513, 224). A subset
    # of the library is included for this example with permission. See http://purl.org/cappellari/software
    # for suggestions of more up-to-date stellar libraries.

    # vazdekis = glob.glob(dir + 'Rbi1.30z*.fits')
    # vazdekis.sort()
    # FWHM_tem = 1.8 # Vazdekis spectra have a resolution FWHM of 1.8A.

    template_spectra = glob.glob(template_fits)
    assert len(template_spectra) > 0, 'Template spectra not found: {}'.format(template_fits)
    FWHM_tem = template_resolution

    # Extract the wavelength range and logarithmically rebin one spectrum to the same velocity scale of the SAURON
    # galaxy spectrum, to determine the size needed for the array which will contain the template spectra.

    # hdu = fits.open(vazdekis[0])
    # ssp = hdu[0].data
    # h2 = hdu[0].header
    # lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    # sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    # templates = np.empty((sspNew.size,len(vazdekis)))

    with fits.open(template_spectra[0]) as temp_hdu:
        ssp = temp_hdu[0].data
        h2 = temp_hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0., h2['CDELT1'] * (h2['NAXIS1'] - 1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    templates = np.empty((sspNew.size, len(template_spectra)))

    # Convolve the whole Vazdekis library of spectral templates with the quadratic difference between the SAURON and the
    # Vazdekis instrumental resolution. Logarithmically rebin and store each template as a column in the array TEMPLATES

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

    # FWHM_dif = np.sqrt(FWHM_gal ** 2 - template_resolution ** 2)
    FWHM_dif = np.sqrt(FWHM_tem ** 2 - FWHM_gal ** 2)
    sigma = FWHM_dif / 2.355 / h2['CDELT1']  # SIGMA DIFFERENCE IN PIXELS, 1.078435697220085

    # Logarithmically rebin the whole Mun library of spectra, and store each template as a column in the array TEMPLATES

    # for j in range(len(vazdekis)):
    # hdu = fits.open(vazdekis[j])
    # ssp = hdu[0].data
    # ssp = ndimage.gaussian_filter1d(ssp, sigma)
    # sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    # templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

    for j in range(len(template_spectra)):
        with fits.open(template_spectra[j]) as temp_hdu_j:
            ssp_j = temp_hdu_j[0].data
        ssp_j = ndimage.gaussian_filter1d(ssp_j, sigma)
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp_j, velscale=velscale)

        templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

    # The galaxy and the template spectra do not have the same starting wavelength. For this reason an extra velocity
    # shift DV has to be applied to the template to fit the galaxy spectrum. We remove this artificial shift by using
    # the keyword VSYST in the call to PPXF below, so that all velocities are measured with respect to DV. This assume
    # the redshift is negligible.In the case of a high-redshift galaxy one should de-redshift its wavelength to the
    # rest frame before using the line below (see above).

    vel_list = []
    sig_list = []
    dV_list = []
    dsigma_list = []

    h3_list = []
    h4_list = []
    dh3_list = []
    dh4_list = []

    for j in range(len(glob.glob(bin_sci.format('*')))):
        b_gal = fits.getdata(bin_sci.format(j), 0)

        b_gal = ndimage.gaussian_filter1d(b_gal, sigma)

        galaxy, logLam1, velscale = util.log_rebin(lam_range, b_gal, velscale=velscale)
        noise = galaxy * 0 + 1  # Assume constant noise per pixel here

        c = 299792.458
        dv = (logLam2[0] - logLam1[0]) * c  # km/s

        # vel = 1500.  # Initial estimate of the galaxy velocity in km/s
        z = np.exp(vel_init / c) - 1  # Relation between velocity and redshift in pPXF

        goodPixels = util.determine_goodpixels(logLam1, lamRange2, z)
        if len(badPixels) > 0:
            indices = []
            for idx, pix in enumerate(goodPixels):
                if pix in badPixels:
                    indices.append(idx)
            goodPixels = np.delete(goodPixels, indices)

        # Here the actual fit starts. The best fit is plotted on the screen.
        # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.

        start = [vel_init, sig_init]  # (km/s), starting guess for [V,sigma]
        t = clock()

        pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodPixels, plot=plot, moments=4,
                  degree=4, vsyst=dv, bias=bias, quiet=False, clean=clean)

        if not quiet:
            print("Formal errors:")
            print("     dV    dsigma   dh3      dh4")
            print("".join("%8.2g" % f for f in pp.error * np.sqrt(pp.chi2)))

        # If the galaxy is at significant redshift z and the wavelength has been de-redshifted with the three lines
        # "z = 1.23..." near the beginning of this procedure, the best-fitting redshift is now given by the following
        # commented line (equation 2 of Cappellari et al. 2009, ApJ, 704, L34;
        # http://adsabs.harvard.edu/abs/2009ApJ...704L..34C)
        # print, 'Best-fitting redshift z:', (z + 1)*(1 + sol[0]/c) - 1

        # Gwen obtains the velocity and sigma information from the SOL parameter
        # moments = 4 so sol = [vel, sig, h3, h4]
        vel_list.append(pp.sol[0])
        sig_list.append(pp.sol[1])
        dV_list.append((pp.error * np.sqrt(pp.chi2))[0])
        dsigma_list.append((pp.error * np.sqrt(pp.chi2))[1])

        h3_list.append(pp.sol[2])
        h4_list.append(pp.sol[3])
        dh3_list.append((pp.error * np.sqrt(pp.chi2))[2])
        dh4_list.append((pp.error * np.sqrt(pp.chi2))[3])

        hdu_best = fits.PrimaryHDU()
        hdu_best.data = pp.bestfit
        hdu_best.header['CD1_1'] = gal_hdr['CD1_1']
        hdu_best.header['CDELT1'] = gal_hdr['CD1_1']
        hdu_best.header['CRPIX1'] = gal_hdr['CRPIX1']
        hdu_best.header['CRVAL1'] = gal_hdr['CRVAL1']
        hdu_best.header['NAXIS1'] = pp.bestfit.size
        hdu_best.header['CTYPE1'] = 'LINEAR'  # corresponds to sampling of values above
        hdu_best.header['DC-FLAG'] = '1'  # 0 = linear, 1= log-linear sampling
        hdu_best.writeto(ppxf_bestfit.format(j), clobber=True)

    print('Elapsed time in PPXF: %.2f s' % (clock() - t))

    np.savetxt(ppxf_file,
               np.column_stack([vel_list, sig_list, h3_list, h4_list, dV_list, dsigma_list, dh3_list, dh4_list]),
               fmt=b'%10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f',
               header='velocity    sigma       h3           h4         dV          dsigma       dh3         dh4')

    return vel_list


def ppxf_simulation(ppxf_bestfit, lam_range, target_sn, bias=0.6, spaxel=0):
    """
    2. Perform a fit of your kinematics *without* penalty (PPXF keyword BIAS=0).
       The solution will be noisy and may be affected by spurious solutions,
       however this step will allow you to check the expected mean ranges in
       the Gauss-Hermite parameters [h3,h4] for the galaxy under study;

        see ppxf_output_sn30_bias0.txt
        mean(h3), std(h3), max(h3) = -0.2555, 0.09090, 0.007468
         # ignoring endpoints, max(h3[1:-1]) = -0.01376
        mean(h4), std(h4), max(h4) = -0.07712, 0.1423, 0.136607
         # ignoring endpoints, max(h4[1,-1]) = 0.13594
        max(dvel), min(dvel), max(dvel)-np.min(dvel) = 119.2918, 4.08643, 115.20544
        mean(vel) = 1543.0359
        max(sig), min(sig) = 180.3, 36.23



    3. Perform a Monte Carlo simulation of your spectra, following e.g. the
       included ppxf_simulation_example.pro routine. Adopt as S/N in the simulation
       the chosen value (S/N)_min and as input [h3,h4] the maximum representative
       values measured in the non-penalized pPXF fit of the previous step;

    4. Choose as penalty (BIAS) the *largest* value such that, for sigma > 3*velScale,
       the mean difference between the output [h3,h4] and the input [h3,h4]
       is well within the rms scatter of the simulated values
       (see e.g. Fig.2 of Emsellem et al. 2004, MNRAS, 352, 721).
    """

    # dir = 'spectra/'
    # file = dir + 'Rbi1.30z+0.00t12.59.fits'
    # hdu = pyfits.open(file)
    # ssp = hdu[0].data
    # h = hdu[0].header

    bestfit_file = ppxf_bestfit.format(spaxel)
    assert os.path.exists(bestfit_file), 'Best fit spectra not found: {}'.format(bestfit_file)

    with fits.open(bestfit_file) as best_hdu:
        ssp = best_hdu[0].data
        h = best_hdu[0].header

    # lamRange = h['CRVAL1'] + np.array([0., h['CDELT1'] * (h['NAXIS1'] - 1)])
    # star, logLam, velscale = util.log_rebin(lamRange, ssp)

    star, logLam, velscale = util.log_rebin(lam_range, ssp)

    # The finite sampling of the observed spectrum is modeled in detail: the galaxy spectrum is obtained by oversampling
    # the actual observed spectrum to a high resolution. This represents the true spectrum, which is later resampled
    # to lower resolution to simulate the observations on the CCD. Similarly, the convolution with a well-sampled LOSVD
    # is done on the high-resolution spectrum, and later resampled to the observed resolution before fitting with PPXF.

    factor = 10  # Oversampling integer factor for an accurate convolution
    starNew = ndimage.interpolation.zoom(star, factor, order=1)  # The underlying spectrum, known at high resolution
    star = rebin(starNew, factor)  # Make sure that the observed spectrum is the integral over the pixels

    # vel = 0.3  # velocity in *pixels* [=V(km/s)/velScale]
    # h3 = 0.1  # Adopted G-H parameters of the LOSVD
    # h4 = -0.1
    # sn = 60.  # Adopted S/N of the Monte Carlo simulation
    # m = 300  # Number of realizations of the simulation
    # sigmaV = np.linspace(0.8, 4, m)  # Range of sigma in *pixels* [=sigma(km/s)/velScale]

    print('bestfit {} : velscale = {}'.format(spaxel, velscale))
    # bestfit 0 : velscale = 57.44245804
    # max(sig), min(sig) = 164.026343, 11.306016 [km/s] =  2.87765514,  0.19835116 [pix]

    vel = 0.3  # mean vel = 1556.4989, I have no idea why this value is 0.3 ...
    h3 = -0.003146
    h4 = -0.003662
    sn = target_sn
    m = 300  # Number of realizations of the simulation
    sigmaV = np.linspace(0.8, 4, m)  # Range of sigma in *pixels* [=sigma(km/s)/velScale]
    # sigmaV = np.linspace(0.198, 2.877, m)

    result = np.zeros((m, 4))  # This will store the results
    t = clock()
    np.random.seed(123)  # for reproducible results

    for j in range(m):
        sigma = sigmaV[j]
        dx = int(abs(vel) + 4.0 * sigma)  # Sample the Gaussian and GH at least to vel+4*sigma
        x = np.linspace(-dx, dx, 2 * dx * factor + 1)  # Evaluate the Gaussian using steps of 1/factor pixels.
        w = (x - vel) / sigma
        w2 = w ** 2
        gauss = np.exp(-0.5 * w2) / (np.sqrt(2. * np.pi) * sigma * factor)  # Normalized total(gauss)=1
        h3poly = w * (2. * w2 - 3.) / np.sqrt(3.)  # H3(y)
        h4poly = (w2 * (4. * w2 - 12.) + 3.) / np.sqrt(24.)  # H4(y)
        losvd = gauss * (1. + h3 * h3poly + h4 * h4poly)

        galaxy = signal.fftconvolve(starNew, losvd, mode="same")  # Convolve the oversampled spectrum
        galaxy = rebin(galaxy, factor)  # Integrate spectrum into original spectral pixels
        noise = galaxy / sn  # 1sigma error spectrum
        galaxy = np.random.normal(galaxy, noise)  # Add noise to the galaxy spectrum
        start = np.array(
            [vel + np.random.random(), sigma * np.random.uniform(0.85, 1.15)]) * velscale  # Convert to km/s

        pp = ppxf(star, galaxy, noise, velscale, start, goodpixels=np.arange(dx, galaxy.size - dx),
                  plot=False, moments=4, bias=bias, quiet=True)
        result[j, :] = pp.sol

    print('Calculation time: %.2f s' % (clock() - t))

    plt.clf()
    plt.subplot(221)
    plt.plot(sigmaV * velscale, result[:, 0] - vel * velscale, '+k')
    plt.plot(sigmaV * velscale, np.ones(len(sigmaV * velscale)) * np.mean(result[:, 0] - vel * velscale), '-b')
    plt.plot(sigmaV * velscale, sigmaV * velscale * 0, '-r')
    plt.ylim(-40, 40)
    plt.xlabel('$\sigma_{in}\ (km\ s^{-1})$')
    plt.ylabel('$V - V_{in}\ (km\ s^{-1}$)')

    plt.subplot(222)
    plt.plot(sigmaV * velscale, result[:, 1] - sigmaV * velscale, '+k')
    plt.plot(sigmaV * velscale, np.ones(len(sigmaV * velscale)) * np.mean(result[:, 1] - sigmaV * velscale), '-b')
    plt.plot(sigmaV * velscale, sigmaV * velscale * 0, '-r')
    plt.ylim(-40, 40)
    plt.xlabel('$\sigma_{in}\ (km\ s^{-1})$')
    plt.ylabel('$\sigma - \sigma_{in}\ (km\ s^{-1}$)')

    plt.subplot(223)
    plt.plot(sigmaV * velscale, result[:, 2], '+k')
    plt.plot(sigmaV * velscale, sigmaV * velscale * 0 + h3, '-r')
    plt.plot(sigmaV * velscale, np.ones(len(sigmaV * velscale)) * np.mean(result[:, 2]), '-b')
    plt.ylim(-0.2 + h3, 0.2 + h3)
    plt.xlabel('$\sigma_{in}\ (km\ s^{-1})$')
    plt.ylabel('$h_3$')

    plt.subplot(224)
    plt.plot(sigmaV * velscale, result[:, 3], '+k')
    plt.plot(sigmaV * velscale, sigmaV * velscale * 0 + h4, '-r')
    plt.plot(sigmaV * velscale, np.ones(len(sigmaV * velscale)) * np.mean(result[:, 3]), '-b')
    plt.ylim(-0.2 + h4, 0.2 + h4)
    plt.xlabel('$\sigma_{in}\ (km\ s^{-1})$')
    plt.ylabel('$h_4$')

    plt.tight_layout()
    plt.show()


def scopy_flux(flux_sci, flux_scopy_fits, flux_scopy_range, flux_scopy_file):
    """
    Combine (average) all spectra (according to bin) in the image for a given spectral range, calculate mean flux
    INPUT: FLUX_SCI, FLUX_SCOPY_FITS, FLUX_SCOPY_RANGE
    OUTPUT: FLUX_SCOPY_FILE
    """

    if os.path.exists(flux_scopy_file):
        print('File {} already exists'.format(flux_scopy_file))
        return

    files_in_dir = glob.glob(flux_sci.format('*'))
    assert len(files_in_dir) > 0, 'No files match {}'.format(flux_sci.format('*'))

    from pyraf import iraf

    iraf.noao()
    iraf.onedspec()

    flux_scopy_fits_i_data_mean = []

    for i in range(len(files_in_dir)):

        flux_sci_i = flux_sci.format(i)
        flux_scopy_fits_i = flux_scopy_fits.format(i)

        if not os.path.exists(flux_scopy_fits_i):
            iraf.scopy(flux_sci_i, flux_scopy_fits_i, w1=flux_scopy_range[0], w2=flux_scopy_range[1])

        flux_scopy_fits_i_data = fits.getdata(flux_scopy_fits_i, 0)
        assert flux_scopy_fits_i_data.ndim != 0, "Scrop'd array is empty"

        flux_scopy_fits_i_data_mean.append(flux_scopy_fits_i_data.mean())

    np.array(flux_scopy_fits_i_data_mean).tofile(flux_scopy_file, sep='\n')


def fxcor(spec, task, template_spec, spec_list_file, fxcor_output_file, fxcor_range="4173-5447", interactive='no'):
    """
    Run fxcor on the binned spectra to determine the relative velocity (with respect to the template spectra)
    INPUT: BIN_SCI, FXCOR_TEMPLATE
    OUTPUT: FXCOR_BIN_LIST, FXCOR_FILE(s)

    Note: make sure to change the parameters for each function based on desired output/use.
    """

    if os.path.exists(fxcor_output_file + '.txt'):
        print('File {} already exists'.format(fxcor_output_file + '.txt'))
        return

    assert task == "abs" or task == "ems", "'task' is neight 'ab' or 'ems'"

    spec_list = []
    for i in range(len(glob.glob(spec.format('*')))):
        spec_list.append(spec.format(i))
    assert len(spec_list) > 0, 'Input files {} do not exist'.format(spec.format('*'))
    np.array(spec_list).tofile(spec_list_file, sep='\n')

    from pyraf import iraf

    iraf.noao()
    iraf.rv()

    if task == 'ems':
        iraf.fxcor('@{}'.format(spec_list_file), spec.format(template_spec), output=fxcor_output_file, continuum="both",
                   interactive=interactive, order=1, high_rej=2, low_rej=2, osample=fxcor_range, rsample=fxcor_range,
                   rebin="smallest", imupdate="no", pixcorr="no", filter="both", f_type="welch", cuton=20, cutoff=1000,
                   fullon=30, fulloff=800, ra="RA", dec="DEC", ut="UTSTART", epoch="EQUINOX", verbose="txtonly")

    elif task == 'abs':
        # Run interactively to make sure not fitting noise features, adapet osample/rsample to avoid such noise
        iraf.fxcor('@{}'.format(spec_list_file), spec.format(template_spec), output=fxcor_output_file, continuum="both",
                   interactive=interactive, order=1, high_rej=2, low_rej=2, osample=fxcor_range, rsample=fxcor_range,
                   rebin="smallest", imupdate="no", pixcorr="no", filter="both", f_type="welch", cuton=20, cutoff=800,
                   fullon=30, fulloff=800, ra="RA", dec="DEC", ut="UTSTART", epoch="EQUINOX", verbose="txtonly")

    assert os.path.exists(fxcor_output_file + '.txt'), 'Error in iraf.fxcor: File {} was not created'.format(
        fxcor_output_file + '.txt')


def rvsao(bin_sci, task, template_spectra, rvsao_file, rvsao_bin_list, interactive="no", linesig=1.5,
          czguess=0., st_lambda="INDEF", end_lambda="INDEF"):
    """
    Use the rvsao task emsao or xcsao to measure relative velocity from emission or absoption spectra
    """

    assert task == 'xcsao' or task == 'emsao', "task is not either 'xcsao' or 'emsao'"

    if os.path.exists(rvsao_file):
        print('File {} already exists'.format(rvsao_file))
        return

    bin_list = []
    for i in range(len(glob.glob(bin_sci.format('*')))):  # to ensure order is 0-61 (not 0, 1, 10, 11, etc)
        bin_list.append(bin_sci.format(i))
    assert len(bin_list) > 0, 'Absorption/emission(?) bin spectra do not exist: {}'.format(em_bin_sci.format('*'))
    np.array(bin_list).tofile(rvsao_bin_list, sep='\n')

    from pyraf import iraf

    iraf.images()
    iraf.rvsao()

    if task == 'xcsao':
        iraf.xcsao('@{}'.format(rvsao_bin_list), templates=bin_sci.format(template_spectra), report_mode=2,
                   logfiles=rvsao_file, displot=interactive, low_bin=10, top_low=20, top_nrun=80, nrun=211,
                   zeropad="yes", nzpass=1, curmode="no", pkmode=2, s_emchop="no", vel_init="guess", czguess=czguess,
                   st_lambda=st_lambda, end_lambda=end_lambda)

    elif task == 'emsao':

        # Run this interactively as the fit can fail often. Depending on S/N you may have to decrease the linesig
        # to detect the lines, or if lines are matched inconsistently play with the st_lambda, end_lambda parameters
        iraf.emsao('@{}'.format(rvsao_bin_list), logfiles=rvsao_file, displot=interactive, report_mode=2,
                   contsub_plot="no", st_lambda=st_lambda, end_lambda=end_lambda, vel_init="guess", czguess=czguess,
                   linesig=linesig)
        # linesig=0.92, vel_init="guess", czguess=1500,


def plot_velfield_setup(vel, v2b_xy_file, flux_scopy_file):
    assert os.path.exists(v2b_xy_file), 'File {} does not exist'.format(v2b_xy_file)
    assert os.path.exists(flux_scopy_file), 'File {} does not exist'.format(flux_scopy_file)
    assert len(vel) > 0, 'Input velocity does not make sense'

    xbar, ybar, xnode, ynode = np.loadtxt(v2b_xy_file, unpack=True, skiprows=1)
    flux = np.loadtxt(flux_scopy_file, unpack=True)

    assert len(xbar) == len(ybar), 'Xbar is not the same length as Ybar'
    assert len(xbar) == len(vel), 'Xbar is not the same length as vel'
    assert len(xbar) == len(flux), 'Xbar is not the same length as flux'

    plt.clf()
    plt.title('Velocity')
    plot_velfield(xbar, ybar, vel, flux=flux, colorbar=True, label='km/s')
    plt.show()


def plot_velfield_bybin(vel, v2b_file, vmin=None, vmax=None):
    """
    Read in the velocity and spatial information and plot with a colourbar
    to change range of colourbar use vmin and max parameters
    """

    if vmin is None:
        vmin = np.min(vel)

    if vmax is None:
        vmax = np.max(vel)

    bin_table = pd.read_table(v2b_file, sep=r"\s*", engine='python', names=["x", "y", "binNum"], skiprows=1)

    vel_arr = np.zeros(shape=len(bin_table.index))

    for i in range(len(vel)):
        query = bin_table.query('binNum == {}'.format(i))
        for idx in query.index:
            vel_arr[idx] = vel[i]

    x_uq = np.unique(bin_table.x.values)
    y_uq = np.unique(bin_table.y.values)

    table = pd.DataFrame(data={'x': bin_table.x.values, 'y': bin_table.y.values, 'vel': vel_arr})

    grid = np.zeros(shape=(y_uq.size, x_uq.size))
    for x in range(x_uq.size):
        for y in range(y_uq.size):
            query = table.query('x == {}'.format(x + x_uq[0])).query('y == {}'.format(y + y_uq[0]))
            if len(query.vel.values) == 0:
                print ('ERROR: coordinate ({},{}) does not have bin number'.format(x + x_uq[0], y + y_uq[0]))
                grid[y, x] = np.mean(table.vel.values)
            else:
                grid[y, x] = query.vel.values[0]

    plt.figure()
    ax = plt.gca()
    im = ax.imshow(grid, vmin=vmin, vmax=vmax)
    plt.axis('off')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.show()


def rebin(x, factor):
    """
    Rebin a one-dimensional vector by averaging
    in groups of "factor" adjacent values

    """
    return np.mean(x.reshape(-1, factor), axis=1)


def read_emsao_output(emaso_input):
    emsao_vel = []
    with open(emaso_input) as infile:
        for line in infile:
            emsao_vel.append(float(line.split()[4]))
    return emsao_vel


def read_xcsao_output(xcsao_input):
    return pd.read_table(XCSAO_FILE, sep=r"\s*", engine='python', usecols=[3], names=["vrel"], squeeze=True).values


def read_fxcor_output(fxcor_file):
    return pd.read_table(fxcor_file + '.txt', sep=r"\s*", engine='python', skiprows=16, usecols=[10],
                         names=["vrel"], squeeze=True).values


def crop_table(xysn_file, xysn_file_cropped, x=[25, 45], y=[20, 40]):
    if os.path.exists(xysn_file_cropped):
        return
    table = pd.read_table(xysn_file, sep=r"\s*", engine='python', names=["y", "x", "sig", "noise"])
    table_cropped = table.query("{} < x < {}".format(x[0], x[1])).query("{} < y < {}".format(y[0], y[1]))

    table_cropped.to_csv(xysn_file_cropped, sep=" ", header=False, index=False)


def ppxf_kinematics_gas(bin_sci, ppxf_file, ppxf_bestfit, template_fits, template_res,
                        vel_init=1500., sig_init=100., bias=0.6, plot=False, quiet=True):
    """
    Follow the pPXF usage example by Michile Cappellari
    INPUT: DIR_SCI_COMB (comb_fits_sci_{S/N}/bin_sci_{S/N}.fits), TEMPLATE_* (spectra/Mun1.30z*.fits)
    OUTPUT: PPXF_FILE (ppxf_output_sn30.txt), OUT_BESTFIT_FITS (ppxf_fits/ppxf_bestfit_sn30.fits)
    """

    if os.path.exists(ppxf_file):
        print('File {} already exists'.format(ppxf_file))
        return

    # Read in the galaxy spectrum, logarithmically rebin
    #
    assert len(glob.glob(bin_sci.format('*'))) > 0, 'Binned spectra {} not found'.format(glob.glob(bin_sci.format('*')))
    with fits.open(bin_sci.format(0)) as gal_hdu:
        gal_data = gal_hdu[0].data
        gal_hdr = gal_hdu[0].header
    lamRange1 = gal_hdr['CRVAL1'] + np.array([1. - gal_hdr['CRPIX1'], gal_hdr['NAXIS1'] - gal_hdr['CRPIX1']]) \
                                    * gal_hdr['CD1_1']

    galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_data)
    galaxy = galaxy / np.median(galaxy)  # Normalize spectrum to avoid numerical issues
    wave = np.exp(logLam1)

    # The noise level is chosen to give Chi^2/DOF=1 without regularization (REGUL=0).
    # A constant error is not a bad approximation in the fitted wavelength
    # range and reduces the noise in the fit.
    #
    noise = galaxy * 0 + 0.01528  # Assume constant noise per pixel here

    c = 299792.458  # speed of light in km/s
    FWHM_gal = 2.3  # GMOS IFU has an instrumental resolution FWHM of 2.3 A

    # stars_templates, lamRange_temp, logLam_temp = setup_spectral_library(velscale, FWHM_gal)

    # Read in template galaxy spectra
    #
    template_spectra = glob.glob(template_fits)
    assert len(template_spectra) > 0, 'Template spectra not found: {}'.format(template_fits)
    with fits.open(template_spectra[0]) as temp_hdu:
        ssp = temp_hdu[0].data
        h2 = temp_hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0., h2['CDELT1'] * (h2['NAXIS1'] - 1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    stars_templates = np.empty((sspNew.size, len(template_spectra)))
    FWHM_tem = template_res

    # FWHM_dif = np.sqrt(FWHM_gal ** 2 - template_res ** 2)
    FWHM_dif = np.sqrt(FWHM_tem ** 2 - FWHM_gal ** 2)
    sigma = FWHM_dif / 2.355 / h2['CDELT1']  # SIGMA DIFFERENCE IN PIXELS, 1.078435697220085

    for j in range(len(template_spectra)):
        with fits.open(template_spectra[j]) as temp_hdu_j:
            ssp_j = temp_hdu_j[0].data
        ssp_j = ndimage.gaussian_filter1d(ssp_j, sigma)
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp_j, velscale=velscale)
        stars_templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

    # we save the original array dimensions, which are needed to specify the regularization dimensions
    #
    reg_dim = stars_templates.shape[1:]

    # See the pPXF documentation for the keyword REGUL,
    # for an explanation of the following two lines
    #
    regul_err = 0.004  # Desired regularization error

    # Construct a set of Gaussian emission line templates.
    # Estimate the wavelength fitted range in the rest frame.
    #
    z = np.exp(vel_init / c) - 1  # Relation between velocity and redshift in pPXF
    lamRange_gal = np.array([lamRange1[0], lamRange1[-1]]) / (1 + z)
    gas_templates, line_names, line_wave = util.emission_lines(logLam2, lamRange_gal, FWHM_gal)

    # Combines the stellar and gaseous templates into a single array
    # during the PPXF fit they will be assigned a different kinematic COMPONENT value
    #
    templates = np.hstack([stars_templates, gas_templates])


    # Here the actual fit starts. The best fit is plotted on the screen.
    #
    # IMPORTANT: Ideally one would like not to use any polynomial in the fit
    # as the continuum shape contains important information on the population.
    # Unfortunately this is often not feasible, due to small calibration
    # uncertainties in the spectral shape. To avoid affecting the line strength of
    # the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
    # multiplicative ones (MDEGREE=10). This is only recommended for population, not
    # for kinematic extraction, where additive polynomials are always recommended.
    #
    # start = [vel, 180.]  # (km/s), starting guess for [V,sigma]
    start = [vel_init, sig_init]  # (km/s), starting guess for [V,sigma]

    # Assign component=0 to the stellar templates and
    # component=1 to the gas emission lines templates.
    # One can easily assign different kinematic components to different gas species
    # e.g. component=1 for the Balmer series, component=2 for the [OIII] doublet, ...)
    # Input a negative MOMENTS value to keep fixed the LOSVD of a component.
    #
    nTemps = stars_templates.shape[1]
    nLines = gas_templates.shape[1]
    component = [0] * nTemps + [1] * nLines
    moments = [4, 2]  # fit (V,sig,h3,h4) for the stars and (V,sig) for the gas
    start = [start, start]  # adopt the same starting value for both gas and stars

    bin_sci_list = glob.glob(bin_sci.format('*'))

    pp_sol = np.zeros([moments[0] + moments[1], len(bin_sci_list)])
    pp_error = np.zeros([moments[0] + moments[1], len(bin_sci_list)])

    print (len(bin_sci_list))
    # for j in range(1, len(bin_sci_list)):
    for j in range(4, 8):

        gal_data = fits.getdata(bin_sci.format(j), 0)
        gal_data_new = ndimage.gaussian_filter1d(gal_data, sigma)

        galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_data_new, velscale=velscale)
        noise = galaxy * 0 + 1  # Assume constant noise per pixel here

        # The galaxy and the template spectra do not have the same starting wavelength.
        # For this reason an extra velocity shift DV has to be applied to the template
        # to fit the galaxy spectrum. We remove this artificial shift by using the
        # keyword VSYST in the call to PPXF below, so that all velocities are
        # measured with respect to DV. This assume the redshift is negligible.
        # In the case of a high-redshift galaxy one should de-redshift its
        # wavelength to the rest frame before using the line below as described
        # in PPXF_KINEMATICS_EXAMPLE_SAURON.
        #
        c = 299792.458
        dv = (logLam2[0] - logLam1[0]) * c  # km/s

        t = clock()

        pp = ppxf(templates, galaxy, noise, velscale, start, plot=plot, moments=moments, degree=-1, mdegree=10,
                  vsyst=dv, clean=False, regul=1. / regul_err, reg_dim=reg_dim, component=component, bias=bias,
                  quiet=quiet)

        # Save the velocity, sigma, h3, h4 information for both stellar and gas to a table
        for k, sol in enumerate(pp.sol[0]):
            pp_sol[k, j] = pp.sol[0][k]
            pp_error[k, j] = pp.error[0][k]
        for k, sol in enumerate(pp.sol[1]):
            pp_sol[k + len(pp.sol[0]), j] = pp.sol[1][k]
            pp_error[k + len(pp.error[0]), j] = pp.error[1][k]

        # Plot fit results for stars and gas
        #
        if plot:
            plt.clf()
            # plt.subplot(211)
            plt.plot(wave, pp.galaxy, 'k')
            plt.plot(wave, pp.bestfit, 'b', linewidth=2)
            plt.xlabel("Observed Wavelength ($\AA$)")
            plt.ylabel("Relative Flux")
            plt.ylim([-0.1, 1.3])
            plt.xlim([np.min(wave), np.max(wave)])
            plt.plot(wave, pp.galaxy - pp.bestfit, 'd', ms=4, color='LimeGreen', mec='LimeGreen')  # fit residuals
            plt.axhline(y=-0, linestyle='--', color='k', linewidth=2)
            stars = pp.matrix[:, :nTemps].dot(pp.weights[:nTemps])
            plt.plot(wave, stars, 'r', linewidth=2)  # overplot stellar templates alone
            gas = pp.matrix[:, -nLines:].dot(pp.weights[-nLines:])
            plt.plot(wave, gas + 0.15, 'b', linewidth=2)  # overplot emission lines alone
            plt.legend()
            plt.show()

            # When the two Delta Chi^2 below are the same, the solution is the smoothest
            # consistent with the observed spectrum.
            #
            print('Desired Delta Chi^2: %.4g' % np.sqrt(2 * galaxy.size))
            print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1) * galaxy.size))
            print('Elapsed time in PPXF: %.2f s' % (clock() - t))

            w = np.where(np.array(component) == 1)[0]  # Extract weights of gas emissions
            print('++++++++++++++++++++++++++++++')
            print('Gas V=%.4g and sigma=%.2g km/s' % (pp.sol[1][0], pp.sol[1][1]))
            print('Emission lines peak intensity:')
            for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:, w].T):
                print('%12s: %.3g' % (name, weight * np.max(line)))
            print('------------------------------')

            # Plot stellar population mass distribution

            # plt.subplot(212)
            # weights = pp.weights[:np.prod(reg_dim)].reshape(reg_dim) / pp.weights.sum()
            # plt.imshow(np.rot90(weights), interpolation='nearest', cmap='gist_heat', aspect='auto', origin='upper',
            # extent=(np.log10(1.0), np.log10(17.7828), -1.9, 0.45))
            # plt.colorbar()
            # plt.title("Mass Fraction")
            # plt.xlabel("log$_{10}$ Age (Gyr)")
            # plt.ylabel("[M/H]")
            # plt.tight_layout()

            plt.legend()
            plt.show()

        if not quiet:
            print("Formal errors:")
            print("     dV    dsigma   dh3      dh4")
            print("".join("%8.2g" % f for f in pp.error * np.sqrt(pp.chi2)))

        hdu_best = fits.PrimaryHDU()
        hdu_best.data = pp.bestfit
        hdu_best.header['CD1_1'] = gal_hdr['CD1_1']
        hdu_best.header['CDELT1'] = gal_hdr['CD1_1']
        hdu_best.header['CRPIX1'] = gal_hdr['CRPIX1']
        hdu_best.header['CRVAL1'] = gal_hdr['CRVAL1']
        hdu_best.header['NAXIS1'] = pp.bestfit.size
        hdu_best.header['CTYPE1'] = 'LINEAR'  # corresponds to sampling of values above
        hdu_best.header['DC-FLAG'] = '1'  # 0 = linear, 1= log-linear sampling
        hdu_best.writeto(ppxf_bestfit.format(j), clobber=True)

    np.savetxt(ppxf_file,
               np.column_stack(
                   [pp_sol[0], pp_error[0], pp_sol[1], pp_error[1], pp_sol[2], pp_error[2],
                    pp_sol[3], pp_error[3], pp_sol[4], pp_error[4], pp_sol[5], pp_error[5]]),
               fmt=b'%10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f %10.4f  %10.4f  %10.4f  %10.4f',
               header='Stellar pop:                                          Gas: \n'
                      'vel  d_vel  sigma  d_sigma  h3  d_h3  h4  d_h4        vel  d_vel  sigma  d_sigma')


def bin_sources(sci_ext, xysn_file, v2b_file, v2b_xy_file, plot=False):
    """
    Create bins specifically about the three sources we observe in the 2d flatted image
    """

    if os.path.exists(v2b_file) & os.path.exists(v2b_xy_file):
        print ('Files {}, {} already exists'.format(v2b_file, v2b_xy_file))
        return

    data = fits.getdata(sci_ext)

    # Find coordinates of each object
    c1 = [35.84, 29.03]  # center nucleas
    c2 = [48.47, 25.32]  # off center nucleus
    c3 = [31.16, 17.48]  # blob features
    # find elliptical parameters of each object
    a1 = 5.
    b1 = 5.44
    a2 = 4.23
    b2 = 5.2
    a3 = 5.
    b3 = 3.6

    cby, cbx = data.shape
    arr = np.zeros(shape=(cby, cbx))

    # Cut background into four squares
    arr[0:int(cby / 2), 0:int(cbx / 2)] = 0.
    arr[int(cby / 2):cby, 0:int(cbx / 2)] = 6.
    arr[0:int(cby / 2), int(cbx / 2):cbx] = 2.
    arr[int(cby / 2):cby, int(cbx / 2):cbx] = 4.

    # Cut each of those squares in half
    m = float(cby) / float(cbx)
    for x in range(int(cbx / 2)):
        for y in range(int(cby / 2)):
            if x > y / m:
                arr[y, x] = 1.
    for x in range(int(cbx / 2), cbx):
        for y in range(int(cby / 2)):
            if x > (y - cby) / (-1 * m):
                arr[y, x] = 3.
    for x in range(int(cbx / 2)):
        for y in range(int(cby / 2), cby):
            if x < (y - cby) / (-1 * m):
                arr[y, x] = 7.
    for x in range(int(cbx / 2), cbx):
        for y in range(int(cby / 2), cby):
            if x < y / m:
                arr[y, x] = 5.

    # Define ellipses over the sources
    arr = bin_by_object.ellipse1(arr, c1[0], c1[1], 6. + a1, 6. + b1, 8.)
    arr = bin_by_object.ellipse2(arr, c2[0], c2[1], 3. + a2, 3. + b2, 9.)
    arr = bin_by_object.ellipse3(arr, c3[0], c3[1], 4. + a3, 4. + b3, 10.)

    xx, yy, signal, noise = np.loadtxt(xysn_file, unpack=True)

    x = []
    y = []
    classe = []
    with open(v2b_file, 'w') as outfile:
        outfile.write('# x  y  binNum\n')
        for j in range(arr.shape[0] - 1):
            for i in range(arr.shape[1] - 1):
                if not (i in xx[noise == 0.]) and not (j in yy[noise == 0.]):
                    outfile.write('{}  {}  {}\n'.format(i, j, int(arr[j, i])))
                    x.append(i)
                    y.append(j)
                    classe.append(int(arr[j, i]))

    x = np.array(x)
    y = np.array(y)
    classe = np.array(classe)

    good = np.unique(classe)
    xNode = ndimage.mean(x, labels=classe, index=good)
    yNode = ndimage.mean(y, labels=classe, index=good)

    if plot:
        weights = [xNode, yNode]
        plt.imshow(arr)
        plt.scatter(*weights)
        plt.show()

    # I have issues here creating the barycentric nodes, so instead i choose to only use the centerpoints already found
    # cvt=True
    # pixelsize=None
    # plot=True
    # quiet=True
    # wvt=True
    # xNode2, yNode2, scale, it = voronoi_2d_binning._cvt_equal_mass(x, y, signal, noise, xNode, yNode, quiet, wvt)
    # classe, xBar, yBar, sn, area = voronoi_2d_binning._compute_useful_bin_quantities(x, y, signal, noise, xNode2, yNode2, scale)

    xBar = xNode
    yBar = yNode
    np.savetxt(v2b_xy_file, np.column_stack([xBar, yBar, xNode, yNode]), header='xBar  yBar  xNode   yNode',
               fmt=b'%10.6f %10.6f %10.6f %10.6f')


if __name__ == '__main__':
    """
    See detailed information at the top of the script
    """

    if not os.path.exists(PROC_PATH):
        os.mkdir(PROC_PATH)
    if not os.path.exists(SCROP_PATH):
        os.mkdir(SCROP_PATH)
    if not os.path.exists(FXCOR_PATH):
        os.makedirs(FXCOR_PATH)
    if not os.path.exists(PPXF_PATH):
        os.mkdir(PPXF_PATH)
    if not os.path.exists(RVSAO_PATH):
        os.mkdir(RVSAO_PATH)
    if not os.path.exists(DIR_SCI_COMB):
        os.mkdir(DIR_SCI_COMB)
    if not os.path.exists(DIR_VAR_COMB):
        os.mkdir(DIR_VAR_COMB)

    # UNCOMMENT AS NEEDED (not the fanciest system I know...)

    # ## Optional:
    # To flatten 3D cube to 2D in specific wavelengths to visualize gas emission maps
    # (Don't use these flattend cubes otherwise)
    #
    # scrop_cube(IMAGE_CUBE, SCROP_RANGE, CUBE_SCROPD)
    # imcombine_flatten(CUBE_SCROPD, SCI_EXT_SCROPD, VAR_EXT_SCROPD)

    # Prepare the input: flatted the cube to 2D to read signal and noise planes, then apply voronoi binning
    #
    # flatten_cube(IMAGE_CUBE, SCI_EXT, VAR_EXT)
    # make_table(IMAGE_CUBE, SCI_EXT, VAR_EXT, XYSN_FILE)
    # voronoi_binning(XYSN_FILE, V2B_FILE, V2B_XY_FILE)

    # ## Optional: define bins around the sources instead of using voronoi binning method, built for S/N of ~ 30-50
    #
    # bin_sources(SCI_EXT, XYSN_FILE, SRCBIN_FILE, SRCBIN_XY_FILE) # ELLIPSE VALUES ARE HARDCODED
    # V2B_FILE = SRCBIN_FILE
    # V2B_XY_FILE = SRCBIN_XY_FILE

    # Combine spectra into thier respective bin, and calculate average flux of each bin (used for plotting):
    #
    # combine_spectra(V2B_FILE, IMAGE_CUBE, BIN_SCI, FLUX_SCI, BIN_VAR, FLUX_VAR)
    # scopy_flux(FLUX_SCI, FLUX_SCOPY_FITS, FLUX_SCOPY_RANGE, FLUX_SCOPY_FILE)

    # ## Optional: clean spectra of large noise features. Could define range of pixes (in log bins) to use with the
    # badPixel array below with pPXF if consistent across all spectra, or use clean=True with pPXF to remove outliers
    # from the fit (make sure to do plot=True when doing this)
    #
    # clean_spec(BIN_SCI, feat_sci, bad_region)

    # ## Do once:
    # To determine optimal penalty (BIAS) first run with BIAS=0 then preform monte carlo simulation
    # See information in ppxf.py (or the readme which comes with the ppxf download) for more details
    # The chosen penalty for IC 225 is BIAS=0.6
    #
    # PPXF_FILE_BIAS = PPXF_FILE.strip('.txt') + '_bias0.txt'
    # PPXF_BESTFIT_BIAS = PPXF_BESTFIT.strip('.fits') + '_bias0.fits'
    # ppxf_IC225_vJT.ppxf_kinematics(BIN_SCI, PPXF_FILE_BIAS, PPXF_BESTFIT_BIAS, TEMPLATE_FITS, TEMPLATE_RES, LAM_RANGE,
    #                                VEL_INIT, SIG_INIT, 0)
    # ppxf_IC225_vJT.ppxf_simulation(PPXF_BESTFIT_BIAS, LAM_RANGE, TARGET_SN, bias=0.5, spaxel=0)

    # Run pPXF to measure stellar population kinematics (from absorption features) and plot results
    # Assign bad pixel range to spurious consistent noise features from imperfect interpolation across the chip gap.
    # If there are strong noise featues which effect the fit that are not constant, use clean=yes instead.
    #
    # badPixels = np.append(np.arange(384., 395.), np.arange(984., 996.))
    # ppxf_vel = ppxf_IC225_vJT.ppxf_kinematics(BIN_SCI, PPXF_FILE, PPXF_BESTFIT, TEMPLATE_FITS, TEMPLATE_RES, LAM_RANGE,
    # VEL_INIT, SIG_INIT, bias=0.6, plot=True, badPixels=badPixels, clean=True)
    # ppxf_vel, ppxf_sig, h3, h4, ppxf_dvel, ppxf_dsig, dh3, dh4 = np.loadtxt(PPXF_FILE, unpack=True)
    # plot_velfield_setup(ppxf_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)
    #
    # ## Optional: If uses self-defined bins, plot with this method instead
    #
    # plot_velfield_bybin(ppxf_vel, V2B_FILE)

    # Make spectra of just emission and just absorption lines, necessary to have isolated absorption spectra for
    # cross-correlation fitting techniques, use plot=yes to check
    # If this does not work (for example, if the emission features are not symmetric, or too close to each other)
    # you can isolate the emission spectra with the function remove_lines.subtract_bestfit and use the bestfit as the
    # absorption spectra. This method will give you a LOG binned spectra though.
    # As the smission lines in our case are quite strong, we can just use the galaxy spectra for the xcor. You may have
    # adjsut the rvsao.emsao.lingsig parameter (default 1.5) to 1.0 (for example) to detect the lines if you're not
    # rebinning to S/N that's high enough for a 1.5sigma detection.
    #
    # remove_lines.fit_absorp_spectra(BIN_SCI, PPXF_BESTFIT, ABS_BIN_SCI, plot=True)
    # remove_lines.fit_emission_spectra(BIN_SCI, PPXF_BESTFIT, EM_BIN_SCI, plot=False)
    #
    # EM_BIN_SCI = EM_BIN_SCI.split('.fits')[0] + '_log.fits'
    ABS_BIN_SCI = PPXF_BESTFIT
    # remove_lines.subtract_besftfit(BIN_SCI, PPXF_BESTFIT, EM_BIN_SCI)

    # Run rv fxcor for stellar population velocites, then gas velocities
    # select an appropriate range to cross-correlation i.e. not including the chip gap feature, run interactively
    # to see if this is an issue or not.
    #
    # fxcor(ABS_BIN_SCI, 'abs', TEMPLATE_SPECTRA, ABS_FXCOR_BIN_LIST, ABS_FXCOR_FILE, "4337-5021", "no")
    # plot_velfield_setup(read_fxcor_output(ABS_FXCOR_FILE), V2B_XY_FILE, FLUX_SCOPY_FILE)
    # plot_velfield_bybin(read_fxcor_output(ABS_FXCOR_FILE), V2B_FILE)
    #
    # fxcor(EM_BIN_SCI, 'ems', TEMPLATE_SPECTRA, EM_FXCOR_BIN_LIST, EM_FXCOR_FILE, "4337-5021", "no")
    # plot_velfield_setup(read_fxcor_output(EM_FXCOR_FILE), V2B_XY_FILE, FLUX_SCOPY_FILE)
    # plot_velfield_bybin(read_fxcor_output(EM_FXCOR_FILE), V2B_FILE)

    # Run rvsao xcsao for stellar population velocites, then gas velocities. If the S/N isn't high enough for a 1.5 sig
    # detection in emsao, adjust parameter linesig.
    # If the emission line fitting procedure does not work, use the original spectra instead.
    # rvsao.xcsao as the option (s_emchop, t_emchpp) to remove emission or absorption lines, but I haven't had much
    # success using it. Instead I use the technique above
    # adjust st_lambda and end_lambda to avoid regions of strong noise that give spurious measurements
    #
    rvsao(ABS_BIN_SCI, 'xcsao', TEMPLATE_SPECTRA, XCSAO_FILE, XCSAO_BIN_LIST, "no", st_lambda=4500., end_lambda=5000.)
    plot_velfield_setup(read_xcsao_output(XCSAO_FILE), V2B_XY_FILE, FLUX_SCOPY_FILE)
    plot_velfield_bybin(read_xcsao_output(XCSAO_FILE), V2B_FILE)
    #
    # rvsao(EM_BIN_SCI, 'emsao', TEMPLATE_SPECTRA, EMSAO_FILE, EMSAO_BIN_LIST, "yes", 1., 1550., 4300.)
    # rvsao(BIN_SCI, 'emsao', TEMPLATE_SPECTRA, EMSAO_FILE, EMSAO_BIN_LIST, "yes", 1., 1550., 4300.)
    # plot_velfield_setup(read_emsao_output(EMSAO_FILE), V2B_XY_FILE, FLUX_SCOPY_FILE)
    # plot_velfield_bybin(read_emsao_output(EMSAO_FILE), V2B_FILE)

    # ## Optional: pPXF can also be used to measure the gas kinematics by fitting a 2 composite model, however the
    # fit of the stellar kinematics does not match what it measured by the method above, and doesn't make much sense.
    # Despite this, the gas kinematics match well the output from fxcor so they look a-okay
    #
    # PPXF_FILE = os.path.join(PPXF_PATH, 'ppxf_output_sn{}_absspec_4-8.txt'.format(TARGET_SN))
    # PPXF_BESTFIT = os.path.join(PPXF_PATH, 'bestfit_{}_absspec.fits')
    # #
    # ppxf_IC225_vJT.ppxf_kinematics_gas_parallel(BIN_SCI, PPXF_FILE, PPXF_BESTFIT, TEMPLATE_FITS, TEMPLATE_RES, VEL_INIT,
    #                                             SIG_INIT, plot=True)
    # ppxf_vel = pd.read_table(PPXF_FILE, sep=r"\s*", skiprows=1, usecols=[0], names=["vel"], squeeze=True).values
    # #
    # plot_velfield_setup(ppxf_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)
    # plot_velfield_setup(ppxf_vel_gas, V2B_XY_FILE, FLUX_SCOPY_FILE)


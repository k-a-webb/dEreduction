# written by Kristi Webb for data of IC225 taken by Bryan Miller
# similar to an earlier method developed by Gwen Rudie in IDL/cl/awk

from __future__ import print_function
from astropy.io import fits
import numpy as np
from time import clock
from voronoi_2d_binning import voronoi_2d_binning
import glob
import pandas as pd
import os
from scipy import ndimage, signal
from ppxf import ppxf
import ppxf_util as util
import glob
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from sauron_colormap import sauron
from cap_plot_velfield import plot_velfield
from scipy.optimize import curve_fit
# import specutils as su

"""
This script covers all the kinematic analysis of a reduced 3D image cube by the method in reduce_and_cubel.cl

To visualize the 2D flattened cube for a given spectral range the steps are as follows:
    - Crop the 3D cube for a given spectral range
            scrop_cube(IMAGE_CUBE, SCROP_RANGE, CUBE_SCROPD)
    - Flatten the cropped 3D cubed
            imcombine_flatten(CUBE_SCROPD, SCI_EXT_SCROPD, VAR_EXT_SCROPD)
*** These flattened images are not used in the kinematic analysis, but are useful for visual confirmation of
    relative flux levels between the nucleui and surrounding gas

The relavent nebular emission lines are:
    - OIII: 4980-5035 A (4959-5700/5368.54488701)
    - H_beta: 4883-4887 A
    - H_gamma: 4360-4362 A
The continuum spectrum is:
    - continuum: 4370-4870 A

To run pPXF:
    - Compile coordinates of signal and noise measurements for the binning
            make_table(IMAGE_CUBE, SCI_EXT, VAR_EXT, XYSN_FILE)
    - Preform voronoi 2d binning to acheive desired signal to noise
            voronoi_binning(XYSN_FILE, V2B_FILE, V2B_XY_FILE)
    - Combine spectra of the same bin into fits files
            combine_spectra(V2B_FILE, IMAGE_CUBE, BIN_SCI, FLUX_SCI, BIN_VAR, FLUX_VAR)
    - Preform pPXF on the binned spectra
            ppxf_kinematics(BIN_SCI, PPXF_FILE, PPXF_BESTFIT, LAM_RANGE, TEMPLATE_FITS, TEMPLATE_RESOLUTION, VEL_INIT, SIG_INIT)
To plot pPXF output:
    - Copy the desired spectral range and read mean flux values to be used for plotting
            scopy_flux(FLUX_SCI, FLUX_SCOPY_FITS, FLUX_SCOPY_RANGE, FLUX_SCOPY_FILE)
    - Access arrays to input into cap_plot_velfield()
            plot_ppxf_velfield(V2B_XY_FILE, PPXF_FILE, FLUX_SCOPY_FILE)

To run fxcor to preform cross correlation to determine relative velocity:
    - Create input file from list of spectral bins and select the bin to use as the template
            fxcor_bins(bin_sci, fxcor_bin_list, fxcor_template, fxcor_file)
To plot the pPXF output:
    - Access arrays to input into cap_plot_velfield()
            plot_fxcor_velfield(V2B_XY_FILE, FXCOR_FILE, FLUX_FILE)


Template spectra for fxcor/rvsoa (bin with high SN)
    SN 15 - spectra 41 has SN=22.3
    SN 30 - spectra 3 has SN=45.8
    SN 20 - spectra 22 has SN=30.4537 (using spectra 3 here gives a weird values at bin 30)

"""

# These are the files to change based on desired use
# --------------------------------------------------

DIR_PATH = '/Users/kwebb/IFU_reduction_wl'  # Working directory (where the 3D cube is)
IMAGE_CUBE = os.path.join(DIR_PATH, 'dchsteqpxbprgN20051205S0006_add_shift.fits')
TARGET_SN = 20  # REMEMBER TO CHANGE TEMPLATE SPECTRA FOR EACH S/N
TEMPLATE_SPECTRA = 6  # bin 38 at SN 20 is about in the middle of the frame

# To create 2D flattened science and variance images of specific spectral range
# This is only necessary for visual comparison of the intensity of regions for a given spectral range
SCROP_RANGE = [4360, 4362]  # wavelength range to scrop the image cube to which will be flattened
SCROP_PATH = os.path.join(DIR_PATH, 'scrop_proc')  # contains the output of all the scrop-type methods
CUBE_SCROPD = os.path.join(SCROP_PATH, 'IC225_3D_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))
SCI_EXT_SCROPD = os.path.join(SCROP_PATH, 'IC225_2D_sci_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))
VAR_EXT_SCROPD = os.path.join(SCROP_PATH, 'IC225_2D_var_{}_{}.fits'.format(SCROP_RANGE[0], SCROP_RANGE[1]))

# These variables define the chosen file structure I use to organise the output
# -----------------------------------------------------------------------------

# Flattened 3D cubes images
SCI_EXT = os.path.join(DIR_PATH, 'IC225_2D_sci_shift.fits')
VAR_EXT = os.path.join(DIR_PATH, 'IC225_2D_var_shift.fits')

PROC_PATH = os.path.join(DIR_PATH, 'doall_proc_{}_shift'.format(TARGET_SN))

# Organise output of Voronoi binning
XYSN_FILE = os.path.join(PROC_PATH, 'x_y_signal_noise.txt')  # output of make_table
V2B_FILE = os.path.join(PROC_PATH, 'v2b_output_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning
V2B_XY_FILE = os.path.join(PROC_PATH, 'v2b_output_xy_sn{}.txt'.format(TARGET_SN))  # output of voronoi binning

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
PPXF_PATH = os.path.join(PROC_PATH, 'ppxf_proc_cleaned')
PPXF_FILE = os.path.join(PPXF_PATH, 'ppxf_output_sn{}.txt'.format(TARGET_SN))
PPXF_BESTFIT = os.path.join(PPXF_PATH, 'bestfit_{}.fits')

# Organise output of fxcor
FXCOR_PATH = os.path.join(PROC_PATH, 'fxcor_proc')
FXCOR_BIN_LIST = os.path.join(FXCOR_PATH, 'bin_sci_sn{}_list.lis'.format(TARGET_SN))
FXCOR_FILE = os.path.join(FXCOR_PATH, 'fxcor_bin_sci_sn{}_tmpl{}'.format(TARGET_SN, TEMPLATE_SPECTRA))
EM_FXCOR_FILE = os.path.join(FXCOR_PATH, 'fxcor_em_bin_sci_sn{}_tmpl{}'.format(TARGET_SN, TEMPLATE_SPECTRA))
EM_FXCOR_BIN_LIST = os.path.join(FXCOR_PATH, 'em_bin_sci_sn{}_list.lis'.format(TARGET_SN))
ABS_FXCOR_FILE = os.path.join(FXCOR_PATH, 'fxcor_abs_bin_sci_sn{}_tmpl{}'.format(TARGET_SN, TEMPLATE_SPECTRA))
ABS_FXCOR_BIN_LIST = os.path.join(FXCOR_PATH, 'abs_bin_sci_sn{}_list.lis'.format(TARGET_SN))

# Organize output of rvsao
RVSAO_PATH = os.path.join(PROC_PATH, 'rvsao_proc')
XCSAO_BIN_LIST = os.path.join(RVSAO_PATH, 'abs_bin_sci_sn{}_list.lis'.format(TARGET_SN))
XCSAO_TEMPLATE = os.path.join(RVSAO_PATH, BIN_SCI.format(TEMPLATE_SPECTRA))
XCSAO_FILE = os.path.join(RVSAO_PATH, 'xcsao_bin_sci_sn{}_tmpl{}.txt'.format(TARGET_SN, TEMPLATE_SPECTRA))
SYSTEMATIC_VEL = 1320
EMSAO_FILE = os.path.join(RVSAO_PATH, 'emsao_bin_sci_sn{}_tmpl{}.txt'.format(TARGET_SN, TEMPLATE_SPECTRA))
EMSAO_BIN_LIST = os.path.join(RVSAO_PATH, 'em_bin_sci_sn{}_list.lis'.format(TARGET_SN))

# pPXF parameters
VEL_INIT = 1500.  # initial guess for velocity
SIG_INIT = 150.  # inital guess of sigma distribution
LAM_RANGE = [4186, 5369]  # wavelength range for logarithmic rebinning (full wavelength range observed [4186:5369])
# template from MILES Library spanning 3540-7410 A, with resolution spectral resolution of 2.54 A (FWHM),
# sigma~64 km/s, R~2000. From Sanchez-Blazquez, et al. (2006)
# (http://www.iac.es/proyecto/miles/pages/stellar-libraries/miles-library.php)
TEMPLATE_FITS = '/Users/kwebb/idl/cappellari/ppxf/spectra/Mun1.30z*.fits'
TEMPLATE_RESOLUTION = 2.54  # FWHM of the template spectra *** IN ANGSTROMS ***

# To create combined (averaged) spectra to determine mean flux of a specific wavelength range to plot
# YOU DO NOT NEED TO CROP THE SPECTRA HERE, it is just a good idea to check that you measure the same velocity for
# any given spectral range
FLUX_SCOPY_RANGE = [4186, 5369]  # [4370, 4870]
FLUX_SCOPY_FITS_SUFFIX = 'flux_{}.fits'  # wavelength range combined sci spectra from flux*
FLUX_SCOPY_FILE = os.path.join(PROC_PATH, 'binned_flux_{}.txt'.format(TARGET_SN))
FLUX_SCOPY_FITS = os.path.join(DIR_SCI_COMB, FLUX_SCOPY_FITS_SUFFIX)


def flatten_cube(image_cube, sci_ext, var_ext):
    """
    equivalent of last step of preprocessing
    """

    from pyraf import iraf

    iraf.images()

    iraf.imcombine(image_cube + '[sci]', sci_ext, project="yes")
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
        cube_header = cube_hdu[1].header
    crval3 = cube_header['CRVAL3']
    crpix3 = cube_header['CRPIX3']
    cd33 = cube_header['CD3_3']
    npix = cube_header['NAXIS3']

    wmin = crval3 + (1. - crpix3) * cd33
    wmax = crval3 + (npix - crpix3) * cd33
    # dwav = cd33

    assert scrop_range[0] >= wmin, 'Cannot crop spectra outside of wavelength range [{},{}]'.format(wmin, wmax)
    assert scrop_range[1] <= wmax, 'Cannot crop spectra outside of wavelength range [{},{}]'.format(wmin, wmax)

    x1 = (scrop_range[0] - crval3) / cd33 + crpix3
    x2 = (scrop_range[1] - crval3) / cd33 + crpix3

    scrop_cube1 = cube_data1[x1:x2, :, :]
    scrop_cube2 = cube_data2[x1:x2, :, :]

    hdu_out = fits.HDUList()
    hdu_out.append(fits.PrimaryHDU(data=cube_data0, header=cube_hdu[0].header))
    hdu_out.append(fits.ImageHDU(data=scrop_cube1, name='SCI'))
    hdu_out.append(fits.ImageHDU(data=scrop_cube2, name='VAR'))
    hdu_out.header = cube_header
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

    x, y, signal, noise = np.loadtxt(xysn_file, unpack=True)  # , skiprows=3)

    ''' CHECK VALIDITY OF THIS '''
    noise[noise == 0] = 0.000001

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

    if os.path.exists(bin_sci.format(1)):
        print('Binned spectra already exist')
        return

    # read in the output of the voronoi binning
    v2b_output = pd.read_table(v2b_file, sep=r"\s*", engine='python', skiprows=1, names=["x", "y", "bin"])

    # for each bin, combine the spectra of the same bin
    for i in range(np.amax(v2b_output['bin'].values) + 1):
        bin_sci_i = bin_sci.format(i)
        bin_var_i = bin_var.format(i)
        flux_sci_i = flux_sci.format(i)
        flux_var_i = flux_var.format(i)

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
            pix_sci.append(image_cube_sci[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1])
            pix_var.append(image_cube_var[k, spaxel_list["x"].values[j - 1] + 1, spaxel_list["y"].values[j - 1] + 1])

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


def ppxf_kinematics(bin_sci, ppxf_file, ppxf_bestfit, template_fits, template_resolution, lam_range=[4186, 5369],
                    vel_init=1500., sig_init=100., bias=0.):
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
    #   galaxy = galaxy/np.median(galaxy) # Normalize spectrum to avoid numerical issues
    #   noise = galaxy*0 + 0.0049           # Assume constant noise per pixel here

    galaxy, logLam1, velscale = util.log_rebin(lam_range, gal_lin)
    galaxy = galaxy / np.median(galaxy)  # Normalize spectrum to avoid numerical issues

    # Read the list of filenames from the Single Stellar Population library by Vazdekis (1999, ApJ, 513, 224). A subset
    # of the library is included for this example with permission. See http://purl.org/cappellari/software
    # for suggestions of more up-to-date stellar libraries.

    #   vazdekis = glob.glob(dir + 'Rbi1.30z*.fits')
    #   vazdekis.sort()
    #   FWHM_tem = 1.8 # Vazdekis spectra have a resolution FWHM of 1.8A.

    template_spectra = glob.glob(template_fits)
    assert len(template_spectra) > 0, 'Template spectra not found: {}'.format(template_fits)

    # Extract the wavelength range and logarithmically rebin one spectrum to the same velocity scale of the SAURON
    # galaxy spectrum, to determine the size needed for the array which will contain the template spectra.

    #   hdu = fits.open(vazdekis[0])
    #   ssp = hdu[0].data
    #   h2 = hdu[0].header
    #   lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    #   sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    #   templates = np.empty((sspNew.size,len(vazdekis)))

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
    FWHM_dif = np.sqrt(template_resolution ** 2 - FWHM_gal ** 2)
    sigma = FWHM_dif / 2.355 / h2['CDELT1']  # SIGMA DIFFERENCE IN PIXELS, 1.078435697220085

    # Logarithmically rebin the whole Mun library of spectra, and store each template as a column in the array TEMPLATES

    #   for j in range(len(vazdekis)):
    #       hdu = fits.open(vazdekis[j])
    #       ssp = hdu[0].data
    #       ssp = ndimage.gaussian_filter1d(ssp, sigma)
    #       sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp, velscale=velscale)
    #       templates[:, j] = sspNew / np.median(sspNew)  # Normalizes templates

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
        # Gwen uses goodpixel range: [range(2,195),range(235,790),range(830,890),range(930,940),range(980,1300)]
        # util.determine_goodpixels selects [range(0,201),range(230,792),range(821,896),range(925,945),range(975,1299)]

        # Here the actual fit starts. The best fit is plotted on the screen.
        # Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.

        start = [vel_init, sig_init]  # (km/s), starting guess for [V,sigma]
        t = clock()

        pp = ppxf(templates, galaxy, noise, velscale, start, goodpixels=goodPixels, plot=False, moments=4,
                  degree=4, vsyst=dv, bias=bias, quiet=False, clean=True)

        # print("Formal errors:")
        # print("     dV    dsigma   dh3      dh4")
        # print("".join("%8.2g" % f for f in pp.error * np.sqrt(pp.chi2)))

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


def fxcor(spec, template_spec, spec_list_file, fxcor_output_file):
    """
    Run fxcor on the binned spectra to determine the relative velocity (with respect to the template spectra)
    INPUT: BIN_SCI, FXCOR_TEMPLATE
    OUTPUT: FXCOR_BIN_LIST, FXCOR_FILE(s)
    """

    if os.path.exists(fxcor_output_file + '.txt'):
        print('File {} already exists'.format(fxcor_output_file + '.txt'))
        return

    spec_list = []
    for i in range(len(glob.glob(spec.format('*')))):
        spec_list.append(spec.format(i))
    assert len(spec_list) > 0, 'Input files {} do not exist'.format(spec.format('*'))
    np.array(spec_list).tofile(spec_list_file, sep='\n')

    from pyraf import iraf

    iraf.noao()
    iraf.rv()

    iraf.fxcor('@{}'.format(spec_list_file), spec.format(template_spec), output=fxcor_output_file, continuum="both",
               interactive="no", order=1, high_rej=2, low_rej=2, osample="4500-5100", rsample="4500-5100",
               rebin="smallest", imupdate="no", pixcorr="no", filter="both", f_type="welch", cuton=20, cutoff=1000,
               fullon=30, fulloff=800, ra="RA", dec="DEC", ut="UTSTART", epoch="EQUINOX", verbose="nogki")

    assert os.path.exists(fxcor_output_file + '.txt'), 'Error in iraf.fxcor: File {} was not created'.format(
        fxcor_output_file + '.txt')


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
    #   star, logLam, velscale = util.log_rebin(lamRange, ssp)

    star, logLam, velscale = util.log_rebin(lam_range, ssp)

    # The finite sampling of the observed spectrum is modeled in detail: the galaxy spectrum is obtained by oversampling
    # the actual observed spectrum to a high resolution. This represents the true spectrum, which is later resampled
    # to lower resolution to simulate the observations on the CCD. Similarly, the convolution with a well-sampled LOSVD
    # is done on the high-resolution spectrum, and later resampled to the observed resolution before fitting with PPXF.

    factor = 10  # Oversampling integer factor for an accurate convolution
    starNew = ndimage.interpolation.zoom(star, factor, order=1)  # The underlying spectrum, known at high resolution
    star = rebin(starNew, factor)  # Make sure that the observed spectrum is the integral over the pixels

    #   vel = 0.3  # velocity in *pixels* [=V(km/s)/velScale]
    #   h3 = 0.1  # Adopted G-H parameters of the LOSVD
    #   h4 = -0.1
    #   sn = 60.  # Adopted S/N of the Monte Carlo simulation
    #   m = 300  # Number of realizations of the simulation
    #   sigmaV = np.linspace(0.8, 4, m)  # Range of sigma in *pixels* [=sigma(km/s)/velScale]

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


def rebin(x, factor):
    """
    Rebin a one-dimensional vector by averaging
    in groups of "factor" adjacent values

    """
    return np.mean(x.reshape(-1, factor), axis=1)


def rvsao_xcsao(bin_sci, ppxf_bestfit, template_spectra, rvsao_file):
    """

    """

    if os.path.exists(rvsao_file):
        print('File {} already exists'.format(rvsao_file))
        return

    from pyraf import iraf

    iraf.images()
    iraf.rvsao()

    bin_files = []
    for i in range(len(glob.glob(bin_sci.format('*')))):  # Do this so order is 0-61 (not 0, 1, 10, 11, etc)
        bin_files.append(bin_sci.format(i))

    bestfit_template = ppxf_bestfit.format(template_spectra)

    bin_header = fits.getheader(bin_sci.format(template_spectra), 0)
    cd11 = fits.getval(bestfit_template, 'CD1_1')
    print('>>>>> CD1_1 is {}'.format(cd11))
    if cd11 is None:
        iraf.hedit(bestfit_template, "CD1_1", bin_header['CD1_1'], add="yes", update="yes", verify="no")
        iraf.hedit(bestfit_template, "CRVAL1", bin_header['CRVAL1'], add="yes", update="yes", verify="no")
        iraf.hedit(bestfit_template, "CRPIX1", bin_header['CRPIX1'], add="yes", update="yes", verify="no")

    for bin_spec in bin_files:
        iraf.xcsao(bin_spec, templates=bestfit_template, report_mode=2, logfiles=rvsao_file, displot='yes',
                   low_bin=10, top_low=20, top_nrun=80, nrun=211, zeropad="yes", nzpass=1, curmode="no", pkmode=2,
                   s_emchop="yes", vel_init="guess", czguess=1500)


def rvsao_xcsao2(bin_sci, template_spectra, xcsao_file, xcsao_bin_list):
    """
    """

    if os.path.exists(xcsao_file):
        print('File {} already exists'.format(xcsao_file))
        return

    bin_list = []
    for i in range(len(glob.glob(bin_sci.format('*')))):  # to ensure order is 0-61 (not 0, 1, 10, 11, etc)
        bin_list.append(bin_sci.format(i))
    assert len(bin_list) > 0, 'Absorption bin spectra do not exist: {}'.format(em_bin_sci.format('*'))
    np.array(bin_list).tofile(xcsao_bin_list, sep='\n')

    from pyraf import iraf

    iraf.images()
    iraf.rvsao()

    iraf.xcsao('@{}'.format(xcsao_bin_list), templates=bin_sci.format(template_spectra), report_mode=2,
               logfiles=xcsao_file, displot='yes', low_bin=10, top_low=20, top_nrun=80, nrun=211, zeropad="yes",
               nzpass=1, curmode="no", pkmode=2, s_emchop="no", vel_init="guess", czguess=1500)


def rvsao_emsao(em_bin_sci, em_bin_list_file, em_rvsao_file):
    """
    """

    if os.path.exists(em_rvsao_file):
        print('File {} already exists'.format(em_rvsao_file))
        return

    em_bin_list = []
    for i in range(len(glob.glob(em_bin_sci.format('*')))):  # to ensure order is 0-61 (not 0, 1, 10, 11, etc)
        em_bin_list.append(em_bin_sci.format(i))
    assert len(em_bin_list) > 0, 'Emission bin spectra do not exist: {}'.format(em_bin_sci.format('*'))
    np.array(em_bin_list).tofile(em_bin_list_file, sep='\n')

    from pyraf import iraf

    iraf.images()
    iraf.rvsao()

    iraf.emsao('@{}'.format(em_bin_list_file), logfiles=em_rvsao_file, vel_init="guess", czguess=1500, linesig=1.,
               displot="no", report_mode=2, contsub_plot="no")


def rebin_log_to_lin(bestfit_template):
    if os.path.exists(bestfit_template.strip('.fits') + '_rebin.fits'):
        print('File {} already exists'.format(bestfit_template.strip('.fits') + '_rebin.fits'))
        return

    with fits.open(bestfit_template) as bf_hdu:
        bf_header = bf_hdu[0].header
        bf_data = bf_hdu[0].data

    lambda_array = su.wavehead(bf_header, axis=1)

    bf_delta_rebin = 1.8
    lamda_rebin = np.arange(np.min(lambda_array), np.max(lambda_array), bf_delta_rebin)

    bf_data_rebin = su.rebin_spec(lambda_array, bf_data, lamda_rebin)

    plt.plot(lambda_array[0:99], bf_data[0:99], color='blue')
    plt.plot(lamda_rebin[0:99], bf_data_rebin[0:99], color='red')
    plt.show()

    # Write output fits file
    if not os.path.isfile(bestfit_template.strip('.fits') + '_rebin.fits'):
        hdu_rebin = fits.HDUList()
        hdu_rebin.append(fits.PrimaryHDU(data=bf_data_rebin, header=bf_header))
        hdu_rebin[0].header['CDELT1'] = bf_delta_rebin
        hdu_rebin[0].header['CD1_1'] = bf_delta_rebin
        hdu_rebin[0].header['CRPIX1'] = bf_header['CRPIX1']
        hdu_rebin[0].header['CRVAL1'] = lamda_rebin[0]
        hdu_rebin[0].header['DC-FLAG'] = 0
        hdu_rebin.writeto(bestfit_template.strip('.fits') + '_rebin.fits')

    return bestfit_template.strip('.fits') + '_rebin.fits'


def remove_absorp_lines(bin_sci, ppxf_bestfit, em_bin_sci):
    """
    """

    if os.path.exists(em_bin_sci.format(1)):
        print('Emission spectra already exists')
        return

    from pyraf import iraf

    iraf.images()

    for i in range(len(glob.glob(bin_sci.format('*')))):
        bin_file = bin_sci.format(i)
        assert os.path.exists(bin_file), 'Input spectra {} does not exist'.format(bin_file)
        bestfit_template = ppxf_bestfit.format(i)
        assert os.path.exists(bestfit_template), 'Bestfit spectra {} does not exist'.format(bestfit_template)
        # bestfit_template_rebin = rebin_log_to_lin(bestfit_template)
        # assert os.path.exists(bestfit_template_rebin), 'Rebinned bestfit spectra {} does not exist'.format(bestfit_template_rebin)

        # Create files with just absorption/just emission spectra (via operations with ppxf bestfit template spectra)
        em_bin_file = em_bin_sci.format(i)

        if not os.path.exists(em_bin_file):
            iraf.imarith(bin_file, '-', bestfit_template, em_bin_file)
            assert os.path.exists(em_bin_file), 'Imarith failed, {} was not created'.format(em_bin_file)
            # if not os.path.exists(bin_file_ablines):
            # iraf.imarith(bin_file, '-', bestfit_template_rebin, bin_file_ablines)
            #    assert os.path.exists(bin_file_ablines), 'Imarith failed, {} was not created'.format(bin_file_ablines)


def make_template():
    ssp_dir = '/Users/kwebb/idl/cappellari/ppxf/spectra/'
    ssp_list = ['bc2003_hr_m22_salp_ssp.145', 'bc2003_hr_m32_salp_ssp.145', 'bc2003_hr_m42_salp_ssp.130',
                'bc2003_hr_m52_salp_ssp.70', 'bc2003_hr_m62_salp_ssp.70', 'bc2003_hr_m72_salp_ssp.70']

    ssp0 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[0]), sep='  ', skiprows=6, engine='python',
                         usecols=[0, 1, 2, 3, 4, 5],
                         names=["lambda", "flux0", "flux1", "flux2", "flux3", "flux4"])
    ssp1 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[1]), sep='  ', skiprows=6, engine='python',
                         usecols=[1, 2, 3, 4, 5],
                         names=["flux5", "flux6", "flux7", "flux8", "flux9"])
    ssp2 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[2]), sep='  ', skiprows=6, engine='python',
                         usecols=[1, 2, 3, 4, 5, 6, 7],
                         names=["flux10", "flux11", "flux12", "flux13", "flux14", "flux15", "flux16"])
    ssp3 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[3]), sep='  ', skiprows=6, engine='python',
                         usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                         names=["flux17", "flux18", "flux19", "flux20", "flux21", "flux22", "flux23", "flux24",
                                "flux25", "flux26", "flux27", "flux28"])
    ssp4 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[4]), sep='  ', skiprows=6, engine='python',
                         usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                         names=["flux29", "flux30", "flux31", "flux32", "flux33", "flux34", "flux35", "flux36",
                                "flux37", "flux38", "flux39", "flux40"])
    ssp5 = pd.read_table('{}/{}'.format(ssp_dir, ssp_list[5]), sep='  ', skiprows=6, engine='python',
                         usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                         names=["flux41", "flux42", "flux43", "flux44", "flux45", "flux46", "flux47", "flux48",
                                "flux49", "flux50", "flux51", "flux52"])

    ssp = pd.concat([ssp0["flux0"], ssp0["flux1"], ssp0["flux2"], ssp0["flux3"], ssp0["flux4"],
                     ssp1, ssp2, ssp3, ssp4, ssp5], axis=1)

    return ssp


def clean_spec_30(bin_sci):
    bad_region = [27, 55]  # use splot to determine this region, the '$' will cahnge scale from wavelength to pixels

    # Use splot to remove the continuum from the spectra - 't' then '-' then 'q' then 'i' and choose output file name
    # To get the continuum (which I will use to fill in the bad region) subtract the continuum removed spectra
    # from the original spectra. i.e. (spectra) - (spectra - continuum) = continuum

    bin_spec = bin_sci.format(30)
    bin_spec_features = bin_spec.split('bin')[0] + 'features_bin' + bin_spec.split('bin')[1]  # my chosen file name
    bin_spec_orig = bin_spec.split('bin')[0] + 'orig_bin' + bin_spec.split('bin')[1]  # save original under new name

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


def read_emsao_output(emaso_input):
    emsao_vel = []
    with open(emaso_input) as infile:
        for line in infile:
            emsao_vel.append(float(line.split()[3]))
    return emsao_vel


def gauss_function(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b


def gauss_fit(data, lni):
    x = np.linspace(0, len(lni) - 2, len(lni) - 1)
    y = data[lni[0]:lni[-1]]

    amp = np.max(y) - np.min(y)
    x0 = np.argmax(y)
    sigma = 2.
    b = np.min(y)

    popt, pcov = curve_fit(gauss_function, x, y, p0=[amp, x0, sigma, b])
    gauss = gauss_function(x, popt[0], popt[1], popt[2], popt[3])

    return np.subtract(gauss, np.min(gauss))


def remove_emission_lines(bin_sci, abs_bin_sci):
    """
    """

    if os.path.exists(abs_bin_sci.format(1)):
        print('Absorption spectra already exists')
        return

    for j in range(len(glob.glob(bin_sci.format('*')))):

        with fits.open(bin_sci.format(j)) as hdu:
            odata = hdu[0].data
            ohdr = hdu[0].header

        galaxy, logLam1, velscale = util.log_rebin([4186, 5369], odata)
        emlns, lnames, lwave = util.emission_lines(logLam1, [4186, 5369], 2.3)

        # three lines will be detected in this wavelength range, Hg, Hb, and a OIII doublet
        # each line will be given it's own axis in the emlns array
        ln1 = emlns[:, 0]
        ln2 = emlns[:, 1]
        ln3 = emlns[:, 2]

        ln1i = []
        ln2i = []
        ln3i = []
        for i in range(ln1.size):
            if ln1[i] > 0:
                ln1i.append(i)
            if ln2[i] > 0:
                ln2i.append(i)
            if ln3[i] > 0:
                ln3i.append(i)

        # split up the OIII doublet
        ln3i1 = ln3i[0:int(len(ln3i) / 2)]
        ln3i2 = ln3i[int(len(ln3i) / 2):]

        # Apply a Gaussian fit to each emission line
        gaussHa = gauss_fit(odata, ln1i)
        gaussHb = gauss_fit(odata, ln2i)
        gaussOIIIa = gauss_fit(odata, ln3i1)
        gaussOIIIb = gauss_fit(odata, ln3i2)

        # Add all the spectral lines together
        ems = []
        for i in range(odata.size):
            if i in np.linspace(0, len(ln1i) - 2, len(ln1i) - 1) + ln1i[0]:
                ems.append(gaussHa[i - ln1i[0]])
            elif i in np.linspace(0, len(ln2i) - 2, len(ln2i) - 1) + ln2i[0]:
                ems.append(gaussHb[i - ln2i[0]])
            elif i in np.linspace(0, len(ln3i1) - 2, len(ln3i1) - 1) + ln3i1[0]:
                ems.append(gaussOIIIa[i - ln3i1[0]])
            elif i in np.linspace(0, len(ln3i2) - 2, len(ln3i2) - 1) + ln3i2[0]:
                ems.append(gaussOIIIb[i - ln3i2[0]])
            else:
                ems.append(0)

        absorp_spec = np.subtract(odata, ems)
        '''
        plt.plot(range(len(odata)), odata, '-b')
        plt.plot(range(len(ems)), ems, '-g')
        plt.plot(range(len(absorp_spec)), absorp_spec, '-r')
        plt.show()
        '''
        ahdu = fits.PrimaryHDU()
        ahdu.data = absorp_spec
        ahdu.header = ohdr
        ahdu.writeto(abs_bin_sci.format(j))


def select_emission_lines(logLam_temp, lamRange_gal, FWHM_gal):
    """
    Generates an array of Gaussian emission lines to be used as templates in PPXF.
    Additional lines can be easily added by editing this procedure.
    - logLam_temp is the natural log of the wavelength of the templates in Angstrom.
      logLam_temp should be the same as that of the stellar templates.
    - lamRange_gal is the estimated rest-frame fitted wavelength range
      Typically lamRange_gal = np.array([np.min(wave), np.max(wave)])/(1 + z),
      where wave is the observed wavelength of the fitted galaxy pixels and
      z is an initial very rough estimate of the galaxy redshift.
    - FWHM_gal is the instrumantal FWHM of the galaxy spectrum under study in
      Angstrom. Here it is assumed constant. It could be a function of wavelength.
    - The [OI], [OIII] and [NII] doublets are fixed at theoretical flux ratio~3.

    """
    lam = np.exp(logLam_temp)
    sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom

    # Balmer Series:      Hdelta   Hgamma    Hbeta   Halpha
    line_wave = np.array([4101.76, 4340.47, 4861.33, 6562.80])
    line_names = np.array(['Hdelta', 'Hgamma', 'Hbeta', 'Halpha'])
    #                 -----[OII]-----    -----[SII]-----
    line_wave.append([3726.03, 3728.82, 6716.47, 6730.85])
    line_names.append(['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731'])
    #                 -----[OIII]-----
    line_wave.append([4958.92, 5006.84])
    line_names.append('[OIII]5007d')  # single template for this doublet
    #                  -----[OI]-----
    line_wave.append([6363.67, 6300.30])
    line_names.append('[OI]6300d')  # single template for this doublet
    #                 -----[NII]-----
    line_wave.append([6548.03, 6583.41])
    line_names.append('[NII]6583d')  # single template for this doublet


    lam = np.exp(logLam_temp)
    sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom

    # Balmer Series:      Hdelta   Hgamma    Hbeta   Halpha
    line_wave = np.array([4101.76, 4340.47, 4861.33, 6562.80])
    line_names = np.array(['Hdelta', 'Hgamma', 'Hbeta', 'Halpha'])
    emission_lines = np.exp(-0.5*((lam[:,np.newaxis] - line_wave)/sigma)**2)

    #                 -----[OII]-----    -----[SII]-----
    lines = np.array([3726.03, 3728.82, 6716.47, 6730.85])
    names = np.array(['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[OIII]-----
    lines = np.array([4958.92, 5006.84])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.35*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    #                  -----[OI]-----
    lines = np.array([6363.67, 6300.30])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.33*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    #                 -----[NII]-----
    lines = np.array([6548.03, 6583.41])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.34*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    # Only include lines falling within the estimated fitted wavelength range.
    # This is important to avoid instabilities in the PPXF system solution
    #
    w = (line_wave > lamRange_gal[0]) & (line_wave < lamRange_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]



    # Only include lines falling within the estimated fitted wavelength range.
    # This is important to avoid instabilities in the PPXF system solution
    #
    w = (line_wave > lamRange_gal[0]) & (line_wave < lamRange_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave


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

    # To flatten 3D cube to 2D in specific wavelengths (Don't use these flattend cubes otherwise)
    # print('>>>>> Flattening 3D cube')
    # scrop_cube(IMAGE_CUBE, SCROP_RANGE, CUBE_SCROPD)
    # imcombine_flatten(CUBE_SCROPD, SCI_EXT_SCROPD, VAR_EXT_SCROPD)

    # Run pPXF
    # flatten_cube(IMAGE_CUBE, SCI_EXT, VAR_EXT)
    #   make_table(IMAGE_CUBE, SCI_EXT, VAR_EXT, XYSN_FILE)
    #   voronoi_binning(XYSN_FILE, V2B_FILE, V2B_XY_FILE)

    #   combine_spectra(V2B_FILE, IMAGE_CUBE, BIN_SCI, FLUX_SCI, BIN_VAR, FLUX_VAR)
    # Calculate average flux
    #   scopy_flux(FLUX_SCI, FLUX_SCOPY_FITS, FLUX_SCOPY_RANGE, FLUX_SCOPY_FILE)

    # To determine optimal penalty (BIAS) first run with BIAS=0 then preform monte carlo simulation
    # See information in ppxf.py (or the readme which comes with the ppxf download) for more details
    # The chosen penalty for IC 225 is BIAS=0.6
    #   ppxf_kinematics(BIN_SCI, PPXF_FILE.strip('.txt') + '_bias0.txt', PPXF_BESTFIT.strip('.fits') + '_bias0.fits',
    #                TEMPLATE_FITS, TEMPLATE_RESOLUTION, LAM_RANGE, VEL_INIT, SIG_INIT, 0)
    #   ppxf_simulation(PPXF_BESTFIT.strip('.fits') + '_bias0.fits', LAM_RANGE, TARGET_SN, bias=0.5, spaxel=0)

    ppxf_vel = ppxf_kinematics(BIN_SCI, PPXF_FILE, PPXF_BESTFIT, TEMPLATE_FITS, TEMPLATE_RESOLUTION, LAM_RANGE,
                               VEL_INIT, SIG_INIT, bias=0.6)
    # Plot results
    #   ppxf_vel, ppxf_sig, h3, h4, ppxf_dvel, ppxf_dsig, dh3, dh4 = np.loadtxt(PPXF_FILE, unpack=True)
    #   ppxf_vel, ppxf_dvel = np.loadtxt(PPXF_FILE.strip('.txt')+'_subsysvel.txt', unpack=True)
    plot_velfield_setup(ppxf_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)

    # Make spectre of just emission and just absorption lines
    #   remove_emission_lines(BIN_SCI, ABS_BIN_SCI)
    #   remove_absorp_lines(BIN_SCI, PPXF_BESTFIT, EM_BIN_SCI)

    # Run rv fxcor
    #   fxcor(BIN_SCI, TEMPLATE_SPECTRA, FXCOR_BIN_LIST, FXCOR_FILE)
    #   fxcor(EM_BIN_SCI, TEMPLATE_SPECTRA, EM_FXCOR_BIN_LIST, EM_FXCOR_FILE)
    #   fxcor(ABS_BIN_SCI, TEMPLATE_SPECTRA, ABS_FXCOR_BIN_LIST, ABS_FXCOR_FILE)

    #   fxcor_vel = pd.read_table(FXCOR_FILE + '.txt', sep=r"\s*", engine='python', skiprows=16, usecols=[10],
    #                          names=["vrel"], squeeze=True).values
    #   em_fxcor_vel = pd.read_table(EM_FXCOR_FILE + '.txt', sep=r"\s*", engine='python', skiprows=16, usecols=[10],
    #                             names=["vrel"], squeeze=True).values
    #   abs_fxcor_vel = pd.read_table(ABS_FXCOR_FILE + '.txt', sep=r"\s*", engine='python', skiprows=16, usecols=[10],
    #                              names=["vrel"], squeeze=True).values
    #   plot_velfield_setup(fxcor_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)
    #   plot_velfield_setup(em_fxcor_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)
    #   plot_velfield_setup(abs_fxcor_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)

    # Run rvsao xcsao
    #   rvsao_xcsao(ABS_BIN_SCI, PPXF_BESTFIT, TEMPLATE_SPECTRA, XCSAO_FILE)
    #   rvsao_xcsao2(ABS_BIN_SCI, TEMPLATE_SPECTRA, XCSAO_FILE, XCSAO_BIN_LIST)
    #   xcsao_vel = pd.read_table(XCSAO_FILE, sep=r"\s*", engine='python', usecols=[3], names=["vrel"], squeeze=True).values
    #   plot_velfield_setup(xcsao_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)
    #   rvsao_emsao(EM_BIN_SCI, EMSAO_BIN_LIST, EMSAO_FILE)
    #   emsao_vel = read_emsao_output(EMSAO_FILE)
    #   plot_velfield_setup(emsao_vel, V2B_XY_FILE, FLUX_SCOPY_FILE)

    """
    EXAMPLES OF THE SCRIPT OUTPUT

    make_table - x, y, signal, noise
    ----------
       0   0   0.0   0.0
       0   1   0.0   0.0
       0   2   0.0   0.0
       0   3   0.0   0.0
       0   4   0.0242297947407   0.0486673600972
       0   5   0.0287907179445   0.0496814213693
       0   6   0.0284730419517   0.0493335090578
       0   7   0.0294728595763   0.0492805354297
       0   8   0.0299422219396   0.0492736101151
       0   9   0.0302106719464   0.049401614815


    Voronoi binning terminal output at S/N 30
    -----------------------------------------
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


    v2b_output_sn30.txt
    -------------------
    # x  y  binNum
      0.000000   0.000000       24
      0.000000   1.000000       24
      0.000000   2.000000       24
      0.000000   3.000000       24
      0.000000   4.000000       24


    v2b_output_xy_sn30.txt
    ----------------------
    # xBar  yBar  xNode   yNode
     28.322327  32.169845  28.272727  32.045455
     23.398545  29.802916  23.290323  29.709677
     23.862846  35.158733  23.700000  35.166667
     27.721677  25.361104  27.711864  24.932203


    pPXF terminal output at S/N 15 - with moments =4
    ------------------------------
    Best Fit:       V     sigma        h3        h4        h5        h6
    comp. 0   1.56e+03      44.5   -0.0144  0.000161
    chi2/DOF: 170.7
    Function evaluations: 45
    Nonzero Templates:  3  /  31
    Formal errors:
         dV    dsigma   dh3      dh4
          58   1e+02    0.85     1.1
    Elapsed time in PPXF: 1.62 s


    ppxf_output_sn30.txt
    --------------------
    # velocity    sigma       h3           h4         dV          dsigma       dh3         dh4
    1555.078045   54.908761   -0.009226    0.017820   16.998730   13.021736    0.228814    0.091529
    1546.774034   50.347287   -0.030425    0.004431    9.008727   45.028757    0.094371    0.477264
    1559.821945   32.639841   -0.013985   -0.016247  126.604720  145.538259    2.047092    1.750520


    binned_cont_flux_30.txt
    -----------------------
    0.26100513339
    0.193801268935
    0.131731316447
    0.192921474576
    0.0340581573546
    0.212753847241
    0.151029586792


    xcsao_bin_sci_sn30.txt
    ----------------------
    bin_sci_0.fits   /Users/kwebb/IFU  77.21    -1.718    1.460 0.966 151.753
    bin_sci_1.fits   /Users/kwebb/IFU  76.20    -3.727    1.463 0.970 150.110
    bin_sci_10.fits  /Users/kwebb/IFU  43.02     5.536    2.440 0.943 143.557
    bin_sci_11.fits  /Users/kwebb/IFU  38.89     6.026    2.695 0.939 143.641
    bin_sci_12.fits  /Users/kwebb/IFU  33.66    12.732    3.091 0.927 142.409
    bin_sci_13.fits  /Users/kwebb/IFU  44.80     2.267    2.340 0.946 142.572

    """
__author__ = 'kwebb'

import math
from astropy.io import fits
import os

DIR_PATH = '/Users/kwebb/IFU_reduction_wl'  # Working directory (where the 3D cube is)
IMAGE_CUBE = os.path.join(DIR_PATH, 'dchsteqpxbprgN20051205S0006_add_shift.fits')
OUT_FILE = os.path.join(DIR_PATH, 'dchsteqpxbprgN20051205S0006_add_shift.fits')

# Find coordinates of each object
C1 = [35.068, 29.764]
C2 = [47.923, 25.605]
C3 = [30.758, 17.211]

# find elliptical parameters of each object
A1 = 5.
B1 = 5.44
A2 = 4.23
B2 = 5.2
A3 = 5.
B3 = 3.6


def ellipse1(array, cx, cy, a, b, value):
    y = np.linspace(cy - b + 0.2, cy + b - 0.2)
    for j in y:
        i0 = (cx + np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        i1 = (cx - np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        x = np.linspace(i0, i1)
        for i in x:

            bound1 = 21.
            bound2 = (j - 21.) / math.tan(math.radians(70.46)) + 39.5
            bound3 = (j - 21.) / math.tan(math.radians(-66.873)) + 39.5

            if (j > bound1) and (i < bound2):
                array[j, i] = value
    return array


def ellipse2(array, cx, cy, a, b, value):
    y = np.linspace(cy - b + 0.2, cy + b - 0.2)
    for j in y:
        i0 = (cx + np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        i1 = (cx - np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        x = np.linspace(i0, i1)
        for i in x:

            bound1 = 21.
            bound2 = (j - 21.) / math.tan(math.radians(70.46)) + 39.5
            bound3 = (j - 21.) / math.tan(math.radians(-66.873)) + 39.5

            if (j > bound1) and (i > bound2):
                array[j, i] = value
            elif (j < bound1) and (i > bound3):
                array[j, i] = value
    return array


def ellipse3(array, cx, cy, a, b, value):
    y = np.linspace(cy - b + 0.2, cy + b - 0.2)
    for j in y:
        i0 = (cx + np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        i1 = (cx - np.sqrt(a ** 2 * (1. - (j - cy) ** 2 / b ** 2)))
        x = np.linspace(i0, i1)
        for i in x:

            bound1 = 21.
            bound2 = (j - 21.) / math.tan(math.radians(70.46)) + 39.5
            bound3 = (j - 21.) / math.tan(math.radians(-66.873)) + 39.5

            if (j < bound1) and (i < bound3):
                array[j, i] = value
    return array


if __name__ == '__main__':

   with fits.open(IMAGE_CUBE) as hdu:
        hdu.info()
        cdata = hdu['SCI'].data

    cbz, cby, cbx = cdata.shape

    arr = np.zeros(shape=(cby, cbx))
    xbounds = [0, cbx]
    ybounds = [0, cby]

    arr = ellipse1(arr, C1[0], C1[1], 8. + A1, 8. + B1, 3.)
    arr = ellipse2(arr, C2[0], C2[1], 8. + A2, 8. + B2, 6.)
    arr = ellipse3(arr, C3[0], C3[1], 8. + A3, 8. + B3, 9.)

    arr = ellipse1(arr, C1[0], C1[1], 4. + A1, 4. + B1, 2.)
    arr = ellipse2(arr, C2[0], C2[1], 4. + A2, 4. + B2, 5.)
    arr = ellipse3(arr, C3[0], C3[1], 4. + A3, 4. + B3, 8.)

    arr = ellipse1(arr, C1[0], C1[1], A1, B1, 1.)
    arr = ellipse2(arr, C2[0], C2[1], A2, B2, 4.)
    arr = ellipse3(arr, C3[0], C3[1], A3, B3, 7.)

    # hdu=fits.PrimaryHDU(data=arr)
    # hdu.writeto('/Users/kwebb/IFU_reduction_wl/arr.fits', clobber=True)

    with open(OUT_FILE, 'w') as outfile:
        outfile.write('# x  y  binNum\n')
    for j in range(arr.shape[0]):
        for i in range(arr.shape[1]):
            outfile.write('{}  {}  {}\n'.format(i, j, arr[j, i]))
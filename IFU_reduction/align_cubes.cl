# Shift to align 3D cubes and then combine
# Original author Gwen Rudie
# Edited by Kristi Webb

# Aliases
# -------
set dd5=/Users/kwebb/IFU_reduction/proc
set caldir=dd5$calib/
set rawdir=/net/sbfmaps/Volumes/Science/bmiller/bdisk/bmiller/GN-2005B-DD-5/raw/

set mygmos=/Users/kwebb/iraf/scripts/gmos/
directory.nc=1
set stdimage=imtgmos

# Scripts
# -------
onedspec
gemini
#gemlocal
gmos
rv
task gscrspec=mygmos$gscrspec.cl
#task specx2w=mygmos$specx2w.cl
task wrbox=mygmos$wrbox.cl
task gspecshift=mygmos$gspecshift.cl
task findgaps=mygmos$findgaps.cl
task fndblocks=mygmos$fndblocks.cl
task qecorr=mygmos$qecorr.cl
task ifuproc=mygmos$ifuproc.cl
task gfreduce=mygmos$gfreduce.cl
task gfextract=mygmos$gfextract.cl
task gftransform=mygmos$gftransform.cl
task gfresponse=mygmos$gfresponse.cl
task gfskysub=mygmos$gfskysub.cl
task gfbkgsub=mygmos$gfbkgsub.cl
task gfdisplay=mygmos$gfdisplay.cl
task gscombine=mygmos$gscombine_new.cl
task chkblocks=mygmos$chkblocks.cl
task gkeywpars=mygmos$gkeywpars.cl
task gfshift=mygmos$gfshift.cl
task gfxcor=mygmos$gfxcor.cl
task gscalibrate=mygmos$gscalibrate.cl
task scrop=mygmos$scrop.cl

# DECIMAL SHIFT (006)
# -------------------

# make temporary folders to organise output
mkdir tmp_cal/; mkdir tmp_cal/6; mkdir tmp_cal/6/science; mkdir tmp_cal/6/var; mkdir tmp_cal/6/dq
mkdir tmp_cal/6/science/big; tmp_cal/6/science/new; tmp_cal/6/science/orig; tmp_cal/6/science/shift
mkdir tmp_cal/6/var/big; mkdir tmp_cal/6/var/new; mkdir tmp_cal/6/var/orig; mkdir tmp_cal/6/var/shift
mkdir tmp_cal/6/dq/big; mkdir tmp_cal/6/dq/new; mkdir tmp_cal/6/dq/orig; mkdir tmp_cal/6/dq/shift

# BEGIN WITH THE SCIENCE EXTENSIONS

# CUT up the cube into 2D planes
imslice mdcslqtexbrgN20051205S0006[sci] tmp_cal/6/science/orig/sci_0006_ 3

#### to make larger files so that the shift does not crop the image force larger image size
#### copy header information as well
#### insert the 2D images into the larger file
# mkimage new_sci_0006_1029.fits make 0 2 "88 69"
# mkheader tmp_cal/6/science/new/new_sci_0006_0001.fits tmp_cal/6/science/orig/sci_0006_001.fits append=yes
# iminsert tmp_cal/6/science/new/new_sci_0006_0001.fits tmp_cal/6/science/orig/sci_0006_001.fits \
#    tmp_cal/6/science/big/big_sci_0006_0001.fits replace mkim_coords.dat

for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0006_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0006_000"+i)//".fits",("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits","replace","mkim_coords.dat")}

# SHIFT the 2D images
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/shift/shift_sci_0006_0000"+i)//".fits",-9.5,0.75)}

# RESTACK the 2D images into a cube
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/science/shift/shift_sci_0006_0000"+i)//".fits",>> "tmp_cal/2D_sci_006.lis")}
imstack @tmp_cal_2D_sci_006.lis tmp/cal/shifted_sci_006.fits


# NOW DO THE SAME FOR THE VARIANCE

imclice mdcslqtexbrgN20051205S0006[var] tmp_cal/6/var/orig/var_0006_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0006_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0006_000"+i)//".fits",("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits",("tmp_cal/6/var/shift/shift_var_0006_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/var/shift/shift_var_0006_0000"+i)//".fits",>> "tmp_cal/2D_var_006.lis")}
imstack @tmp_cal/2D_var_006.lis tmp_cal/shifted_var_006.fits


# AND AGAIN FOR THE DQ

imslice mdcslqtexbrgN20051205S0006[dq] tmp_cal/6/dq/orig/dq_0006_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0006_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0006_000"+i)//".fits",("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/shift/shift_dq_0006_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/dq/shift/shift_dq_0006_0000"+i)//".fits",>> "tmp_cal/2D_dq_006.lis")}
imstack @tmp_cal/2D_dq_006.lis tmp_cal/shifted_dq_006.fits

# Change the DQ plane to integer type rather than real
imreplace tmp_cal/shifted_dq_006.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_006.fits tmp_cal/shifted_dq_006.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0006.fits[0] amdcslqtexbrgN20051205S0006_decimal.fits
tcopy cslqtexbrgN20051205S0006.fits[1] amdcslqtexbrgN20051205S0006_decimal.fits[1]
imcopy tmp_cal/shifted_sci_006.fits  amdcslqtexbrgN20051205S0006_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_006.fits amdcslqtexbrgN20051205S0006_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_006.fits amdcslqtexbrgN20051205S0006_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0006_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0006_decimal[0] nsciext 1 update+ verify-




# INTEGER SHIFT (006)
# -------------

# set new directories to organize output
mkdir tmp_cal/6_integer; mkdir tmp_cal/6_integer/science; mkdir tmp_cal/6_integer/var; mkdir tmp_cal/_integer6/dq
mkdir tmp_cal/6_integer/science/shift; mkdir tmp_cal/6_integer/var/shift; mkdir tmp_cal/6_integer/dq/shift


# BEGIN WITH SCIENCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits",("tmp_cal/6_integer/science/shift/shift_sci_0006_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/science/shift/shift_sci_0006_0000"+i)//".fits",>> "tmp_cal/2D_sci_006_integer.lis")}
imstack @tmp_cal/2D_sci_006_integer.lis tmp_cal/shifted_sci_006_integer.fits

# VARIANCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits",("tmp_cal/6_integer/var/shift/shift_var_0006_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/var/shift/shift_var_0006_0000"+i)//".fits",>> "tmp_cal/2D_var_006_integer.lis")}
imstack @tmp_cal/2D_var_006_integer.lis tmp_cal/shifted_var_006_integer.fits

# DQ
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits",("tmp_cal/6_integer/dq/shift/shift_dq_0006_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/dq/shift/shift_dq_0006_0000"+i)//".fits",>> "tmp_cal/2D_dq_006_integer.lis")}
imstack @tmp_cal/2D_dq_006_integer.lis tmp_cal/shifted_dq_006_integer.fits

# Change the DQ plane to integer type from real
imreplace tmp_cal/shifted_dq_006_integer.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_006_integer.fits tmp_cal/shifted_dq_006_integer.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0006.fits[0] amdcslqtexbrgN20051205S0006_integer.fits
tcopy cslqtexbrgN20051205S0006.fits[1] amdcslqtexbrgN20051205S0006_integer.fits[1]
imcopy tmp_cal/shifted_sci_006_integer.fits  amdcslqtexbrgN20051205S0006_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_006_integer.fits amdcslqtexbrgN20051205S0006_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_006_integer.fits amdcslqtexbrgN20051205S0006_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0006_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0006_integer[0] nsciext 1 update+ verify-


#### THEN REPEAT FOR OTHER SCIENCE FRAMES


# DECIMAL SHIFT (112)
# -------------------

mkdir tmp_cal/; mkdir tmp_cal/112; mkdir tmp_cal/112/science; mkdir tmp_cal/112/var; mkdir tmp_cal/112/dq
mkdir tmp_cal/112/science/big; tmp_cal/112/science/new; tmp_cal/112/science/orig; tmp_cal/112/science/shift
mkdir tmp_cal/112/var/big; mkdir tmp_cal/112/var/new; mkdir tmp_cal/112/var/orig; mkdir tmp_cal/112/var/shift
mkdir tmp_cal/112/dq/big; mkdir tmp_cal/112/dq/new; mkdir tmp_cal/112/dq/orig; mkdir tmp_cal/112/dq/shift

# SCIENCE
imslice mdcslqtexbrgN20051205S0112[sci] tmp_cal/6/science/orig/sci_0112_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/science/new/new_sci_0112_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/science/new/new_sci_0112_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0112_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/science/new/new_sci_0112_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0112_000"+i)//".fits",("tmp_cal/6/science/big/big_sci_0112_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0112_0000"+i)//".fits",("tmp_cal/6/science/shift/shift_sci_0112_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/science/shift/shift_sci_0112_0000"+i)//".fits",>> "tmp_cal/2D_sci_112.lis")}
imstack @tmp_cal_2D_sci_112.lis tmp/cal/shifted_sci_112.fits

# VARIANCE
imclice mdcslqtexbrgN20051205S0112[var] tmp_cal/6/var/orig/var_0112_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/var/new/new_var_0112_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/var/new/new_var_0112_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0112_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/var/new/new_var_0112_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0112_000"+i)//".fits",("tmp_cal/6/var/big/big_var_0112_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0112_0000"+i)//".fits",("tmp_cal/6/var/shift/shift_var_0112_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/var/shift/shift_var_0112_0000"+i)//".fits",>> "tmp_cal/2D_var_112.lis")}
imstack @tmp_cal/2D_var_112.lis tmp_cal/shifted_var_112.fits

# AND AGAIN FOR THE DQ
imslice mdcslqtexbrgN20051205S0112[dq] tmp_cal/6/dq/orig/dq_0112_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/dq/new/new_dq_0112_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/dq/new/new_dq_0112_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0112_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/dq/new/new_dq_0112_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0112_000"+i)//".fits",("tmp_cal/6/dq/big/big_dq_0112_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0112_0000"+i)//".fits",("tmp_cal/6/dq/shift/shift_dq_0112_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/dq/shift/shift_dq_0112_0000"+i)//".fits",>> "tmp_cal/2D_dq_112.lis")}
imstack @tmp_cal/2D_dq_112.lis tmp_cal/shifted_dq_112.fits

# Change the DQ plane to integer type rather than real
imreplace tmp_cal/shifted_dq_112.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_112.fits tmp_cal/shifted_dq_112.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0112.fits[0] amdcslqtexbrgN20051205S0112_decimal.fits
tcopy cslqtexbrgN20051205S0112.fits[1] amdcslqtexbrgN20051205S0112_decimal.fits[1]
imcopy tmp_cal/shifted_sci_112.fits  amdcslqtexbrgN20051205S0112_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_112.fits amdcslqtexbrgN20051205S0112_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_112.fits amdcslqtexbrgN20051205S0112_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0112_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0112_decimal[0] nsciext 1 update+ verify-

# INTEGER SHIFT (112)
# ------------------

mkdir tmp_cal/112_integer; mkdir tmp_cal/112_integer/science; mkdir tmp_cal/112_integer/var; mkdir tmp_cal/_integer112/dq
mkdir tmp_cal/112_integer/science/shift; mkdir tmp_cal/112_integer/var/shift; mkdir tmp_cal/112_integer/dq/shift

# BEGIN WITH SCIENCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0112_0000"+i)//".fits",("tmp_cal/6_integer/science/shift/shift_sci_0112_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/science/shift/shift_sci_0112_0000"+i)//".fits",>> "tmp_cal/2D_sci_112_integer.lis")}
imstack @tmp_cal/2D_sci_112_integer.lis tmp_cal/shifted_sci_112_integer.fits

# VARIANCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0112_0000"+i)//".fits",("tmp_cal/6_integer/var/shift/shift_var_0112_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/var/shift/shift_var_0112_0000"+i)//".fits",>> "tmp_cal/2D_var_112_integer.lis")}
imstack @tmp_cal/2D_var_112_integer.lis tmp_cal/shifted_var_112_integer.fits

# DQ
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0112_0000"+i)//".fits",("tmp_cal/6_integer/dq/shift/shift_dq_0112_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/dq/shift/shift_dq_0112_0000"+i)//".fits",>> "tmp_cal/2D_dq_112_integer.lis")}
imstack @tmp_cal/2D_dq_112_integer.lis tmp_cal/shifted_dq_112_integer.fits

# Change the DQ plane to integer type from real
imreplace tmp_cal/shifted_dq_112_integer.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_112_integer.fits tmp_cal/shifted_dq_112_integer.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0112.fits[0] amdcslqtexbrgN20051205S0112_integer.fits
tcopy cslqtexbrgN20051205S0112.fits[1] amdcslqtexbrgN20051205S0112_integer.fits[1]
imcopy tmp_cal/shifted_sci_112_integer.fits  amdcslqtexbrgN20051205S0112_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_112_integer.fits amdcslqtexbrgN20051205S0112_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_112_integer.fits amdcslqtexbrgN20051205S0112_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0112_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0112_integer[0] nsciext 1 update+ verify-


# DECIMAL SHIFT (121)
# -------------------

mkdir tmp_cal/; mkdir tmp_cal/121; mkdir tmp_cal/121/science; mkdir tmp_cal/121/var; mkdir tmp_cal/121/dq
mkdir tmp_cal/121/science/big; tmp_cal/121/science/new; tmp_cal/121/science/orig; tmp_cal/121/science/shift
mkdir tmp_cal/121/var/big; mkdir tmp_cal/121/var/new; mkdir tmp_cal/121/var/orig; mkdir tmp_cal/121/var/shift
mkdir tmp_cal/121/dq/big; mkdir tmp_cal/121/dq/new; mkdir tmp_cal/121/dq/orig; mkdir tmp_cal/121/dq/shift

# SCIENCE
imslice mdcslqtexbrgN20051205S0121[sci] tmp_cal/6/science/orig/sci_0121_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/science/new/new_sci_0121_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/science/new/new_sci_0121_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0121_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/science/new/new_sci_0121_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0121_000"+i)//".fits",("tmp_cal/6/science/big/big_sci_0121_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0121_0000"+i)//".fits",("tmp_cal/6/science/shift/shift_sci_0121_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/science/shift/shift_sci_0121_0000"+i)//".fits",>> "tmp_cal/2D_sci_121.lis")}
imstack @tmp_cal_2D_sci_121.lis tmp/cal/shifted_sci_121.fits

# VARIANCE
imclice mdcslqtexbrgN20051205S0121[var] tmp_cal/6/var/orig/var_0121_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/var/new/new_var_0121_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/var/new/new_var_0121_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0121_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/var/new/new_var_0121_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0121_000"+i)//".fits",("tmp_cal/6/var/big/big_var_0121_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0121_0000"+i)//".fits",("tmp_cal/6/var/shift/shift_var_0121_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/var/shift/shift_var_0121_0000"+i)//".fits",>> "tmp_cal/2D_var_121.lis")}
imstack @tmp_cal/2D_var_121.lis tmp_cal/shifted_var_121.fits

# AND AGAIN FOR THE DQ
imslice mdcslqtexbrgN20051205S0121[dq] tmp_cal/6/dq/orig/dq_0121_ 3
for (i=1;i<=1301;i+=1) {mkimage(("tmp_cal/6/dq/new/new_dq_0121_0000"+i)//".fits","make",0,2,"88 69")}
for (i=1;i<=1301;i+=1) {mkheader(("tmp_cal/6/dq/new/new_dq_0121_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0121_000"+i)//".fits")}
for (i=1;i<=1301;i+=1) {iminsert(("tmp_cal/6/dq/new/new_dq_0121_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0121_000"+i)//".fits",("tmp_cal/6/dq/big/big_dq_0121_0000"+i)//".fits","replace","mkim_coords.dat")}
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0121_0000"+i)//".fits",("tmp_cal/6/dq/shift/shift_dq_0121_0000"+i)//".fits",-9.5,0.75)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6/dq/shift/shift_dq_0121_0000"+i)//".fits",>> "tmp_cal/2D_dq_121.lis")}
imstack @tmp_cal/2D_dq_121.lis tmp_cal/shifted_dq_121.fits

# Change the DQ plane to integer type rather than real
imreplace tmp_cal/shifted_dq_121.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_121.fits tmp_cal/shifted_dq_121.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0121.fits[0] amdcslqtexbrgN20051205S0121_decimal.fits
tcopy cslqtexbrgN20051205S0121.fits[1] amdcslqtexbrgN20051205S0121_decimal.fits[1]
imcopy tmp_cal/shifted_sci_121.fits  amdcslqtexbrgN20051205S0121_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_121.fits amdcslqtexbrgN20051205S0121_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_121.fits amdcslqtexbrgN20051205S0121_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0121_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0121_decimal[0] nsciext 1 update+ verify-

# INTEGER SHIFT (121)
# ------------------

mkdir tmp_cal/121_integer; mkdir tmp_cal/121_integer/science; mkdir tmp_cal/121_integer/var; mkdir tmp_cal/_integer121/dq
mkdir tmp_cal/121_integer/science/shift; mkdir tmp_cal/121_integer/var/shift; mkdir tmp_cal/121_integer/dq/shift

# BEGIN WITH SCIENCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/science/big/big_sci_0121_0000"+i)//".fits",("tmp_cal/6_integer/science/shift/shift_sci_0121_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/science/shift/shift_sci_0121_0000"+i)//".fits",>> "tmp_cal/2D_sci_121_integer.lis")}
imstack @tmp_cal/2D_sci_121_integer.lis tmp_cal/shifted_sci_121_integer.fits

# VARIANCE
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/var/big/big_var_0121_0000"+i)//".fits",("tmp_cal/6_integer/var/shift/shift_var_0121_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/var/shift/shift_var_0121_0000"+i)//".fits",>> "tmp_cal/2D_var_121_integer.lis")}
imstack @tmp_cal/2D_var_121_integer.lis tmp_cal/shifted_var_121_integer.fits

# DQ
for (i=1;i<=1301;i+=1) {imshift(("tmp_cal/6/dq/big/big_dq_0121_0000"+i)//".fits",("tmp_cal/6_integer/dq/shift/shift_dq_0121_0000"+i)//".fits",-10,1)}
for (i=1;i<=1301;i+=1) {print(("tmp_cal/6_integer/dq/shift/shift_dq_0121_0000"+i)//".fits",>> "tmp_cal/2D_dq_121_integer.lis")}
imstack @tmp_cal/2D_dq_121_integer.lis tmp_cal/shifted_dq_121_integer.fits

# Change the DQ plane to integer type from real
imreplace tmp_cal/shifted_dq_121_integer.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_121_integer.fits tmp_cal/shifted_dq_121_integer.fits "short" verb-

# MERGE all the planes into an aligned MEF
imcopy dcslqtexbrgN20051205S0121.fits[0] amdcslqtexbrgN20051205S0121_integer.fits
tcopy cslqtexbrgN20051205S0121.fits[1] amdcslqtexbrgN20051205S0121_integer.fits[1]
imcopy tmp_cal/shifted_sci_121_integer.fits  amdcslqtexbrgN20051205S0121_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_121_integer.fits amdcslqtexbrgN20051205S0121_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_121_integer.fits amdcslqtexbrgN20051205S0121_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0121_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0121_integer[0] nsciext 1 update+ verify-



# AGAIN, CONVERT DQ from real type to integer

imreplace amdcslqtexbrgN20051205S0006_decimal[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051205S0006_decimal[DQ] amdcslqtexbrgN20051205S0006_decimal.fits[DQ,1,overwrite] "short" verb-
imreplace amdcslqtexbrgN20051205S0006_integer[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051205S0006_integer[DQ] amdcslqtexbrgN20051205S0006_integer.fits[DQ,1,overwrite] "short" verb-

imreplace amdcslqtexbrgN20051222S0108_decimal[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051222S0108_decimal[DQ] amdcslqtexbrgN20051222S0108_decimal.fits[DQ,1,overwrite] "short" verb-
imreplace amdcslqtexbrgN20051222S0108_integer[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051222S0108_integer[DQ] amdcslqtexbrgN20051222S0108_integer.fits[DQ,1,overwrite] "short" verb-

imreplace amdcslqtexbrgN20051222S0112_decimal[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051222S0112_decimal[DQ] amdcslqtexbrgN20051222S0112_decimal.fits[DQ,1,overwrite] "short" verb-
imreplace amdcslqtexbrgN20051222S0112_integer[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051222S0112_integer[DQ] amdcslqtexbrgN20051222S0112_integer.fits[DQ,1,overwrite] "short" verb-

imreplace amdcslqtexbrgN20051223S0121_decimal[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051223S0121_decimal[DQ] amdcslqtexbrgN20051223S0121_decimal.fits[DQ,1,overwrite] "short" verb-
imreplace amdcslqtexbrgN20051223S0121_integer[DQ] 1. lower=0.01 upper=INDEF
chpixtype amdcslqtexbrgN20051223S0121_integer[DQ] amdcslqtexbrgN20051223S0121_integer.fits[DQ,1,overwrite] "short" verb-


# COMBINE the data cubes

print("amdcslqtexbrgN20051205S0006_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051222S0108_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051222S0112_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051223S0121_decimal.fits",>> "calib_cubes_decimal.lis")

print("amdcslqtexbrgN20051205S0006_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051222S0108_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051222S0112_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051223S0121_integer.fits",>> "calib_cubes_integer.lis")

# First apply the decimal shifts then the integer shifts
gscombine @calib_cubes_decimal.lis IC225_cal_decimal.fits logfile="gmos.log" combine="average" lthresh=-9999 fl_vard+
gscombine @calib_cubes_integer.lis IC225_cal_integer.fits logfile="gmos.log" combine="average" lthresh=-9999 fl_vard+




















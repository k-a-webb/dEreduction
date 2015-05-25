# Shifting GN-2005B-DD-5
# Gwen Rudie

# Aliases
# -------
set dd5=/uv00/reu2/GN-2005B-DD-5/
set caldir=dd5$calib/
set rawdir=dd5$raw/

set mygmos=/uv01/reu2/iraf/scripts/
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


# ###################################################
# Shift the Cubes
# ###################################################


# 006 (DECIMAL SHIFT)

mkdir tmp_cal
mkdir tmp_cal/6
cd tmp_cal/6
mkdir science
mkdir var
mkdir dq

cd science
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../var
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../dq
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../../../

# SCIENCE

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051205S0006[sci] tmp_cal/6/science/orig/sci_0006_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_sci_0006_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/6/science/new/new_sci_0006_0001.fits tmp_cal/6/science/orig/sci_0006_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0006_000"+i)//".fits")
}
;


# Insert the 2D images into the larger file

# iminsert tmp_cal/6/science/new/new_sci_0006_0001.fits tmp_cal/6/science/orig/sci_0006_001.fits tmp_cal/6/science/big/big_sci_0006_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/6/science/new/new_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/orig/sci_0006_000"+i)//".fits",("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits",("tmp_cal/6/science/shift/shift_sci_0006_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6/science/shift/shift_sci_0006_0000"+i)//".fits",>> "tmp_cal/2D_sci_006.lis")
}
;

imstack @tmp_cal/2D_sci_006.lis tmp_cal/shifted_sci_006.fits


# VAR

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051205S0006[var] tmp_cal/6/var/orig/var_0006_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_var_0006_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/6/var/new/new_var_0006_0001.fits tmp_cal/6/var/orig/var_0006_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0006_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/6/var/new/new_var_0006_0001.fits tmp_cal/6/var/orig/var_0006_001.fits tmp_cal/6/var/big/big_var_0006_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/6/var/new/new_var_0006_0000"+i)//".fits",("tmp_cal/6/var/orig/var_0006_000"+i)//".fits",("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits",("tmp_cal/6/var/shift/shift_var_0006_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6/var/shift/shift_var_0006_0000"+i)//".fits",>> "tmp_cal/2D_var_006.lis")
}
;

imstack @tmp_cal/2D_var_006.lis tmp_cal/shifted_var_006.fits



# DQ

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051205S0006[dq] tmp_cal/6/dq/orig/dq_0006_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_dq_0006_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/6/dq/new/new_dq_0006_0001.fits tmp_cal/6/dq/orig/dq_0006_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0006_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/6/dq/new/new_dq_0006_0001.fits tmp_cal/6/dq/orig/dq_0006_001.fits tmp_cal/6/dq/big/big_dq_0006_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/6/dq/new/new_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/orig/dq_0006_000"+i)//".fits",("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits",("tmp_cal/6/dq/shift/shift_dq_0006_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6/dq/shift/shift_dq_0006_0000"+i)//".fits",>> "tmp_cal/2D_dq_006.lis")
}
;

imstack @tmp_cal/2D_dq_006.lis tmp_cal/shifted_dq_006.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_006.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_006.fits tmp_cal/shifted_dq_006.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051205S0006.fits[0] amdcslqtexbrgN20051205S0006_decimal.fits
tcopy cslqtexbrgN20051205S0006.fits[1] amdcslqtexbrgN20051205S0006_decimal.fits[1]
imcopy tmp_cal/shifted_sci_006.fits  amdcslqtexbrgN20051205S0006_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_006.fits amdcslqtexbrgN20051205S0006_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_006.fits amdcslqtexbrgN20051205S0006_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0006_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0006_decimal[0] nsciext 1 update+ verify-



# 006 (INTEGER SHIFT)

mkdir tmp_cal/6_integer
cd tmp_cal/6_integer
mkdir science
mkdir var
mkdir dq

cd science
mkdir shift

cd ../var
mkdir shift

cd ../dq
mkdir shift

cd ../../../

# SCIENCE

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/science/big/big_sci_0006_0000"+i)//".fits",("tmp_cal/6_integer/science/shift/shift_sci_0006_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6_integer/science/shift/shift_sci_0006_0000"+i)//".fits",>> "tmp_cal/2D_sci_006_integer.lis")
}
;

imstack @tmp_cal/2D_sci_006_integer.lis tmp_cal/shifted_sci_006_integer.fits


# VAR

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/var/big/big_var_0006_0000"+i)//".fits",("tmp_cal/6_integer/var/shift/shift_var_0006_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6_integer/var/shift/shift_var_0006_0000"+i)//".fits",>> "tmp_cal/2D_var_006_integer.lis")
}
;

imstack @tmp_cal/2D_var_006_integer.lis tmp_cal/shifted_var_006_integer.fits



# DQ

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/6/dq/big/big_dq_0006_0000"+i)//".fits",("tmp_cal/6_integer/dq/shift/shift_dq_0006_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/6_integer/dq/shift/shift_dq_0006_0000"+i)//".fits",>> "tmp_cal/2D_dq_006_integer.lis")
}
;

imstack @tmp_cal/2D_dq_006_integer.lis tmp_cal/shifted_dq_006_integer.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_006_integer.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_006_integer.fits tmp_cal/shifted_dq_006_integer.fits "short" verb-


# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051205S0006.fits[0] amdcslqtexbrgN20051205S0006_integer.fits
tcopy cslqtexbrgN20051205S0006.fits[1] amdcslqtexbrgN20051205S0006_integer.fits[1]
imcopy tmp_cal/shifted_sci_006_integer.fits  amdcslqtexbrgN20051205S0006_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_006_integer.fits amdcslqtexbrgN20051205S0006_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_006_integer.fits amdcslqtexbrgN20051205S0006_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051205S0006_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051205S0006_integer[0] nsciext 1 update+ verify-




# 108 (DECIMAL SHIFT)

mkdir tmp_cal/108
cd tmp_cal/108
mkdir science
mkdir var
mkdir dq

cd science
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../var
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../dq
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../../../

# SCIENCE

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0108[sci] tmp_cal/108/science/orig/sci_0108_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_sci_0108_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/108/science/new/new_sci_0108_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/108/science/new/new_sci_0108_0001.fits tmp_cal/108/science/orig/sci_0108_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/108/science/new/new_sci_0108_0000"+i)//".fits",("tmp_cal/108/science/orig/sci_0108_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/108/science/new/new_sci_0108_0001.fits tmp_cal/108/science/orig/sci_0108_001.fits tmp_cal/108/science/big/big_sci_0108_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/108/science/new/new_sci_0108_0000"+i)//".fits",("tmp_cal/108/science/orig/sci_0108_000"+i)//".fits",("tmp_cal/108/science/big/big_sci_0108_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/science/big/big_sci_0108_0000"+i)//".fits",("tmp_cal/108/science/shift/shift_sci_0108_0000"+i)//".fits",-5.875,-2.625)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108/science/shift/shift_sci_0108_0000"+i)//".fits",>> "tmp_cal/2D_sci_108.lis")
}
;

imstack @tmp_cal/2D_sci_108.lis tmp_cal/shifted_sci_108.fits


# VAR

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0108[var] tmp_cal/108/var/orig/var_0108_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_var_0108_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/108/var/new/new_var_0108_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/108/var/new/new_var_0108_0001.fits tmp_cal/108/var/orig/var_0108_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/108/var/new/new_var_0108_0000"+i)//".fits",("tmp_cal/108/var/orig/var_0108_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/108/var/new/new_var_0108_0001.fits tmp_cal/108/var/orig/var_0108_001.fits tmp_cal/108/var/big/big_var_0108_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/108/var/new/new_var_0108_0000"+i)//".fits",("tmp_cal/108/var/orig/var_0108_000"+i)//".fits",("tmp_cal/108/var/big/big_var_0108_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/var/big/big_var_0108_0000"+i)//".fits",("tmp_cal/108/var/shift/shift_var_0108_0000"+i)//".fits",-5.875,-2.625)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108/var/shift/shift_var_0108_0000"+i)//".fits",>> "tmp_cal/2D_var_108.lis")
}
;

imstack @tmp_cal/2D_var_108.lis tmp_cal/shifted_var_108.fits



# DQ

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0108[dq] tmp_cal/108/dq/orig/dq_0108_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_dq_0108_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/108/dq/new/new_dq_0108_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/108/dq/new/new_dq_0108_0001.fits tmp_cal/108/dq/orig/dq_0108_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/108/dq/new/new_dq_0108_0000"+i)//".fits",("tmp_cal/108/dq/orig/dq_0108_000"+i)//".fits")
}
;


# Insert the 2D images into the larger file

# iminsert tmp_cal/108/dq/new/new_dq_0108_0001.fits tmp_cal/108/dq/orig/dq_0108_001.fits tmp_cal/108/dq/big/big_dq_0108_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/108/dq/new/new_dq_0108_0000"+i)//".fits",("tmp_cal/108/dq/orig/dq_0108_000"+i)//".fits",("tmp_cal/108/dq/big/big_dq_0108_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/dq/big/big_dq_0108_0000"+i)//".fits",("tmp_cal/108/dq/shift/shift_dq_0108_0000"+i)//".fits",-5.875,-2.625)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108/dq/shift/shift_dq_0108_0000"+i)//".fits",>> "tmp_cal/2D_dq_108.lis")
}
;

imstack @tmp_cal/2D_dq_108.lis tmp_cal/shifted_dq_108.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_108.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_108.fits tmp_cal/shifted_dq_108.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051222S0108.fits[0] amdcslqtexbrgN20051222S0108_decimal.fits
tcopy cslqtexbrgN20051222S0108.fits[1] amdcslqtexbrgN20051222S0108_decimal.fits[1]
imcopy tmp_cal/shifted_sci_108.fits  amdcslqtexbrgN20051222S0108_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_108.fits amdcslqtexbrgN20051222S0108_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_108.fits amdcslqtexbrgN20051222S0108_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051222S0108_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051222S0108_decimal[0] nsciext 1 update+ verify-


# 108 (INTEGER SHIFT)

mkdir tmp_cal/108_integer
cd tmp_cal/108_integer
mkdir science
mkdir var
mkdir dq

cd science
mkdir shift

cd ../var
mkdir shift

cd ../dq
mkdir shift

cd ../../../

# SCIENCE

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/science/big/big_sci_0108_0000"+i)//".fits",("tmp_cal/108_integer/science/shift/shift_sci_0108_0000"+i)//".fits",-6,-3)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108_integer/science/shift/shift_sci_0108_0000"+i)//".fits",>> "tmp_cal/2D_sci_108_integer.lis")
}
;

imstack @tmp_cal/2D_sci_108_integer.lis tmp_cal/shifted_sci_108_integer.fits


# VAR


# Copy the header information from the old file onto the new file

# mkheader tmp_cal/108/var/new/new_var_0108_0001.fits tmp_cal/108/var/orig/var_0108_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/108/var/new/new_var_0108_0000"+i)//".fits",("tmp_cal/108/var/orig/var_0108_000"+i)//".fits")
}
;


# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/var/big/big_var_0108_0000"+i)//".fits",("tmp_cal/108_integer/var/shift/shift_var_0108_0000"+i)//".fits",-6,-3)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108_integer/var/shift/shift_var_0108_0000"+i)//".fits",>> "tmp_cal/2D_var_108_integer.lis")
}
;

imstack @tmp_cal/2D_var_108_integer.lis tmp_cal/shifted_var_108_integer.fits



# DQ

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/108/dq/big/big_dq_0108_0000"+i)//".fits",("tmp_cal/108_integer/dq/shift/shift_dq_0108_0000"+i)//".fits",-6,-3)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/108_integer/dq/shift/shift_dq_0108_0000"+i)//".fits",>> "tmp_cal/2D_dq_108_integer.lis")
}
;

imstack @tmp_cal/2D_dq_108_integer.lis tmp_cal/shifted_dq_108_integer.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_108_integer.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_108_integer.fits tmp_cal/shifted_dq_108_integer.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051222S0108.fits[0] amdcslqtexbrgN20051222S0108_integer.fits
tcopy cslqtexbrgN20051222S0108.fits[1] amdcslqtexbrgN20051222S0108_integer.fits[1]
imcopy tmp_cal/shifted_sci_108_integer.fits  amdcslqtexbrgN20051222S0108_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_108_integer.fits amdcslqtexbrgN20051222S0108_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_108_integer.fits amdcslqtexbrgN20051222S0108_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051222S0108_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051222S0108_integer[0] nsciext 1 update+ verify-



# 112 (DECIMAL SHIFT)

mkdir tmp_cal/112
cd tmp_cal/112
mkdir science
mkdir var
mkdir dq

cd science
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../var
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../dq
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../../../

# SCIENCE

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0112[sci] tmp_cal/112/science/orig/sci_0112_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_sci_0112_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/112/science/new/new_sci_0112_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/112/science/new/new_sci_0112_0001.fits tmp_cal/112/science/orig/sci_0112_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/112/science/new/new_sci_0112_0000"+i)//".fits",("tmp_cal/112/science/orig/sci_0112_000"+i)//".fits")
}
;

# Insert the tmp_cal/2D images into the larger file

# iminsert tmp_cal/112/science/new/new_sci_0112_0001.fits tmp_cal/112/science/orig/sci_0112_001.fits tmp_cal/112/science/big/big_sci_0112_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/112/science/new/new_sci_0112_0000"+i)//".fits",("tmp_cal/112/science/orig/sci_0112_000"+i)//".fits",("tmp_cal/112/science/big/big_sci_0112_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/science/big/big_sci_0112_0000"+i)//".fits",("tmp_cal/112/science/shift/shift_sci_0112_0000"+i)//".fits",-9.375,-1.125)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112/science/shift/shift_sci_0112_0000"+i)//".fits",>> "tmp_cal/2D_sci_112.lis")
}
;

imstack @tmp_cal/2D_sci_112.lis tmp_cal/shifted_sci_112.fits


# VAR

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0112[var] tmp_cal/112/var/orig/var_0112_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_var_0112_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/112/var/new/new_var_0112_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/112/var/new/new_var_0112_0001.fits tmp_cal/112/var/orig/var_0112_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/112/var/new/new_var_0112_0000"+i)//".fits",("tmp_cal/112/var/orig/var_0112_000"+i)//".fits")
}
;


# Insert the 2D images into the larger file

# iminsert tmp_cal/112/var/new/new_var_0112_0001.fits tmp_cal/112/var/orig/var_0112_001.fits tmp_cal/112/var/big/big_var_0112_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/112/var/new/new_var_0112_0000"+i)//".fits",("tmp_cal/112/var/orig/var_0112_000"+i)//".fits",("tmp_cal/112/var/big/big_var_0112_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/var/big/big_var_0112_0000"+i)//".fits",("tmp_cal/112/var/shift/shift_var_0112_0000"+i)//".fits",-9.375,-1.125)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112/var/shift/shift_var_0112_0000"+i)//".fits",>> "tmp_cal/2D_var_112.lis")
}
;

imstack @tmp_cal/2D_var_112.lis tmp_cal/shifted_var_112.fits



# DQ

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051222S0112[dq] tmp_cal/112/dq/orig/dq_0112_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_dq_0112_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/112/dq/new/new_dq_0112_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/112/dq/new/new_dq_0112_0001.fits tmp_cal/112/dq/orig/dq_0112_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/112/dq/new/new_dq_0112_0000"+i)//".fits",("tmp_cal/112/dq/orig/dq_0112_000"+i)//".fits")
}
;


# Insert the 2D images into the larger file

# iminsert tmp_cal/112/dq/new/new_dq_0112_0001.fits tmp_cal/112/dq/orig/dq_0112_001.fits tmp_cal/112/dq/big/big_dq_0112_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/112/dq/new/new_dq_0112_0000"+i)//".fits",("tmp_cal/112/dq/orig/dq_0112_000"+i)//".fits",("tmp_cal/112/dq/big/big_dq_0112_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/dq/big/big_dq_0112_0000"+i)//".fits",("tmp_cal/112/dq/shift/shift_dq_0112_0000"+i)//".fits",-9.375,-1.125)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112/dq/shift/shift_dq_0112_0000"+i)//".fits",>> "tmp_cal/2D_dq_112.lis")
}
;

imstack @tmp_cal/2D_dq_112.lis tmp_cal/shifted_dq_112.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_112.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_112.fits tmp_cal/shifted_dq_112.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051222S0112.fits[0] amdcslqtexbrgN20051222S0112_decimal.fits
tcopy cslqtexbrgN20051222S0112.fits[1] amdcslqtexbrgN20051222S0112_decimal.fits[1]
imcopy tmp_cal/shifted_sci_112.fits  amdcslqtexbrgN20051222S0112_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_112.fits amdcslqtexbrgN20051222S0112_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_112.fits amdcslqtexbrgN20051222S0112_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051222S0112_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051222S0112_decimal[0] nsciext 1 update+ verify-



# 112 (INTEGER SHIFT)

mkdir tmp_cal/112_integer
cd tmp_cal/112_integer
mkdir science
mkdir var
mkdir dq

cd science
mkdir shift

cd ../var
mkdir shift

cd ../dq
mkdir shift

cd ../../../

# SCIENCE

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/science/big/big_sci_0112_0000"+i)//".fits",("tmp_cal/112_integer/science/shift/shift_sci_0112_0000"+i)//".fits",-9,-1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112_integer/science/shift/shift_sci_0112_0000"+i)//".fits",>> "tmp_cal/2D_sci_112_integer.lis")
}
;

imstack @tmp_cal/2D_sci_112_integer.lis tmp_cal/shifted_sci_112_integer.fits


# VAR

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/var/big/big_var_0112_0000"+i)//".fits",("tmp_cal/112_integer/var/shift/shift_var_0112_0000"+i)//".fits",-9,-1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112_integer/var/shift/shift_var_0112_0000"+i)//".fits",>> "tmp_cal/2D_var_112_integer.lis")
}
;

imstack @tmp_cal/2D_var_112_integer.lis tmp_cal/shifted_var_112_integer.fits



# DQ

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/112/dq/big/big_dq_0112_0000"+i)//".fits",("tmp_cal/112_integer/dq/shift/shift_dq_0112_0000"+i)//".fits",-9,-1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/112_integer/dq/shift/shift_dq_0112_0000"+i)//".fits",>> "tmp_cal/2D_dq_112_integer.lis")
}
;

imstack @tmp_cal/2D_dq_112_integer.lis tmp_cal/shifted_dq_112_integer.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_112_integer.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_112_integer.fits tmp_cal/shifted_dq_112_integer.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051222S0112.fits[0] amdcslqtexbrgN20051222S0112_integer.fits
tcopy cslqtexbrgN20051222S0112.fits[1] amdcslqtexbrgN20051222S0112_integer.fits[1]
imcopy tmp_cal/shifted_sci_112_integer.fits  amdcslqtexbrgN20051222S0112_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_112_integer.fits amdcslqtexbrgN20051222S0112_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_112_integer.fits amdcslqtexbrgN20051222S0112_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051222S0112_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051222S0112_integer[0] nsciext 1 update+ verify-




# 121 (DECIMAL SHIFT)

mkdir tmp_cal/121
cd tmp_cal/121
mkdir science
mkdir var
mkdir dq

cd science
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../var
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../dq
mkdir big
mkdir new
mkdir orig
mkdir shift

cd ../../../

# SCIENCE

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051223S0121[sci] tmp_cal/121/science/orig/sci_0121_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_sci_0121_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/121/science/new/new_sci_0121_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/121/science/new/new_sci_0121_0001.fits tmp_cal/112/science/orig/sci_0121_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/121/science/new/new_sci_0121_0000"+i)//".fits",("tmp_cal/121/science/orig/sci_0121_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/121/science/new/new_sci_0121_0001.fits tmp_cal/121/science/orig/sci_0121_001.fits tmp_cal/121/science/big/big_sci_0121_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/121/science/new/new_sci_0121_0000"+i)//".fits",("tmp_cal/121/science/orig/sci_0121_000"+i)//".fits",("tmp_cal/121/science/big/big_sci_0121_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/science/big/big_sci_0121_0000"+i)//".fits",("tmp_cal/121/science/shift/shift_sci_0121_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121/science/shift/shift_sci_0121_0000"+i)//".fits",>> "tmp_cal/2D_sci_121.lis")
}
;

imstack @tmp_cal/2D_sci_121.lis tmp_cal/shifted_sci_121.fits


# VAR

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051223S0121[var] tmp_cal/121/var/orig/var_0121_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_var_0121_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/121/var/new/new_var_0121_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/121/var/new/new_var_0121_0001.fits tmp_cal/121/var/orig/var_0121_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/121/var/new/new_var_0121_0000"+i)//".fits",("tmp_cal/121/var/orig/var_0121_000"+i)//".fits")
}
;


# Insert the 2D images into the larger file

# iminsert tmp_cal/121/var/new/new_var_0121_0001.fits tmp_cal/121/var/orig/var_0121_001.fits tmp_cal/121/var/big/big_var_0121_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/121/var/new/new_var_0121_0000"+i)//".fits",("tmp_cal/121/var/orig/var_0121_000"+i)//".fits",("tmp_cal/121/var/big/big_var_0121_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/var/big/big_var_0121_0000"+i)//".fits",("tmp_cal/121/var/shift/shift_var_0121_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121/var/shift/shift_var_0121_0000"+i)//".fits",>> "tmp_cal/2D_var_121.lis")
}
;

imstack @tmp_cal/2D_var_121.lis tmp_cal/shifted_var_121.fits



# DQ

# Cut up the cube into 2D planes

imslice mdcslqtexbrgN20051223S0121[dq] tmp_cal/121/dq/orig/dq_0121_ 3

# Make larger files so that the shift doesnt crop the image

# mkimage new_dq_0121_1029.fits make 0 2 "88 69"
# "make" for make a new file, "0" for the value of the pixels "2" for 2D "'88 69'" are the dimensions

for (i=1;i<=1301;i+=1) {
mkimage(("tmp_cal/121/dq/new/new_dq_0121_0000"+i)//".fits","make",0,2,"88 69")
}
;

# Copy the header information from the old file onto the new file

# mkheader tmp_cal/121/dq/new/new_dq_0121_0001.fits tmp_cal/121/dq/orig/dq_0121_001.fits append=yes

for (i=1;i<=1301;i+=1) {
mkheader(("tmp_cal/121/dq/new/new_dq_0121_0000"+i)//".fits",("tmp_cal/121/dq/orig/dq_0121_000"+i)//".fits")
}
;

# Insert the 2D images into the larger file

# iminsert tmp_cal/121/dq/new/new_dq_0121_0001.fits tmp_cal/121/dq/orig/dq_0121_001.fits tmp_cal/121/dq/big/big_dq_0121_0001.fits replace mkim_coords.dat 

for (i=1;i<=1301;i+=1) {
iminsert(("tmp_cal/121/dq/new/new_dq_0121_0000"+i)//".fits",("tmp_cal/121/dq/orig/dq_0121_000"+i)//".fits",("tmp_cal/121/dq/big/big_dq_0121_0000"+i)//".fits","replace","mkim_coords.dat")
}
;

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/dq/big/big_dq_0121_0000"+i)//".fits",("tmp_cal/121/dq/shift/shift_dq_0121_0000"+i)//".fits",-9.5,0.75)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121/dq/shift/shift_dq_0121_0000"+i)//".fits",>> "tmp_cal/2D_dq_121.lis")
}
;

imstack @tmp_cal/2D_dq_121.lis tmp_cal/shifted_dq_121.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_121.fits 1. lower=0.01 upper=INDEF 
chpixtype tmp_cal/shifted_dq_121.fits tmp_cal/shifted_dq_121.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051223S0121.fits[0] amdcslqtexbrgN20051223S0121_decimal.fits
tcopy cslqtexbrgN20051223S0121.fits[1] amdcslqtexbrgN20051223S0121_decimal.fits[1]
imcopy tmp_cal/shifted_sci_121.fits  amdcslqtexbrgN20051223S0121_decimal[SCI,1,append+]
imcopy tmp_cal/shifted_var_121.fits amdcslqtexbrgN20051223S0121_decimal[VAR,1,append+]
imcopy tmp_cal/shifted_dq_121.fits amdcslqtexbrgN20051223S0121_decimal[DQ,1,append+]
hedit amdcslqtexbrgN20051223S0121_decimal[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051223S0121_decimal[0] nsciext 1 update+ verify-



# 121 (INTEGER SHIFT)

mkdir tmp_cal/121_integer
cd tmp_cal/121_integer
mkdir science
mkdir var
mkdir dq

cd science
mkdir shift

cd ../var
mkdir shift

cd ../dq
mkdir shift

cd ../../../

# SCIENCE

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/science/big/big_sci_0121_0000"+i)//".fits",("tmp_cal/121_integer/science/shift/shift_sci_0121_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121_integer/science/shift/shift_sci_0121_0000"+i)//".fits",>> "tmp_cal/2D_sci_121_integer.lis")
}
;

imstack @tmp_cal/2D_sci_121_integer.lis tmp_cal/shifted_sci_121_integer.fits


# VAR

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/var/big/big_var_0121_0000"+i)//".fits",("tmp_cal/121_integer/var/shift/shift_var_0121_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121_integer/var/shift/shift_var_0121_0000"+i)//".fits",>> "tmp_cal/2D_var_121_integer.lis")
}
;

imstack @tmp_cal/2D_var_121_integer.lis tmp_cal/shifted_var_121_integer.fits



# DQ

# Shift the 2D images

for (i=1;i<=1301;i+=1) {
imshift(("tmp_cal/121/dq/big/big_dq_0121_0000"+i)//".fits",("tmp_cal/121_integer/dq/shift/shift_dq_0121_0000"+i)//".fits",-10,1)
}
;

# Restack the 2D images into a cube

for (i=1;i<=1301;i+=1) {
print(("tmp_cal/121_integer/dq/shift/shift_dq_0121_0000"+i)//".fits",>> "tmp_cal/2D_dq_121_integer.lis")
}
;

imstack @tmp_cal/2D_dq_121_integer.lis tmp_cal/shifted_dq_121_integer.fits

# Changing the DQ plane to integer rather than real format

imreplace tmp_cal/shifted_dq_121_integer.fits 1. lower=0.01 upper=INDEF
chpixtype tmp_cal/shifted_dq_121_integer.fits tmp_cal/shifted_dq_121_integer.fits "short" verb-

# Stuff all the planes together into one aligned MEF


imcopy dcslqtexbrgN20051223S0121.fits[0] amdcslqtexbrgN20051223S0121_integer.fits
tcopy cslqtexbrgN20051223S0121.fits[1] amdcslqtexbrgN20051223S0121_integer.fits[1]
imcopy tmp_cal/shifted_sci_121_integer.fits  amdcslqtexbrgN20051223S0121_integer[SCI,1,append+]
imcopy tmp_cal/shifted_var_121_integer.fits amdcslqtexbrgN20051223S0121_integer[VAR,1,append+]
imcopy tmp_cal/shifted_dq_121_integer.fits amdcslqtexbrgN20051223S0121_integer[DQ,1,append+]
hedit amdcslqtexbrgN20051223S0121_integer[0] nextend 4 update+ verify-
hedit amdcslqtexbrgN20051223S0121_integer[0] nsciext 1 update+ verify-

# Make DQ integer values instead of real

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


# Combine the Data Cubes

print("amdcslqtexbrgN20051205S0006_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051222S0108_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051222S0112_decimal.fits",>> "calib_cubes_decimal.lis")
print("amdcslqtexbrgN20051223S0121_decimal.fits",>> "calib_cubes_decimal.lis")

print("amdcslqtexbrgN20051205S0006_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051222S0108_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051222S0112_integer.fits",>> "calib_cubes_integer.lis")
print("amdcslqtexbrgN20051223S0121_integer.fits",>> "calib_cubes_integer.lis")

# First the decimal shifts

gscombine @calib_cubes_decimal.lis IC225_cal_decimal.fits logfile="gmos.log" combine="average" lthresh=-9999 fl_vard+

# The integer shifts

gscombine @calib_cubes_integer.lis IC225_cal_integer.fits logfile="gmos.log" combine="average" lthresh=-9999 fl_vard+




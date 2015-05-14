# GFCUBE of GN-2005B-DD-5
# Origanal author Gwen Rudie
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


# GFCUBE - reorganize MEF extensions into a 3D cube
gfcube cslqtexbrgN20051205S0006.fits
#gfcube cslqtexbrgN20051222S0108.fits
#gfcube cslqtexbrgN20051222S0112.fits
#gfcube cslqtexbrgN20051223S0121.fits

# create a dummy variable plane for gfcube
imcopy cslqtexbrgN20051205S0006.fits[0] tmp_var0006.fits
tcopy slqtexbrgN20051205S0006.fits[1] tmp_var0006.fits[1]
imcopy cslqtexbrgN20051205S0006.fits[var,1] tmp_var0006.fits[SCI,1,append+]
hedit tmp_var0006.fits[0] nextend 2 update+ verify-
hedit tmp_var0006.fits[0] nsciext 1 update+ verify-
gfcube tmp_var0006.fits

# stack the outputs from pfcube into a single MDF file
imcopy dcslqtexbrgN20051205S0006.fits[0] mdcslqtexbrgN20051205S0006.fits
tcopy cslqtexbrgN20051205S0006.fits[1] mdcslqtexbrgN20051205S0006.fits[1]
imcopy dcslqtexbrgN20051205S0006.fits[sci,1] mdcslqtexbrgN20051205S0006[SCI,1,append+]
imcopy dtmp_var0006.fits[sci,1] mdcslqtexbrgN20051205S0006[VAR,1,append+]
imcopy dtmp_dq0006.fits[sci,1] mdcslqtexbrgN20051205S0006[DQ,1,append+]
hedit mdcslqtexbrgN20051205S0006[0] nextend 4 update+ verify-
hedit mdcslqtexbrgN20051205S0006[0] nsciext 1 update+ verify-


#### REPEAT FOR EACH FILE
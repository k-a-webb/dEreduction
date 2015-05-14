# Edited Example IFU reduction (example written by Bryan Miller)
# Kristi Webb

# 2005dec03
#GN-2005B-DD-5-3    16:20:39    1            239           Twilight               g        B600/478.0  IFU             300            full   
#GN-2005B-DD-5-3    not quite enough counts, redo junk
#GN-2005B-DD-5-3    16:29:20    2            240           Twilight               g        B600/478.0  IFU             400            full   
#GN-2005B-DD-5-3    35000 cts max ok
#GN-2005B-DD-5-3    16:38:00    3            241           Twilight               g        B600/482    IFU             120            full   

# 2005dec04
#GN-CAL20051204-3  16:56:25    11-15        227-231       Bias         none     mirror     none  0.0            full   

# 2005dec05
#GN-2005B-DD-5-2   04:28:42    1            1             IC225        g        mirror      none            10.0           ccd2  
#GN-2005B-DD-5-2   p=27.6 q=-2.01
#GN-2005B-DD-5-2   04:33:08    2            2             IC225        g        mirror      IFU             60             full  
#GN-2005B-DD-5-2   p=0.875, q=-0.33
#GN-2005B-DD-5-2   04:39:37    3            3             IC225        g        mirror      IFU             60             full  
#GN-2005B-DD-5-2   centered on the midpoint btw the 2 nuclei
#GN-2005B-DD-5-1   04:43:50    1-2          4-5           CuAr         g        B600/478.0  IFU             120.0          full  
#GN-2005B-DD-5-1   Junk
#GN-2005B-DD-5-1   04:52:05    3            6             IC225        g        B600/478.0  IFU             3300.0         full  
#GN-2005B-DD-5-1   ok
#GN-2005B-DD-5-1   05:49:52    4            7             GCALflat     g        B600/478.0  IFU             40.0           full  
#GN-2005B-DD-5-1   ok
#GN-2005B-DD-5-1   05:52:12    5            8             CuAr         g        B600/478.0  IFU             120.0          full  
#GN-2005B-DD-5-1   ok

# 2005dec06
#GN-2005B-DD-5-4            11:51:59    1            124           Hiltner600   g        mirror      none  1              ccd2   
#GN-2005B-DD-5-4            testing seeing (unable to guide with guide star, switched to the brighter target for test)
#GN-2005B-DD-5-4            ok

#GN-2005B-DD-5-4            11:55:34    2-3          125-126       Hiltner600   g        mirror      none  10             ccd2   
#GN-2005B-DD-5-4            ok
#GN-2005B-DD-5-4-003       starting acquistion (guiding on correct star)

#GN-2005B-DD-5-4            12:05:49    4            127           Hiltner600   g        mirror      IFU   2              full   
#GN-2005B-DD-5-4            ok

#GN-2005B-DD-5-4            12:12:00    5-7          128-130       Hiltner600   g        mirror      IFU   2              full   
#GN-2005B-DD-5-4            ok
#GN-2005B-DD-5-4-007       do not see  lower fiber bundle, but have centered star

#GN-2005B-DD-5-5            12:21:23    1-2          131-132       Hiltner600   g        B600/478.0  IFU   160.0          full   
#GN-2005B-DD-5-5            junk, lost guiding
#GN-2005B-DD-5-5-002       ok, exptime 160

#GN-2005B-DD-5-5            12:36:39    3            133           Hiltner600   g        B600/478.0  IFU   300            full   
#GN-2005B-DD-5-5            ok, exptime 300

#GN-2005B-DD-5-5            12:46:42    4            134           Hiltner600   g        B600/478.0  IFU   400            full   
#GN-2005B-DD-5-5            ok, exptime 400

#GN-2005B-DD-5-5            12:56:10    5            135           GCALflat     g        B600/478.0  IFU   40.0           full   
#GN-2005B-DD-5-5            ok

#GN-2005B-DD-5-6            13:23:30    1            136           Hiltner600   g        mirror      none  10             ccd2   
#GN-2005B-DD-5-6            ok

#GN-2005B-DD-5-5            14:37:12    6            138           CuAr         g        B600/478.0  IFU   120.0          full   
#GN-2005B-DD-5-5            ok



#Aliases
#-------
set dd5=/Users/kwebb/IFU_reduction/proc
set caldir=dd5$calib/
set rawdir=/net/sbfmaps/Volumes/Science/bmiller/bdisk/bmiller/GN-2005B-DD-5/raw/

set mygmos=/Users/kwebb/iraf/scripts/gmos/
#set mygmos=/data1/kwebb/gemini/iraf/scripts/gmos/      # oringal file structure, changed for MacOSX path names
directory.nc=1
 
#Scripts
#-------
onedspec
gemini
nmisc
#gemlocal   # not required
gmos
task gscrspec=mygmos$gscrspec.cl
#task specx2w=mygmos$specx2w.cl     # not required
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
task gscombine=mygmos$gscombine.cl
task gscalibrate=mygmos$gscalibrate.cl
task chkblocks=mygmos$chkblocks.cl
task gkeywpars=mygmos$gkeywpars.cl
task gfshift=mygmos$gfshift.cl
task gfxcor=mygmos$gfxcor.cl
rv

#task $pqecorr="$foreign"   # if executing python script in qecorr from iraf, else use the execute command below
pyexecute('mygmos$pqecorr_iraf.py')

# Set common parameter values
gfreduce.rawpath="rawdir$"
gfreduce.bias="gN20051204S0227_bias.fits"
gfreduce.bpmfile="gmos$data/chipgaps.dat"
gfreduce.fl_fluxcal=no
gfreduce.fl_gscrrej=no

# use bpms for the appropriate binning, mbpm is created from the others
ifuproc.mbpm="gn_bpm2x1m.fits"
ifuproc.bpm1="mygmos$bpms/gn_bpm_ccd1_2x1f.pl"
ifuproc.bpm2="mygmos$bpms/gn_bpm_ccd2_2x1f.pl"
ifuproc.bpm3="mygmos$bpms/gn_bpm_ccd3_2x1f.pl"
ifuproc.weights="none"
ifuproc.bpmgaps="gmos$data/chipgaps.dat"

# Width of chip gap and feature width
# Defaults are for 1x1 binning, this data is 2x1
ifuproc.gap12=19 
ifuproc.gap23=19 
ifuproc.fwidth=2.



### PREPARE the bias, twilights, and flux standards once
# biases - combine
gemlist N20051204S 227-231 > bias.lis
gbias @bias.lis gN20051204S0227_bias.fits rawpath=rawdir$ fl_over-

# twilight - bias subtract
gemlist N20051203S 239-241 > twi.lis
gfreduce @twi.lis fl_extract- fl_gsappwave- fl_wavtran- fl_skysub- bias=gN20051204S0227_bias.fits

# flux standard - bias subtract
gemlist N20051206S 133,134,135,138 > dec06.lis
gfreduce @dec06.lis fl_extract- fl_gsappwave- fl_wavtran- fl_skysub- bias=gN20051204S0227_bias.fits

# files should now have prefix rg*


# IFUPROC DOES MOST OF THE WORK NOW - process the standards, Hiltner600
##### here gwen processes rgN20051206S0133 first, then rgN20051206S0134
##### gwen: twilight=rgN20051203S0240.fits fl_inter- bkgmask=s0007_blkreg.dat fl_crspec+ fl_qecor+ \
#        mbpm="gn_bpm2x1m.fits" bpm1=caldir$gn_ccd1_bpm2x1f.pl bpm2=caldir$gn_bpm_ccd2_2x1f.pl bpm3=caldir$gn_bpm_ccd3_2x1f.pl \
#        gap12=19 gap23=19 weights=none bpmgaps=gmos$data/chipgaps.dat fwidth=2. dw=.91 nw=1301 w1=4186
#### in ifuproc: twilight different, bkgmask different, flcrspec on \
#        all the same \
#        dw not set default if INDEF, nw not set default INDEF, w1 not set default INDEF
# inputs: image, flat(s), arc(s)
ifuproc rgN20051206S0134 rgN20051206S0135 rgN20051206S0138 twilight=rgN20051203S0239 fl_inter- bkgmask=s0135_blkreg.dat fl_crspec- fl_qecorr+

# files should now have prefix qtexbrg* as bad pixel masked, extracted, wavelength transformed, and then QE corrected ???

# visually confirm proper reconstructed spectra
gfdisplay sqtebrgN20051206S0134.fits ver=1 z1=0 z2=1.e8


#### here gwen measures the velocity difference between the two output slits using arcs
# gtransform ergN20051206S0138.fits wavtraname=ergN20051206S0138
# gfxcor tergN20051206S0138.fits observ="Gemini-North"
# tselect tergN20051206S0138.fits[mdf] s0138_1.fits "NO <= 750"
# tselect tergN20051206S0138.fits[mdf] s0138_2.fits "NO >= 751"
#### calculates the pixel shift
# tstat s0138_1.fits,s0138_2.fits SHIFT lowlim=-9999.0 highlim=10.
## s0138_1.fits  SHIFT
## 749      0.0155420561     0.0169317         0.012        -0.015         0.083
## s0138_2.fits  SHIFT
## 747      0.1295689425     0.0347344         0.133         0.029         0.222
##=0.1295689425-0.0155420561
##0.1140268864
#### corrects velocity difference between two output slits
# gfshift qtexbrgN20051206S0133.fits lqtexbrgN20051206S0133.fits  fl_scal=no shift1=0 shift2=-0.1140268864
# gfshift qtexbrgN20051206S0134.fits lqtexbrgN20051206S0134.fits  fl_scal=no shift1=0 shift2=-0.1140268864
#### sky subtraction
# gfskysub lqtexbrgN20051206S0133.fits expr="XINST>10" fl_inter-
# gfskysub lqtexbrgN20051206S0134.fits expr="XINST>10" fl_inter-


# SUM ALL SPECTRA
#### here gwen also does this for: slqtexbrgN20051206S0133.fits outimages="" expr="default" fl_inter-
#### then combines the two images: gscombine aslqtexbrgN20051206S0133.fits,aslqtexbrgN20051206S0134.fits \
#        Hiltner600_20051206.fits logfile="gmos.log" combine="average" scale="mean" sample="4550:4850" fl_vard-
gfapsum sqtebrgN20051206S0134.fits


# GSSTANDARD - Establish  spectrophotometric  calibration  for GMOS spectra
# correction curve should always be made interactively - avoid regions with strong absorption lines by using keystroke
#       commands to delete and add boxes away from the absorption lines
#   a: adds new bandpass, aa: make a box, d: deletes boxes closest to cusur, q: quit

#### here gwen does this for the combined statdards: gsstandard Hiltner600_20051206.fits fl_inter+ sfile="std" sfunction="sens" \
#       starnam="h600" caldir="onedstds$ctionewcal/" observa="Gemini-North" extinct="gmos$calib/mkoextinct.dat"
gsstandard asqtebrgN20051206S0134.fits starname=h600 caldir=onedstds$ctionewcal/ extinction=gmos$calib/mkoextinct.dat fl_inter-


# GALAXY science images - bias subtract
gemlist N20051205S 5-8 > dec05.lis
gfreduce @dec05.lis fl_extract- fl_gsappwave- fl_wavtran- fl_skysub- bias=gN20051204S0227_bias.fits

# IFUPROC - with one arc only, central wavelength is 478 so flat is 240 (cite gwen)
#### here gwen also does this with both arcs
ifuproc rgN20051205S0006 rgN20051205S0007 rgN20051205S0008  twilight=rgN20051203S0240  fl_inter- bkgmask=s0007_blkreg.dat fl_crspec+ fl_qecorr+ fl_skysub-

# visually inspect reconstructed image
gfdisplay qtexbrgN20051205S0006.fits ver=1 z2=1.e8


# SHIFTS BETWEEN THE SLITS - do the same as gwen above but for a different arc
#### could do for all arcs and calculate average shift as done by gwen
# We have found small pixel offsets between the two pseudo slits that are not removed by the wavelength calibration.
# The reason for this is not understood but here is one workaround using the arc to determine the shift and then
# applying it to the science data.
#gftransform ergN20051205S0008.fits wavtraname=ergN20051205S0008
#gfxcor tergN20051205S0008.fits obs=Gemini-North
#tselect tergN20051205S0008.fits[mdf] slit1.fits "NO <= 750"
#tselect tergN20051205S0008.fits[mdf] slit2.fits "NO >= 751"
#tstat slit1.fits,slit2.fits SHIFT lowlim=-1000.0 highlim=1000.
# slit1.fits  SHIFT
#  749     0.01689319103      0.019838         0.012        -0.023         0.112
# slit2.fits  SHIFT
#  747      0.1320923696     0.0360886         0.134         0.034         0.204
#=0.132-0.017
#0.115
#### calculate the shift - subtract the difference with gfshit and confirm the new shift ~0
#gfshift tergN20051205S0008.fits stergN20051205S0008.fits shift2=-0.115
#gfxcor stergN20051205S0008.fits obs=Gemini-North 
#delete slit?.fits
#tselect stergN20051205S0008.fits[mdf] slit1.fits "NO <= 750"
#tselect stergN20051205S0008.fits[mdf] slit2.fits "NO >= 751"
#tstat slit1.fits,slit2.fits SHIFT lowlim=-1000.0 highlim=1000.
# slit1.fits  SHIFT
#  749     0.01689319103      0.019838         0.012        -0.023         0.112
# slit2.fits  SHIFT
#  747     0.02548995992     0.0365223         0.028        -0.073         0.098

# CORRECT VELOCITY difference between two output slits
gfshift qtexbrgN20051205S0006.fits hqtexbrgN20051205S0006.fits shift2=-0.115

# SKY SUBtraction
gfskysub hqtexbrgN20051205S0006.fits fl_inter-

# CALIBRATE and visually inspect
gscalibrate shqtexbrgN20051205S0006.fits extinction=gmos$calib/mkoextinct.dat fl_ext+ fl_vardq+
gfdispl cshqtexbrgN20051205S0006.fits ver=1 z2=500


#### Repeat for science taken on Dec 22nd, 23rd
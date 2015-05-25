procedure ifuproc_gqecorr (image,flat,arc)

# IFU pipeline script - edited version of ifuproc.cl to include steps similar to 
# James' as well as use gqecorr rather than qecorr
# assumes that the input images are bias and overscan subtracted with gfreduce

string image     {"",prompt="Input object spectra"}
string flat      {"",prompt="List of flats"}
string arc       {"",prompt="List of arcs"}
string twilight  {"",prompt="Twilight flat"}
int    iarc      {1,prompt="Index for arc to use for flat"}
string mbpm      {"mbpm.fits",prompt="Mosaiced bpm file"}
string bpm1      {"ccd1_badpix.pl",prompt="CCD1 bpm"}
string bpm2      {"ccd2_badpix.pl",prompt="CCD2 bpm"}
string bpm3      {"ccd3_badpix.pl",prompt="CCD3 bpm"}
string bpmgaps   {"default"}
string bkgmask   {"",prompt="Mask file for gfbkgsub"}
real   xoffset   {INDEF,prompt="X offset"}
string gap12     {"default",prompt="Gap between chips 1 and 2"}
string gap23     {"default",prompt="Gap between chips 2 and 3"}
bool   fl_xshift {no,prompt="Run gspecshift?"}
bool   fl_crspec {yes,prompt="Run GSCRSPEC to clean cosmic rays?"}
bool   fl_qecorr {yes,prompt="Run QECORR for QE corrections?"}
bool   fl_skysub {yes,prompt="Subtract sky?"}
bool   fl_inter  {yes,prompt="Interactive?"}
real   fwidth    {4.,prompt="Feature width (FWHM) in pixels"}
real   gmsigma   {0.6,min=0.0,prompt="Gaussian sigma for CCD2 smoothing"}
real   gasigma   {0.0,prompt="Gaussian sigma for arc smoothing"}
string weights   {"variance",enum="variance|none",prompt="Weighting during extraction"}
real   wshift    {0.0,prompt="X wavelength shift [nm]"}
real   w1        {INDEF,prompt="Starting wavelength"}
real   w2        {INDEF,prompt="Ending wavelength"}
real   dw        {INDEF,prompt="Wavelength interval per pixel"}
int    nw        {INDEF,prompt="Numbermber of output pixels"}
int    status    {0,prompt="Exit status (0=good)"}
struct *scanfile {"",prompt="Internal use only"}
bool   fl_wrbox  {no,prompt="Mark emission lines in flats?"} # added to be similiar to mosproc, changed default to no

begin

string l_image,l_flat,l_arc,l_twilight,l_weights
string l_mbpm,l_bpm[3],l_bpmgaps,l_ccdsum,l_bkgmask,l_sky
string inst,detect,dettype,detsec,ccdsec,ccdsum,obsmode
string pre,flatlis,fflat,arclis,img,earc[2]
real l_xoffset,l_yoffset,l_apwidth
real l_fwidth,l_gasigma,l_gmsigma
real l_wshift
bool l_fl_inter,l_fl_xshift,l_recenter,l_resize,l_trace
bool l_fl_crspec,l_fl_qecorr,l_fl_skysub
real l_w1,l_w2,l_dw
int l_gap12,l_gap23,nxpix,nypix,l_torder
int xbin,ybin,l_iarc,l_nw
bool l_fl_wrbox

string tmpbpm,l_oflat,l_aflat
real dx
int ii,nflat

l_image=image ; l_flat=flat ; l_arc=arc ; l_twilight=twilight
l_mbpm=mbpm ; l_bpm[1]=bpm1 ; l_bpm[2]=bpm2; l_bpm[3]=bpm3
l_bpmgaps=bpmgaps ; l_bkgmask=bkgmask
l_xoffset=xoffset ; l_yoffset=0.0
l_fl_xshift=fl_xshift ; l_fl_inter=fl_inter
l_fwidth=fwidth ; l_gasigma=gasigma ; l_gmsigma=gmsigma
l_fl_crspec=fl_crspec ; l_fl_qecorr=fl_qecorr ; l_fl_skysub=fl_skysub
l_weights=weights ; l_iarc=iarc
l_w1=w1 ; l_w2=w2 ; l_dw=dw ; l_nw=nw
l_wshift=wshift
status=0
l_fl_wrbox=fl_wrbox

#cache("gspecshift","fparse","imgets")      # commented out as gspecshift was unrecognised

## Get various header parameters 
	hselect(l_image//"[0]","INSTRUME,OBSMODE,DETECTOR,DETTYPE","yes") | scan(inst,obsmode,detect,dettype)
	imgets(l_image//"[0]","DETECTOR")
	detect=imgets.value
	imgets(l_image//"[0]","DETTYPE")
	dettype=imgets.value
	hselect(l_image//"[SCI,1]","i_naxis1,i_naxis2","yes") | scan(nxpix,nypix)
	imgets(l_image//"[SCI,1]","CCDSUM")
	ccdsum=imgets.value
	print(ccdsum) | scan(xbin,ybin)

## Determine gap widths
	if (gap12 != "default") {
		l_gap12=int(real(gap12)+0.5)
	} else {
		l_gap12=int(38.0/real(xbin)+0.5)+2
	}
	if (gap23 != "default") {
		l_gap23=int(real(gap23)+0.5)
	} else {
		l_gap23=l_gap12
	}

## Pick chipgaps file
	if (l_bpmgaps == "default") {
		if (inst == "GMOS-S") {
			l_bpmgaps="mygmos$data/gs_chipgaps.dat"
		} else {
			l_bpmgaps="mygmos$data/gn_chipgaps.dat"
		}
	}
	
#### Not including commented out gmosaic and x offset

## Make bad pixel mask if cannot access the default 'mbpm.fits'
#### to be applied in gfreduce ??
	if (!access(l_mbpm)) {
		mkmbpm (l_mbpm, l_flat, bpm1=l_bpm[1], bpm2=l_bpm[2], bpm3=l_bpm[3],\
			   bpmgaps=l_bpmgaps, fl_wrbox=l_fl_wrbox)
 	}


# Flat, first processing
# ----------------------

	flatlis=mktemp("tmpflats")
	files(l_flat,sort-, > flatlis)
	count(flatlis) | scan(nflat)
	
	l_oflat = ""

	scanfile = flatlis
	while(fscan(scanfile, l_flat) != EOF) {
		gimverify(l_flat)
		l_flat = gimverify.outname//".fits"
		
		if (gimverify.status != 0) {
			print (l_flat//" not found or not a MEF")
			goto error
		}
		if (l_oflat=="") {
			fflat = l_flat
			if (l_bkgmask != "") {
				l_oflat = "teb"//gimverify.outname//"_flat" # i dont think this naming schemem works, should be eqb* eventually
			} else {
				l_oflat = "te"//gimverify.outname//"_flat"
			}	
			if (l_fl_qecorr) {
				l_oflat = "q"//l_oflat
			}
			if (nflat > 1) {
				l_aflat = gimverify.outname//"_avg"
			} else {
				l_aflat=gimverify.outname
			}
		}
	}
print ('>>>>> 2')	
	scanfile = ""
	if (!imaccess("eb"//l_aflat)) {		#### changed in ifuproc from 'e' to 'eb'
	
#		addbpm(l_flat, l_mbpm)
		scanfile = flatlis
		while(fscan(scanfile,l_flat) != EOF) {	
			for (ii=1; ii<=3; ii+=1) {
				nmisc.fixpix(l_flat//"[sci,"//ii//"]",l_bpm[ii], ninterp=1)
#				gemfix(l_flat//"[sci,"//ii//"]", "p"//l_flat//"[sci,"//ii//"]", method="fixpix", fl_inter-)
			}
		}

		scanfile = ""
		if (nflat > 1) {		
			gemcombine("@"//flatlis, l_aflat, combine="average", reject="avsigclip",\
					   offsets="none", scale="none", zero="none", weight="none")
		}

		gfreduce(l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, fl_over-,\
				 fl_gscrrej-, fl_fluxcal-, rawpath="", weights=l_weights,\
				 fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset, fl_inter-)#=l_fl_inter)
				 #### fl_vardq- changed to match james'

# maybe put gemfix here?				 
		## Make a background mask with gffindblocks if cannot access default mask
		if (l_bkgmask != "" && !access(l_bkgmask)) {
		
			gffindblocks(l_aflat, "e"//l_aflat, l_bkgmask)
			if (gffindblocks.status != 0) {
				goto error
			}
		}
		
		## Subtract the scattered light
		if (l_bkgmask != "" && access(l_bkgmask)) {
			
			if (!imaccess("b"//l_aflat)) {
				print("Subtracting scattered light from "//l_aflat//" using ", l_bkgmask)
				
				gfscatsub(l_aflat, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
						  yorder="5,7,5", cross=yes)
			}
# what is this for?			
#			gfreduce("b"//l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
#					 weights=l_weights, fl_extract+, fl_gsappwave+,\
#					 reference="e"//l_aflat, trace-, recenter-, wshift=l_wshift,\
#                    fl_vardq+, fl_qecorr-, xoffset=l_xoffset, rawpath="", fl_inter-)                
        }
    }
	
print ('>>>>> 3')	
# Twilight processing - optional, may not have a significant affect if included
# -------------------

	if (l_bkgmask != "") {
		pre = "b"
	} else {
		pre = ""
	}
	
	if (l_twilight != "" && !imaccess("e"//pre//l_twilight)) {
		for (ii=1; ii<=3; ii+=1) {
			nmisc.fixpix(l_twilight//"[sci,"//ii//"]",l_bpm[ii], ninterp=1)
#			gemfix(l_twilight//"[sci,"//ii//"]", l_bpm[ii], ninterp=1)
		}
		if (l_bkgmask != "" && access(l_bkgmask)) {
			if (!imaccess("b"//l_twilight)) {
				print("Subtracting background from "//l_twilight//" using ", l_bkgmask)
				gfscatsub(l_twilight, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
						 yorder="5,7,5", cross=yes)
			}
		}
		gfreduce(pre//l_twilight, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
				weights=l_weights, fl_extract+, fl_gsappwave+, fl_over-,\
				reference="e"//l_aflat, trace-, recenter-, fl_vardq+, xoffset=l_xoffset, \
				rawpath="", fl_inter-)
	}
	
print ('>>>>> 4')	
# Arc processing
# --------------
	
	arclis=mktemp("tmparclis")
	files(l_arc, sort-, > arclis)
	ii = 0
	earc[1] = ""
	earc[2] = ""
	scanfile = arclis
	while(fscan(scanfile,img) != EOF) {
		ii+=1
		if (!imaccess("e"//img)) {
			gfreduce(img, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, fl_gscrrej-,\
					 rawpath="", fl_inter=no, trace-, recenter-, ref="e"//l_aflat, fl_over-,\
					 weights=l_weights, xoffset=l_xoffset, fl_extract+, fl_gsappwave+)
            gswavelength("e"//img, fwidth=(1.7*l_fwidth), cradius=(1.2*1.7*l_fwidth),\
             			 nlost=10, fl_inter=l_fl_inter, low_reject=2.5, high_reject=2.5, verbose-)
            if (gswavelength.status != 0) {
                goto error
            }
        }
        earc[ii] = "e"//img
    }
    scanfile = ""
    delete(arclis,verify-)
print ('>>>>> 5')     
     
# Determine response function
# ---------------------------

	if (!imaccess(l_oflat)) {
		if (l_bkgmask != "") {
			pre = "b"
		} else {
			pre = ""
		}
		## Works best to run this in pyraf it seems, otherwise multiple sci ext are made?
		gftransform("e"//pre//l_aflat, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
					w2=l_w2, dw=l_dw, nw=l_nw)
		if (l_twilight != "") {
			gftransform("e"//pre//l_twilight, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
					w2=l_w2, dw=l_dw, nw=l_nw)
		}

		if (l_fl_qecorr) {
			gqecorr("b"//l_aflat, outpref="q", refimages="e"//l_arc)
			gqecorr("b"//l_twilight, outpref="q", refimages="e"//l_arc)
			pre = "qb"	
#			gqecorr(pre//l_aflat, outimage="q"//pre//l_aflat, gap12=l_gap12, gap23=l_gap23, fl_fixpix+)
#			gqecorr(pre//l_twilight, outimage="q"//pre//l_twilight, gap12=l_gap12, gap23=l_gap23, fl_fixpix+)
		}
		l_sky = ""
print (l_twilight)
print (pre)		
		if (l_twilight != "") {		
 			l_sky = pre//l_twilight
			gfreduce (l_sky, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, \
				  	  fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", fl_inter=l_fl_inter, \
				  	  weights=l_weights, fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset)
		}
		
		## Re-extract the flat with the scattered-light-subtracted image
		gfreduce ("qb"//l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, \
				  fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", fl_inter=l_fl_inter, \
				  weights=l_weights, fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset)

print ('>>>>> 6')
		## Apply fibre throughput correction			
		gfresponse("eqb"//l_aflat, l_oflat, skyimage="e"//l_sky, order=95,\
				   fl_inter=l_fl_inter, fl_fit-) #### removed iorder=750, sorder=5
		if (gfresponse.status != 0) {
			goto error
		}
	}

print ('>>>>> 7')
# Science processing
# ------------------

### James' order: gfscatsub, gemcrspec, gemfix, gqecorr, gfreduce (sky+, ref, resp, wav, vardq+, sepslits+, w1, w2, dw, nw,), uncorrupt MDF

	gimverify(l_image)
	l_image = gimverify.outname//".fits"
	
	if (l_bkgmask != "") {
		pre = "eb"
	} else {
		pre = "e"
	}
	
	if (!imaccess(pre//l_image)) {
		if (l_bkgmask != "" && access(l_bkgmask)) {	
			if (!imaccess("b"//l_image)) {
				print("Subtracting background from "//l_image//" using ", l_bkgmask)
				gfscatsub(l_image, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
							 yorder="5,7,5", cross=yes)
			}
		}
		if (!imaccess("xb"//l_image)) {
			# although not optimized for IFU, ok with faint data
			gemcrspec ("b"//l_image, "xb"//l_image, fl_vardq+)
		}
		if (!imaccess("qxb"//l_image)) {
		
# To be changed to gemfix eventually (this time fixpix)		
			nmisc.fixpix(l_image//"[sci,2]", l_bpm[2], ninterp=1)
			nmisc.fixpix(l_image//"[sci,3]", l_bpm[3], ninterp=1)
		
			gqecorr ("xb"//l_image, refimages=earc[1], fl_vardq+)
		}

		gfreduce ("qxb"//l_image, fl_sky+, fl_flux-, trace-, recenter-, fl_vardq+, fl_inter-,\
		          reference="e"//l_aflat, response=l_oflat, wavtraname=earc[1], sepslits+,\
		          w1=l_w1, w2=l_w2, dw=l_dw, nw=l_nw, rawpath="")
	
	
#	if (!imaccess(pre//l_image)) {
#		nmisc.fixpix(l_image//"[sci,2]",l_bpm[2],  ninterp=1)
#        nmisc.fixpix(l_image//"[sci,3]", l_bpm[3], ninterp=1)
#
#		if (l_bkgmask != "" && access(l_bkgmask)) {	
#			if (!imaccess("b"//l_image)) {
#				print("Subtracting background from "//l_image//" using ", l_bkgmask)
#				gfscatsub(l_image, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
#							 yorder="5,7,5", cross=yes)
#			}
#		}
#
#		gfreduce("b"//l_image, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran+, fl_skysub=l_fl_skysub, fl_fluxcal-,\
#                 rawpath="", fl_inter-, trace-, weights=l_weights, recenter-, ref="e"//l_aflat,\
#                 response=l_oflat, wavtraname=earc[1], fl_gscrrej=l_fl_crspec, fl_vardq+,\
#                 w1=l_w1, w2=l_w2, dw=l_dw, nw=l_nw, xoffset=l_xoffset)
#		
#		
#		gfreduce ("b"//l_image, fl_vardq+, fl_addmdf-, fl_trim-, fl_bias-, fl_over-,\		
#			      fl_gscrrej=l_fl_crspec, fl_extract+, fl_gsappwave+, fl_skysub=l_fl_skysub,\
#			      fl_wavtran+, fl_fluxcal-, rawpath="", weights=l_weights, reference="e"//l_aflat,\
#			      response=l_oflat, wavtraname=earc[1], trace-, recenter-, xoffset=l_xoffset,\
#			      sepslits+, w1=l_w1, w2=l_w2, dw=l_dw, nw=l_nw, fl_inter-)
                 #### removed: gap12=l_gap12, gap23=l_gap23, wavtranalt=earc[2], datares=l_fwidth, \
                 #### wshift=l_wshift, perlen=l_perlen, wavtranalt=earc[2], mbpmfile=l_mbpm
                 # chnaged addmdf- so won't call greduce, vardq- to amtch above
        if (gfreduce.status != 0) {
            goto error
        }
    }

goto clean

error:
    status=1
    goto clean

clean:
    delete(flatlis,verify-)

end
		
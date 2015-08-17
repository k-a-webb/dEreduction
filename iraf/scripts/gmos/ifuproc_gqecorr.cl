procedure ifuproc_gqecorr (image,flat,arc)

# IFU pipeline script - edited version of ifuproc.cl to include steps similar to 
# James' as well as use gqecorr rather than qecorr
# assumes that the input images are bias and overscan subtracted with gfreduce

# Currently capable of handling up to 4 flats and 2 arcs.
# This is currently designed to work with James Turner's version of gfresponse which includes the varianle
# 'wavtraname'. For sue with the iraf standard package, remove this variable declaration.

string image     {"",prompt="Input object spectra"}
string flat      {"",prompt="List of flats"}
string arc       {"",prompt="List of arcs"}
string twilight  {"",prompt="Twilight flat"}
int    iarc      {1,prompt="Index for arc to use for flat"}
string mbpm      {"mbpm.fits",prompt="Mosaiced bpm file"}
string umbpm      {"umbpm.fits",prompt="Unmosaiced bpm file"}
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
bool   fl_apsum  {yes,prompt="Sum the stellar spectra?"}
bool   fl_inter  {no,prompt="Interactive?"}
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
bool   fl_sepslits {no, prompt="Subtract sky from each slit separately?"}  # For use in final gfreduce to do same as James
bool   fl_fixgaps {no,prompt="Re-interpolate chip gaps after extraction?"}

begin

string l_image,l_flat,l_arc,l_twilight,l_weights
string l_mbpm,l_umbpm,l_bpm[3],l_bpmgaps,l_ccdsum,l_bkgmask,l_sky
string inst,detect,dettype,detsec,ccdsec,ccdsum,obsmode
string pre,flatlis,fflat,arclis,img,earc[2],l_flat_one[4]
real l_xoffset,l_yoffset,l_apwidth
real l_fwidth,l_gasigma,l_gmsigma
real l_wshift
bool l_fl_inter,l_fl_xshift,l_recenter,l_resize,l_trace
bool l_fl_crspec,l_fl_qecorr,l_fl_skysub,l_fl_apsum
real l_w1,l_w2,l_dw
int l_gap12,l_gap23,nxpix,nypix,l_torder
int xbin,ybin,l_iarc,l_nw
bool l_fl_wrbox, l_fl_sepslits, l_fl_fixgaps

string tmpbpm,l_oflat,l_aflat
real dx
int ii,nflat

l_image=image ; l_flat=flat ; l_arc=arc ; l_twilight=twilight
l_mbpm=mbpm ; l_umbpm=umbpm ; l_bpm[1]=bpm1 ; l_bpm[2]=bpm2; l_bpm[3]=bpm3
l_bpmgaps=bpmgaps ; l_bkgmask=bkgmask
l_xoffset=xoffset ; l_yoffset=0.0
l_fl_xshift=fl_xshift ; l_fl_inter=fl_inter
l_fwidth=fwidth ; l_gasigma=gasigma ; l_gmsigma=gmsigma
l_fl_crspec=fl_crspec ; l_fl_qecorr=fl_qecorr ; l_fl_skysub=fl_skysub ; l_fl_apsum=fl_apsum
l_weights=weights ; l_iarc=iarc
l_w1=w1 ; l_w2=w2 ; l_dw=dw ; l_nw=nw
l_wshift=wshift
status=0
l_fl_wrbox=fl_wrbox ; l_fl_sepslits = fl_sepslits ; l_fl_fixgaps=fl_fixgaps

cache("fparse","imgets")      # commented out as gspecshift was unrecognised

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


flatlis=mktemp("tmpflats")
files(l_flat,sort-, > flatlis)
count(flatlis) | scan(nflat)


print ('>>>>> Making bad pixel mask')
## Make bad pixel mask if cannot access the default 'mbpm.fits'
## Use the unmosaiced bad pixel mask with addbpm - l_umbpm
	if (!imaccess(l_mbpm)) {
        ii=0
        scanfile = flatlis
        while(fscan(scanfile,img) != EOF) {
		    ii+=1
            l_flat_one[ii] = img
        }
		mkmbpm (l_mbpm, l_flat_one[1], bpm1=l_bpm[1], bpm2=l_bpm[2], bpm3=l_bpm[3], umbpm=l_umbpm,\
			   bpmgaps=l_bpmgaps, fl_wrbox=l_fl_wrbox)
 	}


# Flat, first processing
# ----------------------
print ('>>>>> Flat - first processing')

	if (l_bkgmask != "") {
		pre = "bp"
	} else {
		pre = "p"
	}
	
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
			if (l_fl_qecorr) {
				l_oflat = "eq"//pre//gimverify.outname//"_flat"
			} else {
				l_oflat = "e"//pre//gimverify.outname//"_flat"
			}
		}	
		
		if (nflat > 1) {
			l_aflat = gimverify.outname//"_avg"
		} else {
			l_aflat=gimverify.outname
		}
	}
	scanfile = ""

	if (nflat > 1) {		
		gemcombine("@"//flatlis, l_aflat, combine="average", reject="avsigclip",\
				   offsets="none", scale="none", zero="none", weight="none", fl_vardq+)
	}

	if (!imaccess("p"//l_aflat)) {
		addbpm(l_aflat, l_umbpm)
		gemfix(l_aflat, "p"//l_aflat, method="fit1d", bitmask=1, order=15, fl_inter-)
	}	
	
	## Extract the spectra
	if (!imaccess("ep"//l_aflat)) { 
		gfreduce("p"//l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
				 fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", weights=l_weights,\
				 fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset, fl_inter=l_fl_inter)
	}
		
	if (!imaccess(pre//l_aflat)) {
		
		## Make background mask from extracted flat spectra					 
		if (l_bkgmask != "" && !access(l_bkgmask)) {
			gffindblocks("p"//l_aflat, "ep"//l_aflat, l_bkgmask)
			if (gffindblocks.status != 0) {
				goto error
			}
		}

		## Subtract the scattered light
		if (l_bkgmask != "" && access(l_bkgmask)) {
			print("Subtracting scattered light from "//l_aflat//" using ", l_bkgmask)
			gfscatsub("p"//l_aflat, l_bkgmask, prefix="b", xorder="5,9,5",\
					  yorder="5,7,5", cross=yes)
		}
    }
    
	
# Twilight processing - optional, may not have a significant affect if included
# -------------------	
print ('>>>>> Twilight - first processing')
	
	if (l_twilight != "" && !imaccess(pre//l_twilight)) {
	
		if (!imaccess("p"//l_twilight)) {
			addbpm(l_twilight, l_umbpm)
			gemfix(l_twilight, "p"//l_twilight, method="fit1d", bitmask=1, order=15, fl_inter-)
		}
		
		if (l_bkgmask != "" && access(l_bkgmask)) {
			print("Subtracting background from "//l_twilight//" using ", l_bkgmask)
			gfscatsub("p"//l_twilight, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
					 yorder="5,7,5", cross=yes)
		}
		if (!imaccess("e"//pre//l_twilight)) {
			gfreduce(pre//l_twilight, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
					 weights=l_weights, fl_extract+, fl_gsappwave+, fl_over-, reference="ep"//l_aflat,\
					 trace-, recenter-, fl_vardq+, xoffset=l_xoffset, rawpath="", fl_inter-)
		}
	}
	#### Note that the order is different here (p > b > then e below) than for the flats (p > e > b (like james'))
	
	
# Arc processing
# --------------
print ('>>>>> Arc processing')	

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
					 rawpath="", fl_inter-, trace-, recenter-, ref="ep"//l_aflat, fl_over-,\
					 weights=l_weights, xoffset=l_xoffset, fl_extract+, fl_gsappwave+, fl_vardq+)
		}			 
		gswavelength("e"//img, fwidth=(1.7*l_fwidth), cradius=(1.2*1.7*l_fwidth),\
					 nlost=10, fl_inter=l_fl_inter, low_reject=2.5, high_reject=2.5, verbose-)
		if (gswavelength.status != 0) {
			goto error
		}
        earc[ii] = "e"//img
    }
    scanfile = ""
    delete(arclis,verify-)
        
     
# Determine response function
# ---------------------------
print ('>>>>> Determining response function')

	if (!imaccess(l_oflat)) {


#### James here only gftransforms the arc as a self consistency check, he does nothing with the flats/twilights		
#		if (!imaccess("te"//pre//l_aflat)) {
#			if (!imaccess("e"//pre//l_aflat)) {
#				gfreduce(pre//l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
#					 fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", weights=l_weights,\
#					 fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset, fl_inter-)#=l_fl_inter)
#			}			
#			gftransform("e"//pre//l_aflat, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
#						w2=l_w2, dw=l_dw, nw=l_nw)
#		}							
#		if (l_twilight != "") {
#			if (!imaccess("te"//pre//l_twilight)) {		
#				gftransform("e"//pre//l_twilight, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
#						w2=l_w2, dw=l_dw, nw=l_nw)
#			}
#		}
	
	
		## Apply correction for differences is QE variation between detectors to flat and twilight
		if (l_fl_qecorr) {
			if (!imaccess("q"//pre//l_aflat)) {
				gqecorr(pre//l_aflat, outpref="q", refimages=earc[l_iarc])
			}
			if (l_twilight != "") {
				if (!imaccess("q"//pre//l_twilight)) {
					gqecorr(pre//l_twilight, outpref="q", refimages=earc[l_iarc])	
				}
			}
			pre = "q"//pre
		}
		
		l_sky = ""
		if (l_twilight != "") {
 			if (!imaccess("e"//pre//l_twilight)) {
				gfreduce (pre//l_twilight, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, \
						  fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", fl_inter=l_fl_inter, \
						  weights=l_weights, fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset)
			}
			l_sky = "e"//pre//l_twilight
		}
		
		## Re-extract the flat with the scattered-light-subtracted image
		if (!imaccess("e"//pre//l_aflat)) {
			gfreduce (pre//l_aflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, \
					  fl_over-, fl_gscrrej-, fl_fluxcal-, rawpath="", fl_inter=l_fl_inter, \
					  weights=l_weights, fl_extract+, fl_gsappwave+, fl_vardq+, xoffset=l_xoffset)
		}

		## Apply fibre throughput correction			
		gfresponse("e"//pre//l_aflat, l_oflat, skyimage=l_sky, order=95, fl_inter=l_fl_inter, fl_fit-, \
		    wavtraname=earc[l_iarc])
		    #### removed iorder=750, sorder=5
		if (gfresponse.status != 0) {
			goto error
		}
	}


# Science processing
# ------------------

### James' order: gfscatsub, gemcrspec, gemfix, gqecorr, gfreduce, uncorrupt MDF

print ('>>>>> Science processing')
	gimverify(l_image)
	l_image = gimverify.outname//".fits"
	
	if (l_bkgmask != "") {
		pre = "bp"
	} else {
		pre = "p"
	}
	
	if (!imaccess("steqpx"//pre//l_image)) {
		if (!imaccess("p"//l_image)) {
			addbpm(l_image, l_umbpm)
			gemfix (l_image, "p"//l_image, method="fit1d", bitmask=1, order=5, fl_inter+) #### check if order is appropriate
		}

		## Remove scattered light from the scient data
		if (l_bkgmask != "" && access(l_bkgmask)) {		
			if (!imaccess("bp"//l_image)) {
				print("Subtracting background from "//l_image//" using ", l_bkgmask)
				gfscatsub("p"//l_image, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
						  yorder="5,7,5", cross=yes)
			}
		}
		
		## Clean CR from science data. Although gemcrspec is not optimized for IFU, ok with faint data
		if (!imaccess("x"//pre//l_image)) {
			gemcrspec (pre//l_image, "x"//pre//l_image, fl_vardq+)
		}
		
		## Apply correction for differences in QE variation to science data
		if (!imaccess("px"//pre//l_image)) {
			gemfix ("x"//pre//l_image, "px"//pre//l_image, method="fixpix", bitmask=8, fl_inter-)
		}

		if (!imaccess("qpx"//pre//l_image)) {
		
			gqecorr ("px"//pre//l_image, refimages=earc[l_iarc], outpref="q", fl_vardq+)
		}

		## Flat field, wavelength calibrate & sky subtract each science dataset
		if (!imaccess("steqpx"//pre//l_image)) {
			gfreduce ("qpx"//pre//l_image, fl_sky+, fl_flux-, trace-, recenter-, fl_vardq+, fl_inter-,\
					  reference="ep"//l_aflat, response=l_oflat, wavtraname=earc[1], \
					  w1=l_w1, w2=l_w2, dw=l_dw, nw=l_nw, rawpath="", sepslits=fl_sepslits, fl_fixgaps=l_fl_fixgaps)
			if (gfreduce.status != 0) {
            	goto error
        	}
		}
    }

    if (fl_apsum && !imaccess("asteqpx"//pre//l_image)) {
        gfapsum("steqpx"//pre//l_image, fl_inter-)
    }

print ('>>>>> IFUPROC DONE <<<<<')

goto clean

error:
    status=1
    goto clean

clean:
    delete(flatlis,verify-)

end
		
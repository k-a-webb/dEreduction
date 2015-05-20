procedure ifuproc_gqecorr (image,flat,arc)

# IFU pipeline script - edited version of ifuproc.cl to include steps similar to 
# James' as well as use gqecorr rather than qecorr
# assumes that the input images are bias and overscan subtracted with gfreduce

string image   {"",prompt="Input object spectra"}
string flat    {"",prompt="List of flats"}
string arc     {"",prompt="List of arcs"}
string twilight {"",prompt="Twilight flat"}
int    iarc     {1,prompt="Index for arc to use for flat"}
string mbpm    {"mbpm.fits",prompt="Mosaiced bpm file"}
string bpm1    {"ccd1_badpix.pl",prompt="CCD1 bpm"}
string bpm2    {"ccd2_badpix.pl",prompt="CCD2 bpm"}
string bpm3    {"ccd3_badpix.pl",prompt="CCD3 bpm"}
string bpmgaps {"default"}
string bkgmask {"",prompt="Mask file for gfbkgsub"}
real   xoffset {INDEF,prompt="X offset"}
string    gap12   {"default",prompt="Gap between chips 1 and 2"}
string    gap23   {"default",prompt="Gap between chips 2 and 3"}
bool   fl_xshift {no,prompt="Run gspecshift?"}
bool   fl_crspec {yes,prompt="Run GSCRSPEC to clean cosmic rays?"}
bool   fl_qecorr {yes,prompt="Run QECORR for QE corrections?"}
bool   fl_skysub {yes,prompt="Subtract sky?"}
bool   fl_inter {yes,prompt="Interactive?"}
real   fwidth   {4.,prompt="Feature width (FWHM) in pixels"}
real   gmsigma  {0.6,min=0.0,prompt="Gaussian sigma for CCD2 smoothing"}
real   gasigma  {0.0,prompt="Gaussian sigma for arc smoothing"}
string weights  {"variance",enum="variance|none",prompt="Weighting during extraction"}
real   perlen   {100.,prompt="Percentage of spectrum to extract."}
real   wshift   {0.0,prompt="X wavelength shift [nm]"}
real   w1       {INDEF,prompt="Starting wavelength"}
real   w2       {INDEF,prompt="Ending wavelength"}
real   dw       {INDEF,prompt="Wavelength interval per pixel"}
int    nw       {INDEF,prompt="Numbermber of output pixels"}
int    status   {0,prompt="Exit status (0=good)"}
struct *scanfile {"",prompt="Internal use only"}
bool   fl_wrbox {no,prompt="Mark emission lines in flats?"} # added to be similiar to mosproc, changed default to no

begin

string l_image,l_flat,l_arc,l_twilight,l_weights
string l_mbpm,l_bpm[3],l_bpmgaps,l_ccdsum,l_bkgmask,l_sky
string inst,detect,dettype,detsec,ccdsec,ccdsum,obsmode
string pre,flatlis,fflat,arclis,img,earc[2]
real l_xoffset,l_yoffset,l_apwidth
real l_fwidth,l_gasigma,l_gmsigma
real l_wshift,l_perlen
bool l_fl_inter,l_fl_xshift,l_recenter,l_resize,l_trace
bool l_fl_crspec,l_fl_qecorr,l_fl_skysub
real l_w1,l_w2,l_dw
int l_gap12,l_gap23,nxpix,nypix,l_torder
int xbin,ybin,l_iarc,l_nw
bool l_fl_wrbox

string tmpbpm,l_oflat,l_avgflat
real dx
int ii,nflat

l_image=image ; l_flat=flat ; l_arc=arc ; l_twilight=twilight
l_mbpm=mbpm ; l_bpm[1]=bpm1 ; l_bpm[2]=bpm2 ; l_bpm[3]=bpm3
l_bpmgaps=bpmgaps ; l_bkgmask=bkgmask
l_xoffset=xoffset ; l_yoffset=0.0
l_fl_xshift=fl_xshift ; l_fl_inter=fl_inter
#l_recenter=recenter ; l_resize=resize 
#l_trace=trace
#l_apwidth=apwidth
#l_gap12=gap12 ; l_gap23=gap23
l_fwidth=fwidth ; l_gasigma=gasigma ; l_gmsigma=gmsigma
#l_torder=torder
l_fl_crspec=fl_crspec ; l_fl_qecorr=fl_qecorr ; l_fl_skysub=fl_skysub
l_weights=weights ; l_iarc=iarc
l_w1=w1; l_w2=w2; l_dw=dw; l_nw=nw
l_perlen=perlen ; l_wshift=wshift
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
	
#### Not includin gcommented out gmosaic and x offset

## Make bad pixel mask if cannot access the default 'mbpm.fits'
#### James manually creates applies a bpm via inspection, addbpm
	if (!access(l_mbpm)) {
		mkbpm (l_mmbpm, l_flat, bpm1=l_bpm[1],bpm2=l_bpm[2],bpm3=l_bpm[3],bpmgaps=l_bpmgaps,fl_wrbox=l_fl_wrbox)
	}

# Flat, first processing
# ----------------------

	flatlis=mktemp("tmpflats")
	files(l_flat,sort-, > flatlis)
	count(flatlis) | scan(nflat)
	
	l_oflat = "" #### Making change here that because I see that this is force so removing if statement below
	
	scanfile = flatlis
	while(fscan(scanfile, l_flat) != EOF) {
		gimverify(l_flat)
		l_flat = gimverify.outname//".fits"
		if (gimverify.status != 0) {
			print (l_flat//" not found or not a MEF")
			goto error
		}
		#if (l_oflat == "") {
			fflat = l_flat
			if (l_bkgmask !+ "") {
				l_oflat = "teb"//gimverify.outname//"_flat"
			} else {
				l_oflat = "te"//gimverify.outname//"_flat"
			}	
			if (l_fl_qecorr) {
				l_oflat = "q"//l_oflat
			}
			if (nflat > 1) {
				l_avgflat = gimverify.outname//"avg"
			} else {
				l_avgflat=gimverify.outname
			}
			
		#}
	}
	
	scanfile = ""
	if (!imaccess("eb"//l_avgflat)) {		#### changed in ifuproc from 'e' to 'eb'
		
		scanfile = flatlis
		while(fscan(scanfile,l_flat) != EOF) {	
			for (ii=1; ii<=3; ii+=1) {
				#### changing use of nmisc.fixpix to gemfix
				gemfix(l_flat//"[sci,"//ii//"]", l_bpm[ii], method="fixpix")
			}
		}
		scanfile = ""
		if (nflat > 1) {
			gemcombine("@"//flatlis, l_avgfalt, combine="average", reject="avsigclip",\
					   offsets="none", scale="none", zero="none", weight="none")
		}
		
		gfreduce(l_avgflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
				 fl_gscrrej-, fl_fluxcal-, rawpath="", fl_inter=l_fl_inter,weights=l_weights,\
				 gsigma=l_gmsigma, fl_extract+, fl_gsappwave+, wshift=l_wshift,\
				 perlen=l_perlen, mbpmfile=l_mbpm, fl_vardq+, fl_qecorr-, xoffset=l_xoffset)
		
				 
		## Apply background mask with gffindblocks if cannot access default mask
		if (l_bkgmask != "" && !access(l_bkgmask)) {
			gffindblocks(l_avgflat, "e"//L_avgflat, l_bkgmask)
			
			if (gffindblocks.status != 0) {
				goto error
			}
		}
		
		## Subtract the scattered light
		if (l_bkgmask != "" && access(l_bkgmask)) {
			
			if (!imaccess("b"//l_avgflat)) {
				print("Subtracting scattered light from "//l_avgflat//" using ", l_bkgmask)
				gfscatsub(l_avgfalt, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
						  yorder="5,7,5", cross=yes)
			}
			gfreduce("b"//l_avgflat, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
					 weights=l_weights, gsigma=l_gmsigma, fl_extract+, fl_gsappwave+,\
					 reference="e"//l_avgflat, trace-, recenter-, wshift=l_wshift, perlen=l_perlen,\
                     mbpmfile=l_mbpm, fl_vardq+, fl_qecorr-, xoffset=l_xoffset, rawpath="")                
        }
    }
	
	
# Twilight processing - optional, may not have a significant affect if included
# -------------------

	if (l_bkgmask != "") {
		pre = "b"
	} else {
		pre = ""
	}
	
	if (l_twilihgt != "" && !imaccess("e"//pre//l_twilight)) {
		for (ii=1; ii<=3; ii+=1) {
			gemfix(l_twilight//"[sci,"ii//"]", l_bpm[ii], ninterp=1)
		}
		if (l_bkgmask != "" && access(l_bkgmask)) {
			if (!imaccess("b"//l_twilight)) {
				print("Subtracting background from "//l_twilight//" using ", l_bkgmask)
				gfcatsub(l_twilight, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
						 yorder="5,7,5", cross=yes)
			}
		}
		gfreduce(pre//l_twilight, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-,\
				weights=l_weights, gsigma=l_gmsigma, fl_extract+, fl_gsappwave+,\
				reference="e"//l_avgflat, trace-, recenter-, wshift=l_wshift, perlen=l_perlen,\
				mbpmfile=l_mbpm, fl_vardq+, fl_qecorr-, xoffset=l_xoffset, rawpath="")
	}
	
	
# Arc processing
# --------------
	
	arclis=mktemp("tmparclis")
	file(l_arc, sort-, > arclis)
	ii = 0
	earc[1] = ""
	earc[2] == ""
	scanfile = arclis
	while(fscan(scanfile,img) != EOF) {
		ii+=1
		if (!imaccess("e"//img)) {
			gfreduce(img, fl_addmdf-, fl_trim-, fl_bias-, fl_wavtran-, fl_skysub-, fl_gscrrej-,\
					 rawpath="", fl_inter=no, trace-, recenter-, ref="e"//l_avgflat,\
					 weights=l_weights, gsigma=l_gmsigma, wshift=l_wshift, perlen=l_perlen,\
					 fl_qecorr-, xoffset=l_xoffset)
            gswavelength("e"//img, fwidth=(1.7*l_fwidth), gsigma=l_gasigma, cradius=(1.2*1.7*l_fwidth),\
             			 nlost=10, fl_inter=l_fl_inter, low_reject=2.5,high_reject=2.5)
            if (gswavelength.status != 0) {
                goto error
            }
        }
        earc[ii] = "e"//img
    }
    scanfile = ""
    delete(arclis,verify-)
     
     
# Determine response function
# ---------------------------

	if (!imaccess(l_oflat)) {
		if (l_bkgmask != "") {
			pre = "b"
		} else {
			pre = ""
		}
		gftransform("e"//pre//l_avgflat, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
					w2=l_w2, dw=l_dw, nw=l_nw)
		if (l_twilight != "") {
			gftransform("e"//pre//l_twilight, wavtraname=earc[l_iarc], fl_vardq+, w1=l_w1,\
					w2=l_w2, dw=l_dw, nw=l_nw)
		}
		pre = "te"//pre
		if (l_fl_qecorr) {
			gqecorr(pre//l_avgflat, outimage="q"//pre//l_avgflat, gap12=l_gap12,\
					gap23=l_gap23, fl_fixpix+)
			gqecorr(pre//l_twilight, outimage="q"//pre//l_twilight, gap12=l_gap12,\
					gap23=l_gap23, fl_fixpix+)
			pre = "q"//pre
		}
		l_sky = ""
		if (l_twilight != "") {
			l_sky = pre//l_twilight
		}
		gfresponse(pre//l_avgflat, l_oflat, skyimage=l_sky, order=95, iorder=750,\
				   fl_inter=l_fl_inter, sorder=5, fl_fit-)
		if (gfresponse.status != 0) {
			goto error
		}
	}


# Science processing
# ------------------

	gimverify(l_image)
	l_image = gimverify.outname//".fits"
	if (l_bkgmask != "") {
		pre = "eb"
	} else {
		pre = "e"
	}
	if (!imaccess(pre//l_image)) {
		gemfix(l_image//"[sci,2]", l_bpm[2], ninterp=1)
		gemfix(l_image//"[sci,3]", l_bpm[3], ninterp=1)
		if (l_bkgmask != "" && access(l_bkgmask)) {	
			if (!imaccess("b"//l_image)) {
				print("Subtracting background from "//l_image//" using ", l_bkgmask)
				gfcatsub(l_image, l_bkgmask, outimage="", prefix="b", xorder="5,9,5",\
							 yorder="5,7,5", cross=yes)
			}
		}
		gfreduce("b"//l_image, fl_addmdf+, fl_trim-, fl_bias-, fl_wavtran+, fl_qecorr=l_fl_qecorr,\
				 fl_skysub=l_fl_skysub, fl_fluxcal-, rawpath="", fl_inter=no, trace-,\
				 weights=l_weights, recenter-, ref="e"//l_avgflat, response=l_oflat, wavtraname=earc[1],\
                 fl_gscrrej=l_fl_crspec, datares=l_fwidth, fl_vardq+, mbpmfile=l_mbpm,\
                 gap12=l_gap12, gap23=l_gap23, gsigma=l_gmsigma, wavtranalt=earc[2], w1=l_w1,\
                 w2=l_w2, dw=l_dw, nw=l_nw, wshift=l_wshift, perlen=l_perlen, xoffset=l_xoffset)
        if (gfreduce.status != 0) {
            goto error
        }
    }

goto clean

error:
    status=1
    goto clean

clean:
vg    delete(flatlis,verify-)

end
		












































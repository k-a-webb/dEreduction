
procedure align_cubes (input, short)

string input   {"",prompt="Input file name"}
string short   {"",prompt="Input file number"}
string prefix  {"shsteqpxbrg",prompt="Input file prefix"}

begin

    string l_input, l_short, l_pre

    l_input=input  # should be of the format N20051205S0006
    l_short=short
    l_pre=prefix

    ## Decimal Shift
    ## -------------


        mkdir ("tmp_cal/"//l_short)
        mkdir ("tmp_cal/"//l_short//"/orig")
        mkdir ("tmp_cal/"//l_short//"/science")
        mkdir ("tmp_cal/"//l_short//"/var")
        mkdir ("tmp_cal/"//l_short//"/dq")
        mkdir ("tmp_cal/"//l_short//"/science/big")
        mkdir ("tmp_cal/"//l_short//"/science/new")
        mkdir ("tmp_cal/"//l_short//"/science/orig")
        mkdir ("tmp_cal/"//l_short//"/science/shift")
        mkdir ("tmp_cal/"//l_short//"/var/big")
        mkdir ("tmp_cal/"//l_short//"/var/new")
        mkdir ("tmp_cal/"//l_short//"/var/orig")
        mkdir ("tmp_cal/"//l_short//"/var/shift")
        mkdir ("tmp_cal/"//l_short//"/dq/big")
        mkdir ("tmp_cal/"//l_short//"/dq/new")
        mkdir ("tmp_cal/"//l_short//"/dq/orig")
        mkdir ("tmp_cal/"//l_short//"/dq/shift")

        ## Cut up the cube into 2D planes
        imslice ("dc"//l_pre//l_input//"[sci]", "tmp_cal/"//l_short//"/science/orig/sci_"//l_short//"_", 3, verbose-)
print ("dc"//l_pre//l_input//"[sci]")
        ## To make larger files so that the shift does not crop the image force larger image size
        ## Copy header information as well
        ## Insert the 2D images into the larger file
        # mkimage tmp_cal/6/science/new/new_sci_6_0001.fits make 0 2 "88 69"
        # mkheader tmp_cal/6/science/new/new_sci_6_0001.fits tmp_cal/6/science/orig/sci_6_001.fits append=yes
        # iminsert tmp_cal/6/science/new/new_sci_6_0001.fits tmp_cal/6/science/orig/sci_6_001.fits \
        #    tmp_cal/6/science/big/big_sci_6_0001.fits replace mkim_coords.dat

        for (i=1; i<=1301; i+=1) {
            mkimage(("tmp_cal/"//l_short//"/science/new/new_sci_"//l_short//"_0000"+i)//".fits", "make", 0, 2, "88 69")
            mkheader(("tmp_cal/"//l_short//"/science/new/new_sci_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/science/orig/sci_"//l_short//"_000"+i)//".fits")
            iminsert(("tmp_cal/"//l_short//"/science/new/new_sci_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/science/orig/sci_"//l_short//"_000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/science/big/big_sci_"//l_short//"_0000"+i)//".fits","replace","mkim_coords.dat")
        }
print ('>>>>> 1')
        ## Shift the 2D images
        for (i=1;i<=1301;i+=1) {
            imshift(("tmp_cal/"//l_short//"/science/big/big_sci_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"/science/shift/shift_sci_"//l_short//"_0000"+i)//".fits",-9.5,0.75)
        }
print ('>>>>> 2')
        ## Restack the 2D images into a cube
        for (i=1;i<=1301;i+=1) {
            print(("tmp_cal/"//l_short//"/science/shift/shift_sci_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_sci_"//l_short//".lis")
        }
        imstack ("@tmp_cal/2D_sci_"//l_short//".lis", "tmp_cal/shifted_sci_"//l_short//".fits")
print ('>>>>> 3')

        ## Repeat for the variance

        imslice ("dc"//l_pre//l_input//"[var]", "tmp_cal/"//l_short//"/var/orig/var_"//l_short//"_", 3, verbose-)
        for (i=1; i<=1301; i+=1) {
            mkimage(("tmp_cal/"//l_short//"/var/new/new_var_"//l_short//"_0000"+i)//".fits","make",0,2,"88 69")
            mkheader(("tmp_cal/"//l_short//"/var/new/new_var_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/var/orig/var_"//l_short//"_000"+i)//".fits")
            iminsert(("tmp_cal/"//l_short//"/var/new/new_var_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/var/orig/var_"//l_short//"_000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/var/big/big_var_"//l_short//"_0000"+i)//".fits","replace","mkim_coords.dat")
            imshift(("tmp_cal/"//l_short//"/var/big/big_var_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"/var/shift/shift_var_"//l_short//"_0000"+i)//".fits",-9.5,0.75)
            print(("tmp_cal/"//l_short//"/var/shift/shift_var_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_var_"//l_short//".lis")
        }
        imstack ("@tmp_cal/2D_var_"//l_short//".lis", "tmp_cal/shifted_var_"//l_short//".fits")
print ('>>>>> 4')

        ## And again for the DQ

        imslice ("dc"//l_pre//l_input//"[dq]", "tmp_cal/"//l_short//"/dq/orig/dq_"//l_short//"_", 3, verbose-)
        for (i=1; i<=1301; i+=1) {
            mkimage(("tmp_cal/"//l_short//"/dq/new/new_dq_"//l_short//"_0000"+i)//".fits","make",0,2,"88 69")
            mkheader(("tmp_cal/"//l_short//"/dq/new/new_dq_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/dq/orig/dq_"//l_short//"_000"+i)//".fits")
            iminsert(("tmp_cal/"//l_short//"/dq/new/new_dq_"//l_short//"_0000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/dq/orig/dq_"//l_short//"_000"+i)//".fits",\
                     ("tmp_cal/"//l_short//"/dq/big/big_dq_"//l_short//"_0000"+i)//".fits","replace","mkim_coords.dat")
            imshift(("tmp_cal/"//l_short//"/dq/big/big_dq_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"/dq/shift/shift_dq_"//l_short//"_0000"+i)//".fits",-9.5,0.75)
            print(("tmp_cal/"//l_short//"/dq/shift/shift_dq_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_dq_"//l_short//".lis")
        }
        imstack ("@tmp_cal/2D_dq_"//l_short//".lis", "tmp_cal/shifted_dq_"//l_short//".fits")
print ('>>>>> 5')

        ## Change the DQ plane to integer type rather than real
        imreplace ("tmp_cal/shifted_dq_"//l_short//".fits", 1., lower=0.01, upper=INDEF)
        chpixtype ("tmp_cal/shifted_dq_"//l_short//".fits", "tmp_cal/shifted_dq_"//l_short//".fits", "short", verb-)
print ('>>>>> 6')

        ## Merge all the planes into an aligned MEF

        imcopy ("dc"//l_pre//l_input//".fits[0]", "adc"//l_pre//l_input//"_decimal.fits")
        tcopy ("c"//l_pre//l_input//".fits[1]", "adc"//l_pre//l_input//"_decimal.fits[1]")
        imcopy ("tmp_cal/shifted_sci_"//l_short//".fits", "adc"//l_pre//l_input//"_decimal[SCI,1,append+]")
        imcopy ("tmp_cal/shifted_var_"//l_short//".fits", "adc"//l_pre//l_input//"_decimal[VAR,1,append+]")
        imcopy ("tmp_cal/shifted_dq_"//l_short//".fits", "adc"//l_pre//l_input//"_decimal[DQ,1,append+]")
        hedit ("adc"//l_pre//l_input//"_decimal[0]", "nextend", 4, update+, verify-)
        hedit ("adc"//l_pre//l_input//"_decimal[0]", "nsciext", 1, update+, verify-)
print ('>>>>> 7')

    ## Integer shift
    ## -------------

        mkdir ("tmp_cal/"//l_short//"_integer")
        mkdir ("tmp_cal/"//l_short//"_integer/science")
        mkdir ("tmp_cal/"//l_short//"_integer/var")
        mkdir ("tmp_cal/"//l_short//"_integer/dq")
        mkdir ("tmp_cal/"//l_short//"_integer/science/shift")
        mkdir ("tmp_cal/"//l_short//"_integer/var/shift")
        mkdir ("tmp_cal/"//l_short//"_integer/dq/shift")

        for (i=1; i<=1301; i+=1) {
            imshift(("tmp_cal/"//l_short//"/science/big/big_sci_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"_integer/science/shift/shift_sci_"//l_short//"_0000"+i)//".fits",-9.5,0.75)
            print(("tmp_cal/"//l_short//"_integer/science/shift/shift_sci_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_sci_"//l_short//"_integer.lis")
        }
        imstack ("@tmp_cal/2D_sci_"//l_short//"_integer.lis", "tmp_cal/shifted_sci_"//l_short//"_integer.fits")
print ('>>>>> 8')
        for (i=1; i<=1301; i+=1) {
            imshift(("tmp_cal/"//l_short//"/var/big/big_var_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"_integer/var/shift/shift_var_"//l_short//"_0000"+i)//".fits",-10,1)
            print(("tmp_cal/"//l_short//"_integer/var/shift/shift_var_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_var_"//l_short//"_integer.lis")
        }
        imstack ("@tmp_cal/2D_var_"//l_short//"_integer.lis", "tmp_cal/shifted_var_"//l_short//"_integer.fits")
print ('>>>>> 9')
        for (i=1; i<=1301; i+=1) {
            imshift(("tmp_cal/"//l_short//"/dq/big/big_dq_"//l_short//"_0000"+i)//".fits",\
                    ("tmp_cal/"//l_short//"_integer/dq/shift/shift_dq_"//l_short//"_0000"+i)//".fits",-10,1)
            print(("tmp_cal/"//l_short//"/dq/shift/shift_dq_"//l_short//"_0000"+i)//".fits",>> "tmp_cal/2D_dq_"//l_short//"_integer.lis")
        }

        imstack ("@tmp_cal/2D_dq_"//l_short//"_integer.lis", "tmp_cal/shifted_dq_"//l_short//"_integer.fits")
print ('>>>>> 10')
        imreplace ("tmp_cal/shifted_dq_"//l_short//"_integer.fits", 1., lower=0.01, upper=INDEF)
        chpixtype ("tmp_cal/shifted_dq_"//l_short//"_integer.fits", "tmp_cal/shifted_dq_"//l_short//"_integer.fits", "short", verb-)

        imcopy ("dc"//l_pre//l_input//".fits[0]", "adc"//l_pre//l_input//"_integer.fits")
        tcopy ("c"//l_pre//l_input//".fits[1]", "adc"//l_pre//l_input//"_integer.fits[1]")
        imcopy ("tmp_cal/shifted_sci_"//l_short//"_integer.fits", "adc"//l_pre//l_input//"_integer[SCI,1,append+]")
        imcopy ("tmp_cal/shifted_var_"//l_short//"_integer.fits", "adc"//l_pre//l_input//"_integer[VAR,1,append+]")
        imcopy ("tmp_cal/shifted_dq_"//l_short//"_integer.fits", "adc"//l_pre//l_input//"_integer[DQ,1,append+]")
        hedit ("adc"//l_pre//l_input//"_integer[0]", "nextend", 4, update+, verify-)
        hedit ("adc"//l_pre//l_input//"_integer[0]", "nsciext", 1, update+, verify-)
print ('>>>>> 11')
end
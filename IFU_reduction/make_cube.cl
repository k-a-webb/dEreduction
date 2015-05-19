procedure make_cube (input, short)

string input   {"",prompt="Input file name"}
int short      {"",prompt="Input file number"}

begin

    string l_input
    l_input=input
    int l_short
    l_short=short
    string l_file

    l_file="shqtexbrg"//l_input

    gfcube (l_file)

    # create dunny 'var' plane
    imcopy ("c"//l_file//"[0]", "tmp_var0"//l_short//".fits")
    tcopy (l_file//"[1]", "tmp_var0"//l_short//".fits[0]")
    imcopy ("c"//l_file//"[var,1]", "tmp_var0"//l_short//".fits[SCI,1,append+]")
    hedit ("tmp_var0"//l_short//".fits[0]", "nextend", 2, update+, verify-)
    hedit ("tmp_var0"//l_short//".fits[0]", "nsciext", 1, update+, verify-)
    gfcube ("tmp_var0"//l_short//".fits")

    # create a dummy 'dq' plane
    imcopy ("c"//l_file//"[0]", "tmp_dq0"//l_short//".fits")
    tcopy (l_file//"[1]", "tmp_dq0"//l_short//".fits[0]")
    imcopy ("c"//l_file//"[dq,1]", "tmp_dq0"//l_short//".fits[SCI,1,append+]")
    hedit ("tmp_dq0"//l_short//".fits[0]", "nextend", 2, update+, verify-)
    hedit ("tmp_dq0"//l_short//".fits[0]", "nsciext", 1, update+, verify-)
    gfcube ("tmp_dq0"//l_short//".fits")

    # now stack the outputs from gfcube into a single MDF file
    imcopy ("d"//l_file//"[0]", "mdc"//l_file)
    tcopy ("c"//l_file//"[1]", "mdc"//l_file//"[1]")
    imcopy ("d"//l_file//"[sci,1]", "mdc"//l_file//"[SCI,1,append+]")
    imcopy ("dtmp_var0"//l_short//"[sci,1]", "mdc"//l_file//"[VAR,1,append+]")
    imcopy ("dtmp_dq0"//l_short//"[sci,1]", "mdc"//l_file//"[DQ,1,append+]")
    hedit ("mdc"//l_file//"[0]", "nextend", 4, update+, verify-)
    hedit ("mdc"//l_file//"[0]", "nextend", 1, update+, verify-)


end
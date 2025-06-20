# Subject:  Demonstrate the structure of the .PTBin binary files written by MIKE SHE.
#           This script does not actually do anything! It has not output whatsoever!
#           It could be extended by doing something usefull with the data read.
# Requires: numpy
# Usage:    Set hard coded file name and start
# Guarantee: none!
#            This script was deloped for MIKE SHE .PTBin file-version 5. It won't work with any other version!
# License:  DHI
# Author:   dhi\uha
# date:     06/2021
# note:     PTBin files are little endian - if you try to run this script on a system where this is not
#           the default you have to explicitly set it (e.g. use type '<i4' instead of just 'i4')!

import numpy as np
import os
off = 0 # keeping track of how many bytes have been read to assert it's all correct - and to stop reading from the file at the right moment
fname = r"C:\Work\Main\Products\Source\MSHE\RegTest\WQ\Karup\Karup_PT.she - Result Files\Karup_PT.PTBin"

with open(fname, 'r') as f:
    dt = np.dtype('i4, f8, i4')
    # MOpnBin 148
    (itype, time, dum) = np.fromfile(f, dtype=dt, count=1)[0]
    off += dt.itemsize # 16
    if(itype != -999):
        print("unexpected itype: {0}".format(itype))
        exit(-1)

    (file_type, iverno, le_file_text, le_user_text) = np.fromfile(f, dtype=np.int32, count=4)    
    off += 4 * 4 # 32
    if(file_type != 4004):
        print('Cannot read file type "{0}" - only PTBin-files (id 4004) are supported!'.format(file_type))
        exit(-1)

    (itype, time, dum) = np.fromfile(f, dtype=dt, count=1)[0]
    off += dt.itemsize # 48
    if(itype != -998):
        print("unexpected itype2: {0}".format(itype))
        exit(-1)

    # MOpnBin 186
    file_text = np.fromfile(f, dtype=np.byte, count=le_file_text)
    off += le_file_text

    if(le_user_text > 0):
        (itype, time, dum) = np.fromfile(f, dtype=dt, count=1)[0]
        off += dt.itemsize
        if(itype != -997):
            print("unexpected itype3: {0}".format(itype))
            exit(-1)        
        user_text = np.fromfile(f, dtype=np.byte, count=le_user_text)
        off += le_user_text

    # END MOpnBin

    (itype, time, dum) = np.fromfile(f, dtype=dt, count=1)[0]
    off += dt.itemsize
    pt_bin_ver = np.fromfile(f, dtype=np.int32, count=1)
    off += 4
    if(pt_bin_ver != 5): # If you find 6 or above, there is no point in continuing: the file format has changed. Reach out to DHI if necessary.
        print("Unsupported pt_bin_ver: {0}! Only supported version is 5!".format(pt_bin_ver))
        exit(-1)

    dt = np.dtype('i4,i4,5i4,5i4')
    # length of projection string, flag if user requests NO conversion of units in output, start-/end date of simulation:
    (le_map_proj, lUseHardCodedResultUnits, (start_y, start_m, start_d, start_h, start_min), (end_y, end_m, end_d, end_h, end_min)) = np.fromfile(f, dtype=dt, count=1)[0]
    off += dt.itemsize
    # print("le_map_proj: {0}".format(le_map_proj)) # 1024 expected
    
    # projection string
    map_proj = np.fromfile(f, dtype=np.byte, count=le_map_proj)
    off += le_map_proj

    ################ Actual data starts here!
    
    dt_outer = np.dtype('i4, f8')
    # 0: number of particles in this time step
    # 1: simulation time in s from start

    dt_inner = np.dtype('i4, f8, f8, f4, i4, i1, f4, f4')
    # 0: particle id
    # 1: x position
    # 2: y position
    # 3: z position
    # 4: layer index
    # 5: not used
    # 6: particle birth time in seconds from start
    # 7: particle path length

    file_size = os.path.getsize(fname)
    if(file_size <= off):
        print("Calculated offset ist at/beyond end of file before any actual data has been read - quitting!")
        exit(-1)
    sample_print = True
    # repeat until end of file (the number of time steps is not stored inside the file!):
    while(off < file_size):
        (n_part, timesec) = np.fromfile(f, dtype=dt_outer, count=1)[0]
        off += dt_outer.itemsize
        # reading all particle data for this time step at once. If that's still too much you could loop and read one particle at a time
        inners = np.fromfile(f, dtype=dt_inner, count=n_part)
        off += n_part * dt_inner.itemsize
        if sample_print:
            print(inners[-1])
            sample_print = False
    
    if(file_size != off):
        print("Calculated offset not equal to file size after reading all data - verify reading correct types/sizes and calculating offset correctly!")
        exit(-1)

print("Success!")
exit(0)

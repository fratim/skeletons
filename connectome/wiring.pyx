import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



from skeletons.utilities import dataIO



cdef extern from 'cpp-wiring.h':
    void CppUpdateResolution(float resolution[3])
    void CppUpdateBlocksize(long blocksize_inp[3])
    void CppUpdateVolumesize(long volumesize[3])
    void CppUpdateBlockindices(long block_z, long block_y, long block_x)
    void CppSkeletonGeneration(const char *prefix, const char *lookup_table_directory, const char *synapses_directory, long *inp_labels)

# extract the wiring diagram for this prefix and segment_ID
def GenerateSkeleton(prefix, block_z, block_y, block_x):
    # start running time statistics
    start_time = time.time()

    #set directory of synapses
    pointclouds_directory = dataIO.PointcloudsDirectory(prefix)


    # everything needs to be long ints to work with c++
    fileName = dataIO.SegmentationsDirectory(prefix)+prefix+"/Zebrafinch-"+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    data = dataIO.ReadH5File(fileName)
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_inp_labels =  np.ascontiguousarray(data, dtype=ctypes.c_int64)

    # get ans set Resolution from meta file
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)
    CppUpdateResolution(&(cpp_resolution[0]))

    # get, check and set grid size from meta file
    blocksize_inp = dataIO.Blocksize(prefix)
    if data.shape[0]!=blocksize_inp[0] or data.shape[1]!=blocksize_inp[1] or data.shape[2]!=blocksize_inp[2]:
        raise ValueError("Data chunk size ("+str(data.shape)+") not equal to size specified in meta file"+str(blocksize_inp)+"!")
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_blocksize_inp = np.ascontiguousarray(blocksize_inp, dtype=ctypes.c_int64)
    CppUpdateBlocksize(&(cpp_blocksize_inp[0]))

    # set volumesize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_volumesize = np.ascontiguousarray(dataIO.Volumesize(prefix), dtype=ctypes.c_int64)
    # change to float to long
    CppUpdateVolumesize(&(cpp_volumesize[0]))

    # set block indices of the current block
    CppUpdateBlockindices(block_z, block_y, block_x)

    # the look up table is in the synapseaware/connectome folder
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppSkeletonGeneration(prefix.encode('utf-8'), lut_directory.encode('utf-8'), pointclouds_directory.encode('utf-8'), &(cpp_inp_labels[0,0,0]))

    # print out statistics for wiring extraction
    print ('Generated skeletons in {:0.2f} seconds'.format(time.time() - start_time))

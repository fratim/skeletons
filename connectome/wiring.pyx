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
    void CppSkeletonGeneration(const char *prefix, const char *lookup_table_directory, long *inp_labels)
    void CppUpdateDirectories(const char* synapses_dir, const char* somae_dir, const char* skeleton_dir)

# extract the wiring diagram for this prefix and segment_ID
def GenerateSkeleton(prefix, block_z, block_y, block_x):
    # start running time statistics
    start_time = time.time()

    # everything needs to be long ints to work with c++
    fileName = dataIO.InputlabelsDirectory(prefix)+"/"+prefix+"/Zebrafinch-"+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
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

    # set block indices of the current block
    CppUpdateDirectories(dataIO.SynapsesDirectory(prefix).encode('utf-8'),dataIO.SomaeDirectory(prefix).encode('utf-8'),dataIO.SkeletonDirectory(prefix).encode('utf-8'))

    # the look up table is in the synapseaware/connectome folder
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppSkeletonGeneration(prefix.encode('utf-8'), lut_directory.encode('utf-8'), &(cpp_inp_labels[0,0,0]))

    # print out statistics for wiring extraction
    print ('Generated skeletons in {:0.2f} seconds'.format(time.time() - start_time))

import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



from skeletons.utilities import dataIO



cdef extern from 'cpp-wiring.h':
    void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3],
            const char* synapses_dir, const char* somae_dir, const char* output_dir, const char* output_dir_z_inp, const char* output_dir_y_inp, const char* output_dir_x_inp);
    void CppSkeletonRefinement(const char *prefix, long segment_ID_query, long block_ind_inp[3]);

# extract the wiring diagram for this prefix and segment_ID
def GenerateSkeleton(prefix, output_folder, block_z, block_y, block_x):
    # start running time statistics
    start_time = time.time()

    # everything needs to be long ints to work with c++
    fileName = dataIO.InputlabelsDirectory(prefix)+"/"+prefix+"/Zebrafinch-input_labels-"+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    data = dataIO.ReadH5File(fileName)
    cdef np.ndarray[long, ndim=3, mode='c'] cpp_inp_labels =  np.ascontiguousarray(data, dtype=ctypes.c_int64)

    # get Resolution from meta file
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)

    # get and check grid size from meta file
    blocksize_inp = dataIO.Blocksize(prefix)
    if data.shape[0]!=blocksize_inp[0] or data.shape[1]!=blocksize_inp[1] or data.shape[2]!=blocksize_inp[2]:
        raise ValueError("Data chunk size ("+str(data.shape)+") not equal to size specified in meta file"+str(blocksize_inp)+"!")
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_blocksize_inp = np.ascontiguousarray(blocksize_inp, dtype=ctypes.c_int64)

    # get volumesize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_volumesize = np.ascontiguousarray(dataIO.Volumesize(prefix), dtype=ctypes.c_int64)

    # get block indices
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind = np.ascontiguousarray(np.array([block_z, block_y, block_x]), dtype=ctypes.c_int64)

    # the look up table is in the synapseaware/connectome folder
    lut_directory = os.path.dirname(__file__)

    output_folder_z_minus = dataIO.OutputDirectory(prefix)+'output-'+str(block_z-1).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'
    output_folder_y_minus = dataIO.OutputDirectory(prefix)+'output-'+str(block_z).zfill(4)+'z-'+str(block_y-1).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'
    output_folder_x_minus = dataIO.OutputDirectory(prefix)+'output-'+str(block_z).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x-1).zfill(4)+'x'+'/'

    # create c++ datablock and set all variables
    CPPcreateDataBlock(     prefix.encode('utf-8'), lut_directory.encode('utf-8'), &(cpp_inp_labels[0,0,0]),
                            &(cpp_resolution[0]), &(cpp_blocksize_inp[0]), &(cpp_volumesize[0]), &(cpp_block_ind[0]),
                            dataIO.SynapsesDirectory(prefix).encode('utf-8'),dataIO.SomaeDirectory(prefix).encode('utf-8'),
                            output_folder.encode('utf-8'),output_folder_z_minus.encode('utf-8'),output_folder_y_minus.encode('utf-8'),output_folder_x_minus.encode('utf-8'))


    # # call the topological skeleton algorithm
    # CppSkeletonGeneration(prefix.encode('utf-8'), lut_directory.encode('utf-8'), &(cpp_inp_labels[0,0,0]))



    # print out statistics for wiring extraction
    print ('Generated skeletons in {:0.2f} seconds'.format(time.time() - start_time))

def RefineSkeleton(prefix, segment_ID_query, block_z, block_y, block_x):

    # get block indices
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind = np.ascontiguousarray(np.array([block_z, block_y, block_x]), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_query_ID = np.ascontiguousarray(np.array([segment_ID_query]), dtype=ctypes.c_int64)
    CppSkeletonRefinement(prefix.encode('utf-8'), cpp_query_ID, &(cpp_block_ind[0]))

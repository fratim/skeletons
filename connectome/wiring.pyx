import os
import time


cimport cython
cimport numpy as np
import ctypes
import numpy as np



from skeletons.utilities import dataIO



cdef extern from 'cpp-wiring.h':
    void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3],
            long inp_blocksize[3], long volume_size[3], long block_ind_inp[3], long block_ind_start_inp[3], long block_ind_end_inp[3],
            const char* synapses_dir, const char* somae_dir, const char* output_dir);
    void CppSkeletonRefinement(const char *prefix, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_begin[3],
            long block_ind_end[3], const char* output_dir);
    void ComputeAnchorPoints(const char *prefix, const char* output_dir, long inp_blocksize_inp[3], long blockind_inp[3], long block_ind_start_inp[3], long block_ind_end_inp[3],
            long volumesize_in_inp[3], long *z_min_wall, long *z_max_wall, long *y_min_wall, long *y_max_wall, long *x_min_wall, long *x_max_wall);
# save walls
def SaveWalls(prefix, output_folder, block_z, block_y, block_x):
    fileName = dataIO.InputlabelsDirectory(prefix)+"/"+prefix+"/Zebrafinch-input_labels-"+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    data = dataIO.ReadH5File(fileName)
    walls_folder = output_folder+'walls/'+prefix+'/'
    dataset = 'main'
    zmin_fp = walls_folder+'zMinWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    zmax_fp = walls_folder+'zMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    ymin_fp = walls_folder+'yMinWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    ymax_fp = walls_folder+'yMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    xmin_fp = walls_folder+'xMinWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    xmax_fp = walls_folder+'xMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    dataIO.WriteH5File(data[0,:,:], zmin_fp , dataset)
    dataIO.WriteH5File(data[-1,:,:],zmax_fp, dataset)
    dataIO.WriteH5File(data[:,0,:],ymin_fp, dataset)
    dataIO.WriteH5File(data[:,-1,:],ymax_fp, dataset)
    dataIO.WriteH5File(data[:,:,0],xmin_fp, dataset)
    dataIO.WriteH5File(data[:,:,-1], xmax_fp, dataset)

# save walls
def MakeAnchorpoints(prefix, output_folder, block_z, block_y, block_x):
    # get indices of the last block
    end_block_ind = np.array(dataIO.StartBlocks(prefix))+np.array(dataIO.NBlocks(prefix))-np.array((1,1,1))

    walls_folder_current = output_folder+'output-'+str(block_z).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'+'walls/'+prefix+'/'
    zmax_fp = walls_folder_current+'zMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    ymax_fp = walls_folder_current+'yMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    xmax_fp = walls_folder_current+'xMaxWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"

    walls_folder_z = output_folder+'output-'+str(block_z+1).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'+'walls/'+prefix+'/'
    walls_folder_y = output_folder+'output-'+str(block_z).zfill(4)+'z-'+str(block_y+1).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'+'walls/'+prefix+'/'
    walls_folder_x = output_folder+'output-'+str(block_z).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x+1).zfill(4)+'x'+'/'+'walls/'+prefix+'/'

    zmin_fp = walls_folder_z+'zMinWall'+str(block_z+1).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    ymin_fp = walls_folder_y+'yMinWall'+str(block_z).zfill(4)+"z-"+str(block_y+1).zfill(4)+"y-"+str(block_x).zfill(4)+"x"+".h5"
    xmin_fp = walls_folder_x+'xMinWall'+str(block_z).zfill(4)+"z-"+str(block_y).zfill(4)+"y-"+str(block_x+1).zfill(4)+"x"+".h5"

    cdef np.ndarray[long, ndim=2, mode='c'] cpp_zmax_wall = np.ascontiguousarray(dataIO.ReadH5File(zmax_fp), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=2, mode='c'] cpp_ymax_wall = np.ascontiguousarray(dataIO.ReadH5File(ymax_fp), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=2, mode='c'] cpp_xmax_wall = np.ascontiguousarray(dataIO.ReadH5File(xmax_fp), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=2, mode='c'] cpp_zmin_wall = np.ascontiguousarray(np.zeros((1,1)), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=2, mode='c'] cpp_ymin_wall = np.ascontiguousarray(np.zeros((1,1)), dtype=ctypes.c_int64)
    cdef np.ndarray[long, ndim=2, mode='c'] cpp_xmin_wall = np.ascontiguousarray(np.zeros((1,1)), dtype=ctypes.c_int64)

    if block_z < end_block_ind[0]:
        cpp_zmin_wall = np.ascontiguousarray(dataIO.ReadH5File(zmin_fp), dtype=ctypes.c_int64)

    if block_y < end_block_ind[1]:
        cpp_ymin_wall = np.ascontiguousarray(dataIO.ReadH5File(ymin_fp), dtype=ctypes.c_int64)

    if block_x < end_block_ind[2]:
        cpp_xmin_wall = np.ascontiguousarray(dataIO.ReadH5File(xmin_fp), dtype=ctypes.c_int64)

    # get blocksize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_blocksize = np.ascontiguousarray(dataIO.Blocksize(prefix), dtype=ctypes.c_int64)

    # get block indices (start)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind = np.ascontiguousarray(np.array([block_z, block_y, block_x]), dtype=ctypes.c_int64)

    # get block indices of start blocks
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_start = np.ascontiguousarray(np.array(dataIO.StartBlocks(prefix)), dtype=ctypes.c_int64)

    # end blocks
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_end = np.ascontiguousarray(np.array(end_block_ind), dtype=ctypes.c_int64)

    # get volumesize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_volumesize = np.ascontiguousarray(dataIO.Volumesize(prefix), dtype=ctypes.c_int64)

    ComputeAnchorPoints(prefix.encode('utf-8'), output_folder.encode('utf-8'), &(cpp_blocksize[0]), &(cpp_block_ind[0]),
                &(cpp_block_ind_start[0]), &(cpp_block_ind_end[0]), &(cpp_volumesize[0]), &(cpp_zmin_wall[0,0]), &(cpp_zmax_wall[0,0]), &(cpp_ymin_wall[0,0]), &(cpp_ymax_wall[0,0]), &(cpp_xmin_wall[0,0]), &(cpp_xmax_wall[0,0]));

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

    # get block indices of start blocks
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_start = np.ascontiguousarray(np.array(dataIO.StartBlocks(prefix)), dtype=ctypes.c_int64)

    # end blocks
    end_block_ind = np.array(dataIO.StartBlocks(prefix))+np.array(dataIO.NBlocks(prefix))-np.array((1,1,1))
    print("block ind end:" + str(end_block_ind))
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_end = np.ascontiguousarray(np.array(end_block_ind), dtype=ctypes.c_int64)

    # the look up table is in the synapseaware/connectome folder
    lut_directory = os.path.dirname(__file__)

    # create c++ datablock and set all variables
    CPPcreateDataBlock(     prefix.encode('utf-8'), lut_directory.encode('utf-8'), &(cpp_inp_labels[0,0,0]), &(cpp_resolution[0]), &(cpp_blocksize_inp[0]),
                            &(cpp_volumesize[0]), &(cpp_block_ind[0]), &(cpp_block_ind_start[0]), &(cpp_block_ind_end[0]),
                            dataIO.SynapsesDirectory(prefix).encode('utf-8'),dataIO.SomaeDirectory(prefix).encode('utf-8'),
                            output_folder.encode('utf-8'))

    # print out statistics for wiring extraction
    print ('Generated skeletons in {:0.2f} seconds'.format(time.time() - start_time))

# run refinement on skeleton
def RefineSkeleton(prefix, output_folder, block_z_start, block_y_start, block_x_start, block_z_end, block_y_end, block_x_end):
    # get blocksize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_blocksize = np.ascontiguousarray(dataIO.Blocksize(prefix), dtype=ctypes.c_int64)
    # get volumesize
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_volumesize = np.ascontiguousarray(dataIO.Volumesize(prefix), dtype=ctypes.c_int64)
    # get Resolution from meta file
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_resolution = np.ascontiguousarray(dataIO.Resolution(prefix)).astype(np.float32)
    # get block indices (start)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_begin = np.ascontiguousarray(np.array([block_z_start, block_y_start, block_x_start]), dtype=ctypes.c_int64)
    # get block indices (end)
    cdef np.ndarray[long, ndim=1, mode='c'] cpp_block_ind_end = np.ascontiguousarray(np.array([block_z_end, block_y_end, block_x_end]), dtype=ctypes.c_int64)

    CppSkeletonRefinement(prefix.encode('utf-8'), &(cpp_resolution[0]), &(cpp_blocksize[0]), &(cpp_volumesize[0]), &(cpp_block_ind_begin[0]), &(cpp_block_ind_end[0]), output_folder.encode('utf-8'))

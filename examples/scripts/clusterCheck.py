import numpy
import time
import numpy as np

def WriteH5File(data, filename, dataset):
    with h5py.File(filename, 'w') as hf:
        # should cover all cases of affinities/images
        hf.create_dataset(dataset, data=data, compression='gzip')

np.random.seed(11)

for _ in range(10):

    start_time = time.time()

    start_time_alloc = time.time()
    random_labels = np.random.randint(0,50000,size=(256,256,256),dtype=np.uint64)
    time_alloc = time.time()-start_time_alloc

    start_time_comp = time.time()
    for iz in range(0,random_labels.shape[0]-2):
        for iy in range(0,random_labels.shape[1]-2):
            for ix in range(0,random_labels.shape[2]-2):

                random_labels[iz,iy,ix] = random_labels[iz,iy,ix]*random_labels[iz+1,iy+1,ix+1]+random_labels[iz+2,iy+2,ix+2]
    time_comp = time.time()-start_time_comp

    start_time_write = time.time()
    WriteH5File(random_labels, "dummy_output"+(str(time.time()))+".h5","main")
    time_write = time.time()-start_time_write

    total_time = time.time()-start_time

    print("total_time, time_alloc, time_comp, time_write " + str((total_time, time_alloc, time_comp, time_write)))

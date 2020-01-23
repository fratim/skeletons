import math
import struct
import os

OR_Z = 0
OR_Y = 1
OR_X = 2

prefix = 'JWR'
max_label = 136
block_size = [1156, 1156, 1156]
volume_size = [1156*3,1156*3,1156*3]

# get the number of blocks
nzblocks = int(math.ceil(volume_size[OR_Z] // block_size[OR_Z]))
nyblocks = int(math.ceil(volume_size[OR_Y] // block_size[OR_Y]))
nxblocks = int(math.ceil(volume_size[OR_X] // block_size[OR_X]))

# create dictionaries for blocks to synapses
global_synapses_per_block = {}
local_synapses_per_block = {}

# go through all neuron ids
for label in range(1, max_label):
    # read all synapses
    synapse_filename = 'synapses_raw/'+str(prefix)+'/cell'+str(label)+'_t_0_z1_1.txt'

    if not os.path.exists(synapse_filename):
        print("file does not exist: " + synapse_filename)
        continue

    with open(synapse_filename, 'r') as fd:
        for synapse in fd:
            # remove spacing
            synapse = synapse.strip().split()

            global_ix = int(synapse[0])*4
            global_iy = int(synapse[1])*4
            global_iz = int(synapse[2])

            global_iv = global_iz * volume_size[OR_Y] * volume_size[OR_X] + global_iy * volume_size[OR_X] + global_ix

            zblock = global_iz // block_size[OR_Z]
            yblock = global_iy // block_size[OR_Y]
            xblock = global_ix // block_size[OR_X]

            local_iz = global_iz % block_size[OR_Z]
            local_iy = global_iy % block_size[OR_Y]
            local_ix = global_ix % block_size[OR_X]

            local_iv = local_iz * block_size[OR_Y] * block_size[OR_X] + local_iy * block_size[OR_X] + local_ix
            block_key = (zblock, yblock, xblock)

            if not block_key in global_synapses_per_block:
                global_synapses_per_block[block_key] = {}
                local_synapses_per_block[block_key] = {}

            if not label in global_synapses_per_block[block_key]:
                global_synapses_per_block[block_key][label] = []
                local_synapses_per_block[block_key][label] = []

            global_synapses_per_block[block_key][label].append(global_iv)
            local_synapses_per_block[block_key][label].append(local_iv)


for zblock in range(nzblocks):
    for yblock in range(nyblocks):
        for xblock in range(nxblocks):
            block_key = (zblock, yblock, xblock)

            assert (len(global_synapses_per_block[block_key]) == len(local_synapses_per_block[block_key]))

            # get the labels for this block
            neuron_ids = sorted(global_synapses_per_block[block_key].keys())
            nneurons = len(neuron_ids)

            # write synapses for this block
            synapse_filename = 'synapses_{}x{}x{}/{}/{}-synapses-{:04d}z-{:04d}y-{:04d}x.pts'.format(block_size[2], block_size[1], block_size[0], prefix,prefix, zblock, yblock, xblock)

            with open(synapse_filename, 'wb') as fd:
                fd.write(struct.pack('qqqqqq', volume_size[OR_Z], volume_size[OR_Y], volume_size[OR_X], block_size[OR_Z], block_size[OR_Y], block_size[OR_X]))

                # write the number of labels
                fd.write(struct.pack('q', nneurons))

                # write all synapses for each neuron
                for neuron_id in neuron_ids:
                    global_synapses = sorted(global_synapses_per_block[block_key][neuron_id])
                    local_synapses = sorted(local_synapses_per_block[block_key][neuron_id])

                    assert (len(global_synapses) == len(local_synapses))
                    nsynapses = len(global_synapses)

                    # write the label and the number of synapses
                    fd.write(struct.pack('qq', neuron_id, nsynapses))

                    # write the global synapses and then the local ones
                    fd.write(struct.pack('%sq' % nsynapses, *global_synapses))
                    fd.write(struct.pack('%sq' % nsynapses, *local_synapses))




print(volume_size)

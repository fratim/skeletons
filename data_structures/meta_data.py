import sys

from skeletons.utilities.constants import *


class MetaData:
    def __init__(self, prefix):
        # open the meta data txt file
        filename = 'meta/{}.meta'.format(prefix)
        with open(filename, 'r') as fd:
            lines = fd.readlines()

            for ix in range(0, len(lines), 2):
                # get the comment and the corresponding value
                comment = lines[ix].strip()
                value = lines[ix + 1].strip()

                if comment == '# resolution in nm':
                    # separate into individual dimensions
                    samples = value.split('x')
                    # need to use 2, 1, and 0 here since the outward facing convention is x,y,z, not z, y, x
                    self.resolution = (float(samples[2]), float(samples[1]), float(samples[0]))
                elif comment == '# volume size':
                    # separate into individual dimensions
                    samples = value.split('x')
                    # need to use 2, 1, and 0 here since the outward facing convention is x,y,z, not z, y, x
                    self.volume_size = (int(samples[2]), int(samples[1]), int(samples[0]))
                elif comment == '# block size':
                    # separate into individual dimensions
                    samples = value.split('x')
                    # need to use 2, 1, and 0 here since the outward facing convention is x,y,z, not z, y, x
                    self.block_size = (int(samples[2]), int(samples[1]), int(samples[0]))
                elif comment == '# segmentations directory':
                    self.segmentations_directory = str(value)
                elif comment == '# pointclouds directory':
                    self.pointclouds_directory = str(value)
                else:
                    raise ValueError("Unknown keyword in metafile")

    def Resolution(self):
        return self.resolution

    def Volumesize(self):
        return self.volume_size

    def Blocksize(self):
        return self.block_size

    def SegmentationsDirectory(self):
        return self.segmentations_directory

    def PointcloudsDirectory(self):
        return self.pointclouds_directory

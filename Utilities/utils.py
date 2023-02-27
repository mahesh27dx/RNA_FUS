import os
import gsd.hoomd
import numpy as np

def load_data(file):
    data = gsd.hoomd.open(file)
    read_dict = np.load(file[:len(file)-4] + '_radius_dict.npy', allow_pickle='True').item()
    return data, read_dict


class GSDUtils:

    @staticmethod
    def saveSnapshot(filename, gsdSnapshot, hoomdSnapshot):
        """
        WOW: This is super wierd, I couldn't find a good way to copy a
        `gsd.hoomd.Snapshot`. The only difference seems like, there is a
        validate method for the properties and the class itself.
        """
        gsdSnapshot.particles.N = hoomdSnapshot.particles.N
        gsdSnapshot.particles.position = hoomdSnapshot.particles.position
        gsdSnapshot.particles.orientation = hoomdSnapshot.particles.orientation
        gsdSnapshot.particles.typeid = hoomdSnapshot.particles.typeid
        gsdSnapshot.particles.types = hoomdSnapshot.particles.types
        gsdSnapshot.configuration.box = hoomdSnapshot.configuration.box

        # Note: Saving the snapshot
        currentPath = os.path.dirname(__file__)
        snapshotPath = SysUtils.generateSnapshotPathFromUtils(currentPath=currentPath, filename=filename)
        with gsd.hoomd.open(name=snapshotPath, mode='xb') as file:
            file.append(gsdSnapshot)

class SysUtils:

    @staticmethod
    def generateSnapshotPath(currentPath, filename):
        # HACK: Generating the savePath

        if (currentPath == ""):
            file_path = "/snapshots/" + filename
        elif (currentPath != ""):
            file_path = currentPath + '/snapshots/' + filename

        return file_path

    @staticmethod
    def generateSnapshotPathFromUtils(currentPath, filename):

        if (currentPath == ""):
            file_path = "./../snapshots/" + filename
        elif (currentPath != ""):
            file_path = currentPath + '/../snapshots/' + filename

        return file_path

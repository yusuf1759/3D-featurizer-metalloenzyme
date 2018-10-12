import h5py
import numpy as np
from grid_features_metalloenzyme import *
def single_run(dimen,space,finput,center,thresh):
        features_grid_single = GridOccupancy(dimen,space,finput,center,thresh)
        return features_grid_single.FillGrid()
def multi_run(dimen, space, centerFile, thresh,pathtodir):
            f = h5py.File("MetalloenzymeDataset.hdf5", "w")
            with open(centerFile, "r") as fc:
                lines_of_centers = fc.readlines()
            count = 0
            NumFiles = len(lines_of_centers)
            X = np.zeros((NumFiles,7,dimen,dimen,dimen),dtype="f")
            Y = np.zeros((NumFiles,2),dtype=np.int8)
            for m in range(1,len(lines_of_centers)):
                line_split = lines_of_centers[m].split()
                site_id = line_split[0]
                pdb_name = line_split[1]
                Label = line_split[2]
                filePath = str(pathtodir) +'/'+ str(pdb_name) + '_Relaxed.pdb'
                center_of_metal = point(float(line_split[3]),float(line_split[4]),float(line_split[5]))

                X[count] = single_run(dimen,space,filePath,center_of_metal,thresh)
                #print test_count_atomtypes(filePath,center_of_metal, thresh)

                print X[count].sum(axis=(1,2,3))




                if (Label == 'False'):
                        #One-Hot encoding of the label
                        Y[count] = [0,1]

                elif (Label == 'True'):

                        Y[count] = [1,0]

                count += 1
                print (str(count) + " files processed")
            dtrain = f.create_dataset("X", data=X, dtype='f')
            dtrain = f.create_dataset("Y", data=Y, dtype='f')
            f.close()


cFile='/Users/yusuf/Downloads/deeplearning/SampleGridTest/SampleSITECenters.txt'
dirIn='/Users/yusuf/Downloads/deeplearning/SampleGridTest'
multi_run(20,1,cFile,6,dirIn)

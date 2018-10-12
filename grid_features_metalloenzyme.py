#! /usr/bin/env python
#__author__=='Yusuf Adeshina'
#Institution: University of Kansas & Fox Chase Cancer Center, Philadelphia
##Idea of PDB parsing from NNScore, McCammon et al(2011) J. Mol. Graph Model
#Scope:
    #This program takes in a file containing the list of PDB files and the respective center of mass of their binding site metal/metals


import numpy as np
import h5py
import fnmatch
import math
import os
import sys
import textwrap
from collections import Counter
class point:
#Class defining X,Y,Z coordinate
	X=88888.0
    	Y=88888.0
    	Z=88888.0
    
    	def __init__ (self, X, Y ,Z):
        	self.X = X
        	self.Y = Y
        	self.Z = Z

    	def copy_of(self):
        	return point(self.X, self.Y, self.Z)

    	def print_coord(self):
        	print str(self.X)+"\t"+str(self.Y)+"\t"+str(self.Z)
        
        
    	def distance_to(self,apoint):
        	return math.sqrt(math.pow(self.X - apoint.X,2) + math.pow(self.Y - apoint.Y,2) + math.pow(self.Z- apoint.Z,2))

    	def description(self):
        	return str(self.X) + " " + str(self.Y) + " " + str(self.Z)

    	def Magnitude(self):
        	return self.distance_to(point(0,0,0))


class atom:
#Atom class definition
    	def __init__ (self):
		self.RecordName = ""
        	self.AtomName = ""
        	self.ResidueName = ""
        	self.Coordinates = point(88888.0,88888.0,88888.0)
        	self.ElemSym = ""
        	self.PDBIndex = ""
        	self.Charge = 0
        	self.ResId = 0
        	self.ChainId = ""
        	self.Structure = ""
        	self.Comment = ""
        
    	def copy_of(self):
        	atomcopy = atom()
		atomcopy.RecordName = self.RecordName
       		atomcopy.AtomName = self.AtomName 
        	atomcopy.ResidueName = self.ResidueName
        	atomcopy.Coordinates = self.Coordinates.copy_of()
        	atomcopy.ElemSym = self.ElemSym 
        	atomcopy.PDBIndex = self.PDBIndex 
       	 	atomcopy.line= self.line
        	atomcopy.Charge = self.Charge 
        	atomcopy.ResId = self.ResId 
        	atomcopy.ChainId = self.ChainId 
        	atomcopy.Structure = self.Structure 
        	atomcopy.Comment = self.Comment
        
       		return atomcopy

        #Function to read lines from PDB files
	def ReadPDBLine(self, Line):
        	self.line = Line
        	self.AtomName = Line[12:16].strip()      
        	self.Coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))

                if self.ElemSym == "": # try to guess at element from name
            		first_two_letters = self.AtomName[0:2].strip().upper()
	    		if first_two_letters=='BR':
                		self.ElemSym='BR'
	    		elif first_two_letters=='AL':
                		self.ElemSym='AL'
            		elif first_two_letters=='CL':
                		self.ElemSym='CL'
            		elif first_two_letters=='BI':
                		self.ElemSym='BI'
            		elif first_two_letters=='AS':
                		self.ElemSym='AS'
            		elif first_two_letters=='AG':
                		self.ElemSym='AG'
            		elif first_two_letters=='LI':
                		self.ElemSym='LI'
            		elif first_two_letters=='MG':
              			self.ElemSym='MG'
            		elif first_two_letters=='MN':
                		self.ElemSym='MN'
           		elif first_two_letters=='RH':
                		self.ElemSym='RH'
            		elif first_two_letters=='ZN':
                		self.ElemSym='ZN'
            		elif first_two_letters=='FE':
                		self.ElemSym='FE'

            		else: #So,we use just the first letter.
                	# Remove any number from ElemSym
                		self.ElemSym = self.AtomName
                		self.ElemSym = self.ElemSym.replace('0','')
                		self.ElemSym = self.ElemSym.replace('1','')
                		self.ElemSym = self.ElemSym.replace('2','')
                		self.ElemSym = self.ElemSym.replace('3','')
                		self.ElemSym = self.ElemSym.replace('4','')
                		self.ElemSym = self.ElemSym.replace('5','')
                		self.ElemSym = self.ElemSym.replace('6','')
                		self.ElemSym = self.ElemSym.replace('7','')
                		self.ElemSym = self.ElemSym.replace('8','')
                		self.ElemSym = self.ElemSym.replace('9','')
                		self.ElemSym = self.ElemSym.replace('@','')

                		self.ElemSym = self.ElemSym[0:1].strip().upper()

		self.PDBIndex = Line[6:11].strip()
        	self.ResidueName = Line[17:20]
        	self.ResidueName = " " + self.ResidueName[-3:] # uses rightmost three characters, essentially removing unique rotamer identification
        
        	try: self.ResId = int(Line[22:26]) # because it's possible the pdbqt might not have any resid entries.
        	except: pass
        
      	 	self.ChainId = Line[21]
        	if self.ResidueName.strip() == "": self.ResidueName = " MOL"
                self.RecordName = Line[0:6]

class Loader:

    	def __init__ (self):
        	self.AllAtoms = {}
		self.TranslatedPoints= {}
		self.TranslatedPointsOnGrid = {}
        	self.NonProteinAtoms = {}
        	self.max_x = -8888.88
        	self.min_x = 8888.88
        	self.max_y = -8888.88
        	self.min_y = 8888.88
        	self.max_z = -8888.88
        	self.min_z = 8888.88

        	self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
     
	def PDBLoad(self, FileName, dim, spacing,metal_com,cutoff, min_x=-8888.88, max_x=8888.88, min_y=-8888.88, max_y=8888.88, min_z=-8888.88, max_z=8888.88):#cutoff

        	autoindex = 1

        	self.__init__()
        
        	# Now load the file into a list
        	file = open(FileName,"r")
        	lines = file.readlines()
        	file.close()
        
        	atom_already_loaded = [] # going to keep track of atomname_resid_chain pairs, to make sure redundants aren't loaded. This basically
                                 # gets rid of rotomers, I think.
		print FileName
		for t in range(0,len(lines)):
			#print "OK"
            		line=lines[t]
            
            		if line[:3] == "END":
				#print "OK"
                		t = textwrap.wrap("WARNING: END or ENDMDL encountered in " + FileName + ". Everyline after this will be ignored. If your PDB file has multiple pose or is an ensemble structure, split them up into individual pose or model", 80)
                		print "\n".join(t) + "\n"
                		print line
               			break
                    
            		if len(line) >= 7:
				#print "OK"
                		if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    			TempAtom = atom()
                   	 		TempAtom.ReadPDBLine(line)
                        
                                        
                                        if (TempAtom.Coordinates.distance_to(metal_com) <= cutoff):
                                                #Load only atoms withing the cutoff range
                                                if self.max_x < TempAtom.Coordinates.X: self.max_x = TempAtom.Coordinates.X
                        			if self.max_y < TempAtom.Coordinates.Y: self.max_y = TempAtom.Coordinates.Y
                        			if self.max_z < TempAtom.Coordinates.Z: self.max_z = TempAtom.Coordinates.Z
                        
                        			if self.min_x > TempAtom.Coordinates.X: self.min_x = TempAtom.Coordinates.X
                        			if self.min_y > TempAtom.Coordinates.Y: self.min_y = TempAtom.Coordinates.Y
                        			if self.min_z > TempAtom.Coordinates.Z: self.min_z = TempAtom.Coordinates.Z

                        			key = TempAtom.AtomName.strip() + "_" + str(TempAtom.ResId) + "_" + TempAtom.ResidueName.strip() + "_" + TempAtom.ChainId.strip() # this string unique identifies each atom
                        
                        			if key in atom_already_loaded and TempAtom.ResidueName.strip() in self.protein_resnames: # so this is a protein atom that has already been loaded once
                            				self.printout("Warning: Duplicate protein atom detected: \"" + TempAtom.line.strip() + "\". Not loading this duplicate.")
                            				print ""
                        
                        			if not key in atom_already_loaded or not TempAtom.ResidueName.strip() in self.protein_resnames: # so either the atom hasn't been loaded, or else it's a non-protein atom
                                                                                                            # so note that non-protein atoms can have redundant names, but protein atoms cannot.
                                                                                                            # This is because protein residues often contain rotamers
                            				atom_already_loaded.append(key) # so each atom can only be loaded once. No rotamers.
                            				self.AllAtoms[autoindex] = TempAtom # So you're actually reindexing everything here.
                            				if not TempAtom.ResidueName[-3:] in self.protein_resnames: self.NonProteinAtoms[autoindex] = TempAtom;#print "OK"
                            				
                           				autoindex = autoindex + 1
                print autoindex                            

		for iatom in self.AllAtoms:
			iTempPoint = point(88888.0,88888.0,88888.0) #Note:Always create instance inside the loop, if created outside the value will be updated to the last value and the all will be the same.
			iTempPoint.X = self.AllAtoms[iatom].Coordinates.X - metal_com.X
			iTempPoint.Y = self.AllAtoms[iatom].Coordinates.Y - metal_com.Y
			iTempPoint.Z = self.AllAtoms[iatom].Coordinates.Z - metal_com.Z
			self.TranslatedPoints[iatom] = iTempPoint
	    
		for katom in self.TranslatedPoints:
			#print self.TranslatedPoints[katom].X
			kTempPoint = point(88888.0,88888.0,88888.0)
     			#kTempPoint.X = (float(self.TranslatedPoints[katom].X) + 0.5*(dim -1)*spacing) - metal_com.X
                        kTempPoint.X = float(self.TranslatedPoints[katom].X)
		        #kTempPoint.Y = (float(self.TranslatedPoints[katom].Y) + 0.5*(dim -1)*spacing) - metal_com.Y
                        kTempPoint.Y = float(self.TranslatedPoints[katom].Y)
                        #kTempPoint.Z = (float(self.TranslatedPoints[katom].Z) + 0.5*(dim -1)*spacing) - metal_com.Z
                        kTempPoint.Z = float(self.TranslatedPoints[katom].Z)
			self.TranslatedPointsOnGrid[katom] = kTempPoint
                        #self.TranslatedPointsOnGrid[katom] = TempPoint
		
	def printout(self, thestring):
        	lines = textwrap.wrap(thestring, 80)
        	for line in lines:
            		print line


class GridOccupancy:
#Build the grid and put the atoms in it. I used 7 atomtypes C,N,O,S,P,Metals,Others

	def __init__ (self,dim,spacing,complex_pdb,CoM,threshold):
		self.dim = dim
		self.spacing = spacing
		self.complex_pdb = complex_pdb
                self.threshold = threshold
		#self.GridData = np.zeros((self.dim,self.dim,self.dim,6))
                self.GridData = np.zeros((7,self.dim,self.dim,self.dim))
		self.Complex_PDB = Loader()
		self.Complex_PDB.PDBLoad(complex_pdb,self.dim,self.spacing,CoM,self.threshold )


	def FillGrid(self):

		regrid = [i - self.dim/2 for i in range (self.dim) ]
    		print regrid
   		
    
		for latom in self.Complex_PDB.AllAtoms:
		
        		if (math.floor(np.abs(self.Complex_PDB.TranslatedPointsOnGrid[latom].X)) < self.dim/2 and math.floor(np.abs(self.Complex_PDB.TranslatedPointsOnGrid[latom].Y)) < self.dim/2  and math.floor(np.abs(self.Complex_PDB.TranslatedPointsOnGrid[latom].Z)) < self.dim/2 ):
            				x_map = math.floor(self.Complex_PDB.TranslatedPointsOnGrid[latom].X)
           	 			y_map = math.floor(self.Complex_PDB.TranslatedPointsOnGrid[latom].Y)
            				z_map = math.floor(self.Complex_PDB.TranslatedPointsOnGrid[latom].Z)
					#print "OK"
					#print self.Complex_PDB.AllAtoms[latom].ElemSym
            				if (self.Complex_PDB.AllAtoms[latom].ElemSym == "C"):
                                                #print 'ok1'
                				#self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),0] = 1
                                                self.GridData[0,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
						#print self.Complex_PDB.AllAtoms[latom].ElemSym
            				elif (self.Complex_PDB.AllAtoms[latom].ElemSym == "N"):
                                                #print 'ok2'
                				#self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),1] = 1
                                                self.GridData[1,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
            				elif (self.Complex_PDB.AllAtoms[latom].ElemSym == "O"):
                                                #print 'ok3'
                				#self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),2] = 1
                                                self.GridData[2,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
                                        elif (self.Complex_PDB.AllAtoms[latom].ElemSym == "S"):
                                                #print 'ok4'
                                                #self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),3] = 1
                                                self.GridData[3,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
                                        elif (self.Complex_PDB.AllAtoms[latom].ElemSym == "P"):
                                                #print 'ok5'
                                                #self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),4] = 1
                                                self.GridData[4,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
                                        elif (self.Complex_PDB.AllAtoms[latom].ElemSym in ["ZN","FE","MN","MG","AU","AG","CO","LI","BI"]):
                                                #print 'ok6'
                                                #self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),4] = 1
                                                self.GridData[5,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
            				#elif (self.Complex_PDB.AllAtoms[latom].ElemSym == "Other"):
                				#self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),3] = 1
            				else:
                				#print ("Atom type:", self.Complex_PDB.AllAtoms[latom].ElemSym, "not recognized. Please fix it, or update the atomtype list")
                				#quit()
						#self.GridData[regrid.index(x_map),regrid.index(y_map),regrid.index(z_map),5] = 1
                                                #print 'ok7'
                                                self.GridData[6,regrid.index(x_map),regrid.index(y_map),regrid.index(z_map)] = 1
		return self.GridData
'''   
#Function to featurize a single protein
def single_run(dimen,space,finput,center,thresh):
	features_grid_single = GridOccupancy(dimen,space,finput,center,thresh)
	return features_grid_single.FillGrid()

#Function to featurize multiple protein
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
        
def test_count_atomtypes(pdbfile,cent,cutof):
    with open(pdbfile, 'r') as pdb:
        lines = pdb.readlines()
    atm_count =[]
    for line in lines:
        if line[0:4]=="ATOM" or line[0:6]=="HETATM": 
            tAtom = atom()
            tAtom.ReadPDBLine(line)
            
            if (tAtom.Coordinates.distance_to(cent) <= cutof):
                atm_count.append(tAtom.ElemSym)
                #print tAtom.Coordinates.distance_to(cent)
    return Counter(atm_count)
#cFile='/Users/yusuf/Downloads/deeplearning/SampleGridTest/SampleSITECenters.txt'
#dirIn='/Users/yusuf/Downloads/deeplearning/SampleGridTest'
#multi_run(20,1,cFile,6,dirIn)
'''

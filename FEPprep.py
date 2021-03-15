#!/usr/bin/env python

import argparse
import os
import sys
import copy
import random
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS, AllChem, rdShapeHelpers, rdmolops, rdForceFieldHelpers
from rdkit.Chem.Fingerprints import FingerprintMols



def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        print("Error: filename %s could not be found: exiting." % (string))



## Function to read in the sdf- or smiles-file
def readMol(filePath):
	fileName = filePath.split('/')[-1]
	fileNameExt = fileName.split('.')[-1]
	fileNameExtLower = fileNameExt.lower()[:3]
    
	if (fileNameExtLower == "sdf"):
		try:
			mol = Chem.MolFromMolFile(filePath, removeHs=False)
		except:
			print("Error reading molecule %s" % (filePath))
		if mol is None or mol == "":
			print("Bad molecule %s" % (filePath))
			sys.exit(1)
	elif (fileNameExtLower == "smi"):
		smilesFile = open(filePath, 'r')
		line = smilesFile.readlines()[0]
		line = line.strip()
		if line is None or line == "":
			print("Empty first line in file %s" % (filePath))
			sys.exit(1)
		SMILES = line.split()[0]
		mol = Chem.MolFromSmiles(SMILES)
		if mol is None or mol == "":
			print("Bad molecule %s" % (filePath))
			sys.exit(1)
		mol = Chem.AddHs(mol)
		smilesFile.close()       
	else:
		print("File format should be of type smiles or sdf")
		sys.exit(1) 
	return mol



# Function to read the reference structure
def readRefMol(filePath):
    fileName = filePath.split('/')[-1]
    fileNameExt = fileName.split('.')[-1]
    fileNameExtLower = fileNameExt.lower()[:3]
    
    if (fileNameExtLower == "sdf"):
        sdfFile = Chem.SDMolSupplier(filePath, removeHs=False)
        if len(sdfFile) > 0: mol = (sdfFile[0])
        if mol is None or mol == "":
            print("Invalid reference structure")
            sys.exit(1)   
    else:
        print("File format should be of type sdf")
        sys.exit(1)   
    return mol
 


# Argparse code and config
parser = argparse.ArgumentParser(description="A script to superimpose molecules onto a reference structure using a MCSS search")
parser.add_argument("-r", "--ref", 
					metavar="", 
					help="Reference file (.sdf)", 
					required=True, 
					type=file_path)
parser.add_argument("-i", "--input",
					metavar="", 
					help="Input file (.sdf or .smi or .smiles)", 
					required=True, 
					type=file_path)
parser.add_argument("-o", "--out", 
					metavar="", 
					help="Output file (.sdf)", 
					required=True, 
					type=str)
parser.add_argument("-t", "--time", 
					metavar="", 
					help="Maximum time (in sec) allowed to find the MCSS [default: 3 sec]", 
					required=False, 
					type=int, 
					default=3)
parser.add_argument("-m", "--mcss", 
					metavar="", 
					help="Minimum number of MCSS atoms required [default: 3]", 
					required=False, 
					type=int,
					default=3)

args = parser.parse_args()

refPath = args.ref
targetPath = args.input
finalOutPath = args.out
timeoutTime = args.time
MCSSThreshold = args.mcss


## Reading and prepping in the reference molecule and the target
refMol = readRefMol(refPath)
targetMol = readMol(targetPath)
SDwriter = Chem.SDWriter(finalOutPath)

		
# Generate conformation for target (SMILES)
if os.path.splitext(targetPath)[1][1:] == "smi" or os.path.splitext(targetPath)[1][1:] == "smiles":
    AllChem.EmbedMolecule(targetMol)
AllChem.EmbedMolecule(targetMol)


# Molecule name
if not '.smi' in targetPath.lower():
    if targetMol.GetProp("_Name") is None or targetMol.GetProp("_Name") == "": 
        targetName = Chem.MolToSmiles(targetMol)
    else:
        targetName = targetMol.GetProp("_Name")
else:
    targetName = Chem.MolToSmiles(targetMol)
 
    
# Searching substructure and converting it to mol object
print("Processing molecule: " + targetName)
print("Searching maximum common substructure... ")
MCSresult = rdFMCS.FindMCS([Chem.RemoveHs(refMol), targetMol], 
                           ringMatchesRingOnly=False,  
                           matchValences=True, 
                           completeRingsOnly=False, 
                           timeout=timeoutTime, 
                           bondCompare=rdFMCS.BondCompare.CompareOrder, 
                           atomCompare=rdFMCS.AtomCompare.CompareElements,
						   matchChiralTag=True,
                           maximizeBonds=False)
MCSsmarts = MCSresult.smartsString
print(MCSsmarts)
MCSmol = Chem.MolFromSmarts(MCSsmarts)
if MCSSThreshold >= MCSresult.numAtoms:
    print("Skipping this molecule, MCSS below the cutoff of %d atoms. Exit." % (MCSSThreshold))
    sys.exit()
   
    
# Matching atom Idx to complementary structure
refMatchIdx = refMol.GetSubstructMatch(MCSmol, useChirality=True)
targetMatchIdx = targetMol.GetSubstructMatch(MCSmol, useChirality=True)
mcssMatchIdx = MCSmol.GetSubstructMatch(MCSmol, useChirality=True)
print(len(refMatchIdx), len(targetMatchIdx), len(mcssMatchIdx))
atomMap = []
for i in range(len(refMatchIdx)): atomMap.append([targetMatchIdx[i], refMatchIdx[i]])


# Transfer coordinates of matching atoms
AllChem.EmbedMolecule(MCSmol)
for ri, ti, si in zip(refMatchIdx, targetMatchIdx, mcssMatchIdx):   
	targetMol.GetConformer().SetAtomPosition(ti, refMol.GetConformer().GetAtomPosition(ri))
	MCSmol.GetConformer().SetAtomPosition(si, refMol.GetConformer().GetAtomPosition(ri))


# Remove hydrogens and add them back
targetMol = Chem.AddHs(Chem.RemoveHs(targetMol), addCoords=True)
refMatchIdx = refMol.GetSubstructMatch(MCSmol, useChirality=True)
targetMatchIdx = targetMol.GetSubstructMatch(MCSmol, useChirality=True)
atomMap = []
for i in range(len(refMatchIdx)): atomMap.append([targetMatchIdx[i], refMatchIdx[i]])


# Restrained minimization
ffProps = AllChem.MMFFGetMoleculeProperties(targetMol, mmffVariant='MMFF94')
ff = AllChem.MMFFGetMoleculeForceField(targetMol, ffProps)
for i in targetMatchIdx: ff.AddFixedPoint(i)
if ff.Minimize(1000) != 0:	  
	print("Error in minimizing: exit")
	sys.exit(1)


# Stereochemistry perceptions
Chem.rdmolops.AssignStereochemistryFrom3D(targetMol, replaceExistingTags=True)


# Write out
SDwriter.write(targetMol)
SDwriter.close()
print("Written file to " + finalOutPath)

sys.exit(0)   

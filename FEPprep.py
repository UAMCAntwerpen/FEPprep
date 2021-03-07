#!/usr/bin/env python

import argparse
import os
import sys
import copy
import random

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFMCS, AllChem, rdShapeHelpers, rdmolops
from rdkit.Chem.Fingerprints import FingerprintMols



def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{string} could not be found")



## Function to read in the sdf- or smiles-file
def readMol(filePath):
    fileName = filePath.split('/')[-1]
    fileNameExt = fileName.split('.')[-1]
    fileNameExtLower = fileNameExt.lower()[:3]
    molInputList = []
    
    if (fileNameExtLower == "sdf"):
        sdfFile = Chem.SDMolSupplier(filePath, removeHs=False)
        for mol in sdfFile:
            if mol is None or mol == "": continue
            molInputList.append(mol)
    elif (fileNameExtLower == "smi"):
        smilesFile = open(filePath, 'r')
        smilesLines = smilesFile.readlines()
        for line in smilesLines:
            line = line.strip()
            if line is None or line == "": continue
            SMILES = line.split()[0]
            mol = Chem.MolFromSmiles(SMILES)
            if mol is None or mol == "": continue
            mol = Chem.AddHs(mol)
            molInputList.append(mol)
        smilesFile.close()       
    else:
        print("File format should be of type smiles or sdf")
        sys.exit(1)  
    return molInputList



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
parser.add_argument("-c", "--conf", 
					metavar="", 
					help="Number of conformers to generate [default: 50]", 
					required=False, 
					type=int, 
					default=50)
parser.add_argument("-t", "--time", 
					metavar="", 
					help="Maximum time (in sec) allowed to find the MCSS [default: 5 sec]", 
					required=False, 
					type=int, 
					default=5)
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
numConfs = args.conf
timeoutTime = args.time
MCSSThreshold = args.mcss



## Reading and prepping in the reference molecule and the targetSMILES
refMol = readRefMol(refPath)
molInputList = readMol(targetPath)

SDwriter = Chem.SDWriter(finalOutPath)

for targetMol in molInputList:
    
    embedError = False

    if targetMol.GetProp("_Name") is None or targetMol.GetProp("_Name") == "": 
        targetName = Chem.MolToSmiles(targetMol)
    else:
        targetName = targetMol.GetProp("_Name")
    
    targetMol = Chem.AddHs(targetMol, addCoords=True)
    molList = [refMol, targetMol]
    
    ## Searching substructure and converting it to mol object
    print("Processing molecule: " + targetName)
    print("Searching maximum common substructure... ")
    MCSresult = rdFMCS.FindMCS(molList, ringMatchesRingOnly=True, matchValences=True, completeRingsOnly=True, timeout=timeoutTime)
    complSMARTS = MCSresult.smartsString
    complMol = Chem.MolFromSmarts(complSMARTS)
    complConf = Chem.Conformer(complMol.GetNumAtoms())
    
    print("Found a maximum common substructure (MCSS) with " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds!")
    if MCSSThreshold >= complMol.GetNumAtoms() :
        print("skipping this molecule, MCSS below the cutoff of %d atoms!" % (MCSSThreshold))
        print("---------------------------------------------------------------------------------\n")
        continue
		
    ## Generate Conformation for target (SMILES)
    AllChem.EmbedMolecule(targetMol)
    
    ## Matching atom Idx to complementary structure
    refMatch = refMol.GetSubstructMatch(complMol)
    targetMatch = targetMol.GetSubstructMatch(complMol)
    
    ## Make a list of complMol atomIdxs to be able to loop over
    complIdxs = []
    for atom in complMol.GetAtoms():
        complIdxs.append(atom.GetIdx())
    
    ## Loop through the atomIdxs matching the substructure of targetMol and complMol, 
    ## to set their positions equal to the exact position of the corresponding atom of the refMol
    print("Generating " + str(numConfs) + " conformations.")
    for refIdx, targetIdx, complIdx in zip(refMatch, targetMatch, complIdxs):   
          refCoords = refMol.GetConformer().GetAtomPosition(refIdx)
          targetMol.GetConformer().SetAtomPosition(targetIdx, refCoords)
          complConf.SetAtomPosition(complIdx, refCoords)
    complMol.AddConformer(complConf) 
    
    ## Constrained embed numConf amount of times, with restraint core being the MCSS
    molList = []
    
    for i in range(numConfs):
        seed = random.randint(0, 9999999)
        tempMol = copy.copy(targetMol)
        try:
            molConf = AllChem.ConstrainedEmbed(tempMol, complMol, useTethers=False, randomseed=seed) 
        except:
            print("Error occured while embedding " + targetName + ", moving on to the next molecule.")
            print("---------------------------------------------------------------------------------\n")
            embedError= True
            break  
        
    if embedError == True:
        continue
    
    molList.append(molConf)
    
    
    shapeList = []
    
    ## Re-adjust the MCSS coordinates to exactly match the refMol before calculating shape metric
    print("Calculating best match...")
    for i, mol in enumerate(molList):
        
        for refIdx, targetIdx in zip(refMatch, targetMatch): 
            
            refCoords = refMol.GetConformer().GetAtomPosition(refIdx)
            mol.GetConformer().SetAtomPosition(targetIdx, refCoords)
            
        ## Calculate volume metric (shape protrude distance) of conf vs refmol and keep track of the smallest value
        prtrd = rdShapeHelpers.ShapeProtrudeDist(mol, refMol)
    
        if len(shapeList) < 1:
            shapeList = [i, prtrd]
        else:
            if shapeList[1] < prtrd:
                continue
            else:
                shapeList = [i, prtrd]
    
    ## Write the file to the output location
    molIndex = shapeList[0]
    finalMol = molList[molIndex]
    
    #Calculating tanimoto distance and formal charge
    fpRef = FingerprintMols.FingerprintMol(refMol)
    fpFinal = FingerprintMols.FingerprintMol(finalMol)
    taniSim = DataStructs.FingerprintSimilarity(fpRef, fpFinal)
    
    formalChargeFinal = rdmolops.GetFormalCharge(finalMol)
    formalChargeRef = rdmolops.GetFormalCharge(refMol)
    
    ## Give additional information as output
    
    print("The final output molecule has a formal charge of " + str(formalChargeFinal) + ", and contains " + str(finalMol.GetNumAtoms()) +" atoms and " +  str(finalMol.GetNumBonds())+ " bonds.")
    print("The reference molecule has a formal charge of " + str(formalChargeRef) + ", and contains " + str(refMol.GetNumAtoms()) +" atoms and " +  str(refMol.GetNumBonds()) + " bonds.")
    print("The Tanimoto Similarity between the final structure and the reference compound is " + str(taniSim))
    print("The maximum common substructure (MCSS) contained " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds")
    print("---------------------------------------------------------------------------------\n")
    
    finalMol.SetProp("_Name", targetName)
    finalMol.SetProp("Tanimoto", str(taniSim))
    
    SDwriter.write(finalMol)

SDwriter.close()
print("Written file to " + finalOutPath)
sys.exit(0)

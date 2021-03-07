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



## ArgParse Code and config
def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{string} is not a valid path")

def dir_path(string):
    
    dirName = os.path.dirname(string)
    fileName = os.path.basename(string)

    if os.path.isdir(dirName):
        if "." in fileName:
            
            return string
        else:
            if fileName.strip() == '':
                raise argparse.ArgumentTypeError("The output path needs to contain a filename")
            else:   
                fileName = fileName + ".sdf"
                string = os.path.join(dirName, fileName)
                return(string)
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{string} is not a valid path")
  

parser = argparse.ArgumentParser(description="A script to ")
parser.add_argument("-r","--ref", metavar="", help="Path of the reference file", required=True, type=file_path)
parser.add_argument("-s","--smiles", "--sdf", "--mol", metavar="", help="Path of the input SMILES, mol or sdf-file", required=True, type=file_path)
parser.add_argument("-o","--out", metavar="", help="Path of the output file", required=True, type=dir_path)
parser.add_argument("-c","--conf", metavar="", help="Amount of conformers to generate. Default is 50", required=False, type=int, default=50)
parser.add_argument("-t","--time", metavar="", help="Amount of time spent to find the maximum common substructure. Default is 5 seconds", required=False, type=int, default=5)
parser.add_argument("-m","--mcss", "--MCSS", metavar="", help="Amount of atoms that need to overlap before continuing the claculations", required=False, type=int, default=3)

args = parser.parse_args()

refPath = args.ref
targetPath = args.smiles
finalOutPath = args.out
numConfs = args.conf
timeoutTime = args.time
MCSSThreshold = args.mcss


## Function to read in a mol- or SDF-file
def readMol(filePath):
    fileName = filePath.split('/')[-1]
    fileNameExt = fileName.split('.')[-1]
    fileNameExtLower = fileNameExt.lower()
    molInputList = []
    
    if(fileNameExtLower == "mol"):
        mol = Chem.MolFromMolFile(filePath, removeHs=False)
        molInputList.append(mol)
    elif(fileNameExtLower == "sdf"):
        sdfFile = Chem.SDMolSupplier(filePath, removeHs=False)
        mol = sdfFile
    elif(fileNameExtLower == "smi" or fileNameExtLower == "smiles"):
        smilesFile = open(filePath, 'r')
        smilesLines = smilesFile.readlines()
        
        for line in smilesLines:
            SMILES = line.split()[0]
            molInputList.append(Chem.MolFromSmiles(SMILES))
        smilesFile.close()       
    else:
        print("did not recognize file format")
        sys.exit()  
    return molInputList

def readRefMol(filePath):
    fileName = filePath.split('/')[-1]
    fileNameExt = fileName.split('.')[-1]
    fileNameExtLower = fileNameExt.lower()
    
    if(fileNameExtLower == "mol"):
        mol = (Chem.MolFromMolFile(filePath))
    elif(fileNameExtLower == "sdf"):
        sdfFile = Chem.SDMolSupplier((filePath))
        mol = (sdfFile[0])
    else:
        print("did not recognize file format")
        sys.exit()   
    return mol

## Reading and prepping in the reference molecule and the targetSMILES
refMol = readRefMol(refPath)
molInputList = readMol(targetPath)

SDwriter = Chem.SDWriter(finalOutPath)

for targetMol in molInputList:
    
    embedError = False
    
    targetName = Chem.MolToSmiles(targetMol)
    
    targetMol = Chem.AddHs(targetMol, addCoords=True)
    molList = [refMol, targetMol]
    
    ## Searching substructure and converting it to mol object
    print("Processing molecule with SMILES: " + targetName)
    print("Searching maximum common substructure... ")
    MCSresult = rdFMCS.FindMCS(molList, ringMatchesRingOnly=True, matchValences=True, completeRingsOnly=True, timeout=timeoutTime)
    complSMARTS = MCSresult.smartsString
    complMol = Chem.MolFromSmarts(complSMARTS)
    complConf = Chem.Conformer(complMol.GetNumAtoms())
    
    print("Found a maximum common substructure (MCSS) with " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds!")
    if MCSSThreshold >= complMol.GetNumAtoms() :
        print("skipping this molecule, MCSS too little!")
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
    print("Generating " + str(numConfs) + " conformations...")
    
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

import argparse
import os
import sys
import copy
import random

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, rdShapeHelpers, rdmolops

## ArgParse Code and config
def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{string} is not a valid path")

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{string} is not a valid path")
  
def file_name(string):
    if ".mol" in string.lower() or ".sdf" in string.lower():
        return string
    else:
        string = string + '.sdf'
        return string

parser = argparse.ArgumentParser(description="A script to ")
parser.add_argument("-r","--ref", metavar="", help="Path of the reference file", required=True, type=file_path)
parser.add_argument("-t","--target", metavar="", help="SMILES input of the target molecule", required=True)
parser.add_argument("-o","--out", metavar="", help="Path of the output file", required=True, type=dir_path)
parser.add_argument("-n","--name", metavar="", help="Name of the output file", required=True, type=file_name)
parser.add_argument("-c","--conf", metavar="", help="Amount of conformers to generate. Default is 50", required=False, type=int, default=50)

args = parser.parse_args()

refPath = args.ref
targetSMILES = args.target
outPath = args.out
fileName = args.name
finalOutPath = os.path.join(outPath, fileName)
numConfs = args.conf

## Function to read in a mol- or SDF-file
def readMol(filePath):
    fileName = filePath.split('/')[-1]
    fileNameExt = fileName.split('.')[-1]
    fileNameExtLower = fileNameExt.lower()
    
    if(fileNameExtLower == "mol"):
        mol = Chem.MolFromMolFile(filePath, removeHs=False)
    elif(fileNameExtLower == "sdf"):
        sdfFile = Chem.SDMolSupplier(filePath, removeHs=False)
        mol = sdfFile[0]
    else:
        print("did not recognize file format")
        sys.exit()     
    return mol

## Reading and prepping in the reference molecule and the targetSMILES
refMol = readMol(refPath)
targetMol = Chem.MolFromSmiles(targetSMILES)
targetMol = Chem.AddHs(targetMol, addCoords=True)
molList = [refMol, targetMol]

## Searching substructure and converting it to mol object
print("Searching largest common subsctructure...")
MCSresult = rdFMCS.FindMCS(molList, ringMatchesRingOnly=True, matchValences=True, completeRingsOnly=True, timeout=30)
complSMARTS = MCSresult.smartsString

complMol = Chem.MolFromSmarts(complSMARTS)
complConf = Chem.Conformer(complMol.GetNumAtoms())

print("Found a maximum common substructure with " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds!")

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
    molConf = AllChem.ConstrainedEmbed(tempMol, complMol, useTethers=False, randomseed=seed)  
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

SDwriter = Chem.SDWriter(finalOutPath)
SDwriter.write(finalMol)
SDwriter.close()
print("Written file to " + finalOutPath)

## Give additional information as output
formalCharge = rdmolops.GetFormalCharge(finalMol)
print("The output molecule has a formal charge of " + str(formalCharge) + ", and contains " + str(finalMol.GetNumAtoms()) +" atoms and " +  str(finalMol.GetNumBonds())+ " bonds.")
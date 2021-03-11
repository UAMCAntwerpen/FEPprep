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
		print("Error: filename %s could not be found: exiting." % (string))



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
 


# Function to get all torsions in a molecule
def enumerateTorsions(mol):
    torsionSmarts = '[!$(*#*)&!D1]~[!$(*#*)&!D1]'
    torsionQuery = Chem.MolFromSmarts(torsionSmarts)
    matches = mol.GetSubstructMatches(torsionQuery)
    torsionList = []
    for match in matches:
        idx2 = match[0]
        idx3 = match[1]
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        jAtom = mol.GetAtomWithIdx(idx2)
        kAtom = mol.GetAtomWithIdx(idx3)
        if (((jAtom.GetHybridization() != Chem.HybridizationType.SP2)
            and (jAtom.GetHybridization() != Chem.HybridizationType.SP3))
            or ((kAtom.GetHybridization() != Chem.HybridizationType.SP2)
            and (kAtom.GetHybridization() != Chem.HybridizationType.SP3))): continue
        for b1 in jAtom.GetBonds():
            if (b1.GetIdx() == bond.GetIdx()): continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in kAtom.GetBonds():
                if ((b2.GetIdx() == bond.GetIdx()) or (b2.GetIdx() == b1.GetIdx())): continue
                idx4 = b2.GetOtherAtomIdx(idx3)
                # skip 3-membered rings
                if (idx4 == idx1): continue
                torsionList.append((idx1, idx2, idx3, idx4))
    return torsionList
 


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
    if not '.smi' in targetPath.lower():
        if targetMol.GetProp("_Name") is None or targetMol.GetProp("_Name") == "": 
            targetName = Chem.MolToSmiles(targetMol)
        else:
            targetName = targetMol.GetProp("_Name")
    else:
        targetName = Chem.MolToSmiles(targetMol)

    targetMol = Chem.AddHs(targetMol, addCoords=True)
    
    ## Searching substructure and converting it to mol object
    print("Processing molecule: " + targetName)
    print("Searching maximum common substructure... ")
    MCSresult = rdFMCS.FindMCS([Chem.RemoveHs(refMol), targetMol], 
                               ringMatchesRingOnly=False,  
                               matchValences=True, 
                               completeRingsOnly=False, 
                               timeout=timeoutTime, 
                               bondCompare=rdFMCS.BondCompare.CompareOrder, 
                               atomCompare=rdFMCS.AtomCompare.CompareElements,
                               maximizeBonds=False)
    complSMARTS = MCSresult.smartsString
    complMol = Chem.MolFromSmarts(complSMARTS)
    complConf = Chem.Conformer(complMol.GetNumAtoms())

    if MCSSThreshold >= complMol.GetNumAtoms() :
        print("skipping this molecule, MCSS below the cutoff of %d atoms!" % (MCSSThreshold))
        continue
		
    ## Generate Conformation for target (SMILES)
    if os.path.splitext(targetPath)[1][1:] == "smi" or os.path.splitext(targetPath)[1][1:] == "smiles":
        AllChem.EmbedMolecule(targetMol)
    
    ## Matching atom Idx to complementary structure
    refMatch = refMol.GetSubstructMatch(complMol)
    targetMatch = targetMol.GetSubstructMatch(complMol)
    
    # Align the torsions of the target molecule with those of the reference
    refTorsions = []
    torsions = enumerateTorsions(refMol)
    for torsion in torsions:
        if torsion[0] in refMatch and torsion[1] in refMatch and torsion[2] in refMatch and torsion[3] in refMatch:
            if refMol.GetAtomWithIdx(torsion[1]).IsInRing() and refMol.GetAtomWithIdx(torsion[2]).IsInRing(): continue
            x = list(torsion)
            refTorsions.append(x)
 
    targetTorsions = []
    torsions = enumerateTorsions(targetMol)
    for torsion in torsions:
        if torsion[0] in targetMatch and torsion[1] in targetMatch and torsion[2] in targetMatch and torsion[3] in targetMatch:
            x = list(torsion)
            targetTorsions.append(x)
            targetTorsions.append([i for i in reversed(x)])
    
    t2r = {}
    atomMap = []
    for i in range(len(refMatch)):
        t2r[targetMatch[i]] = refMatch[i]
        atomMap.append([targetMatch[i], refMatch[i]])

    for targetTorsion in targetTorsions:
        a = t2r[targetTorsion[0]]
        b = t2r[targetTorsion[1]]
        c = t2r[targetTorsion[2]]
        d = t2r[targetTorsion[3]]
        for refTorsion in refTorsions:
            if a == refTorsion[0] and b == refTorsion[1] and c == refTorsion[2] and d == refTorsion[3]:
                angle = AllChem.GetDihedralDeg(refMol.GetConformer(), a, b, c, d)
                AllChem.SetDihedralDeg(targetMol.GetConformer(), 
                                       targetTorsion[0],
                                       targetTorsion[1],
                                       targetTorsion[2],
                                       targetTorsion[3], angle)

    # Align the target onto the reference
    AllChem.AlignMol(targetMol, refMol, atomMap = atomMap)    

    # Make a list of complMol atomIdxs to be able to loop over
    complIdxs = []
    for atom in complMol.GetAtoms(): complIdxs.append(atom.GetIdx())
    
    # Loop through the atomIdxs matching the substructure of targetMol and complMol, 
    # to set their positions equal to the exact position of the corresponding atom of the refMol
    for refIdx, targetIdx, complIdx in zip(refMatch, targetMatch, complIdxs):   
          refCoords = refMol.GetConformer().GetAtomPosition(refIdx)
          #targetMol.GetConformer().SetAtomPosition(targetIdx, refCoords)
          complConf.SetAtomPosition(complIdx, refCoords)
    complMol.AddConformer(complConf) 
    
    # Constrained embed numConf amount of times, with restraint core being the MCSS
    molList = []
    
    print("Generating " + str(numConfs) + " conformations.")
    for i in range(numConfs):
        seed = random.randint(0, 9999999)
        tempMol = copy.copy(targetMol)
        try:
            molConf = AllChem.ConstrainedEmbed(tempMol, complMol, useTethers=False, randomseed=seed) 
        except:
            print("Error occured while embedding " + targetName)
            embedError= True
            continue
        molList.append(molConf)
        
    if embedError == True: sys.exit()
    
    shapeList = []
    
    # Re-adjust the MCSS coordinates to exactly match the refMol before calculating shape metric
    print("Calculating best match...")
    for i, mol in enumerate(molList):
        
        for refIdx, targetIdx in zip(refMatch, targetMatch): 
            
            refCoords = refMol.GetConformer().GetAtomPosition(refIdx)
            mol.GetConformer().SetAtomPosition(targetIdx, refCoords)
            
        # Calculate volume metric (shape protrude distance) of conf vs refmol and keep track of the smallest value
        prtrd = rdShapeHelpers.ShapeProtrudeDist(mol, refMol)
    
        if len(shapeList) < 1:
            shapeList = [i, prtrd]
        else:
            if shapeList[1] < prtrd:
                continue
            else:
                shapeList = [i, prtrd]
    
    # Write the file to the output location
    molIndex = shapeList[0]
    finalMol = molList[molIndex]

    # Do a final match to make sure that all coordinates, including those of matching hydrogens, are identical
    MCSresult = rdFMCS.FindMCS([refMol, finalMol], 
                               ringMatchesRingOnly=False,  
                               matchValences=True, 
                               completeRingsOnly=False, 
                               timeout=timeoutTime, 
                               bondCompare=rdFMCS.BondCompare.CompareOrderExact, 
                               atomCompare=rdFMCS.AtomCompare.CompareElements,
                               maximizeBonds=False)
    complSMARTS = MCSresult.smartsString
    complMol = Chem.MolFromSmarts(complSMARTS)
    complConf = Chem.Conformer(complMol.GetNumAtoms())
    print("Found a MCSS with " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds!")

    refMatch = refMol.GetSubstructMatch(complMol)
    targetMatch = finalMol.GetSubstructMatch(complMol)
    atomMap = []
    for i in range(len(refMatch)): atomMap.append([targetMatch[i], refMatch[i]])
    AllChem.AlignMol(finalMol, refMol, atomMap = atomMap)    

    complIdxs = []
    for atom in complMol.GetAtoms(): complIdxs.append(atom.GetIdx())
    
    for refIdx, targetIdx, complIdx in zip(refMatch, targetMatch, complIdxs):   
          refCoords = refMol.GetConformer().GetAtomPosition(refIdx)
          finalMol.GetConformer().SetAtomPosition(targetIdx, refCoords)
    
    # Calculating tanimoto distance and formal charge
    fpRef = FingerprintMols.FingerprintMol(refMol)
    fpFinal = FingerprintMols.FingerprintMol(finalMol)
    taniSim = DataStructs.FingerprintSimilarity(fpRef, fpFinal)
    
    formalChargeFinal = rdmolops.GetFormalCharge(finalMol)
    formalChargeRef = rdmolops.GetFormalCharge(refMol)
    
    # Give additional information as output
    
    print("The final output molecule has a formal charge of " + str(formalChargeFinal) + ", and contains " + str(finalMol.GetNumAtoms()) +" atoms and " +  str(finalMol.GetNumBonds())+ " bonds.")
    print("The reference molecule has a formal charge of " + str(formalChargeRef) + ", and contains " + str(refMol.GetNumAtoms()) +" atoms and " +  str(refMol.GetNumBonds()) + " bonds.")
    print("The Tanimoto Similarity between the final structure and the reference compound is " + str(taniSim))
    print("The maximum common substructure (MCSS) contained " + str(complMol.GetNumAtoms()) + " atoms and " + str(complMol.GetNumBonds()) + " bonds")
    
    finalMol.SetProp("_Name", targetName)
    finalMol.SetProp("Tanimoto", str(taniSim))
    
    SDwriter.write(finalMol)

SDwriter.close()
print("Written file to " + finalOutPath)
sys.exit(0)

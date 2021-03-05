# FEPprep
Python-based alignment tool for free energy perturbation (FEP) preparation

### What the script does:
- After reading in the input information, the script tries to find the maximum common substructure (MCSS) between the reference compound and the input SMILES. It performs this search for 5 seconds as a default, but this can be configured using the -t flag. 
- Subsequently, the positions of the atoms that are part of this MCSS in the input SMILES are set to same coordinates from the matching atoms of the reference molecule. 
- Afterwards, a constrained embed is performed, generating different conformers. The MCSS core is kept constrained, while the atoms that are not part it can move freely. In default 50 conformers are generated, but this is configurable with the -c flag. 
- Next, the atoms that are part of the MCSS are sligthly readjusted to perfectly match the atoms of the reference molecule again. 
- Finally, the volume overlap between each conformer and the reference compound is calculated. From this, the one that has the biggest overlap with the reference molecule can be determined, and this result is subsequently written as a file to the specified output path.
### Command flags:
- -r or --ref: The filepath to the reference sdf or mol file
- -s or --smiles: The SMILES-string of a molecule you want allign with the reference molecule. This argument should be wrapped like a string i.e.: "NC(C(=O)N1Cc2c(C1)cccc2)CC(=O)N1CCN(CC1)C(c1ccc(cc1)F)c1ccc(cc1)F"
- -o or --out: The filepath that indicates where the outputfile should be written to
- -c or --conf: [optional] The amount of conformations that should be generated. The default value is 50.
- -t or --time: [optional] The amount of time spent searching for the maximum common substructure (MCSS). Usually, a result is quickly found, but the algoritm tries to exhaustively check every possibility. After the set time (in seconds) is up, the algorithm will continue with the best found result. The default value is 5
  

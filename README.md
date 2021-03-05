# FEPprep
Python-based alignment tool for free energy perturbation (FEP) preparation

### Command flags:
- -r or --ref: The filepath to the reference sdf or mol file
- -s or --smiles: The SMILES-string of a molecule you want allign with the reference molecule (this argument should be wrapped like a string like this: "NC(C(=O)N1Cc2c(C1)cccc2)CC(=O)N1CCN(CC1)C(c1ccc(cc1)F)c1ccc(cc1)F")
- -o or --out: The filepath that indicates where the outputfile should be written to
- -c or --conf: [optional] The amount of conformations that should be generated. The default value is 50.
- -t or --time: [optional] The amount of time spent searching for the maximum common substructure (MCSS). Usually, a result is quickly found, but the algoritm tries to exhaustively check every possibility. After the set time (in seconds) is up, the algorithm will continue with the best found result. The default value is 5
  

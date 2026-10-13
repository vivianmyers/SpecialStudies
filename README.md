# Pipeline
PDB mmCIF -> Add Hydrogens -> Find standard (covalent) bonds -> Find H-bonds -> Error filters -> Export

# UCSF ChimeraX1.12
- open mmCIF
- addh
- save file as PBD_H format mmCIF
- ```runscript "/full-path/export_covalent_bonds.py"```

# Citations
- UCSF ChimeraX: Tools for structure building and analysis. Meng EC, Goddard TD, Pettersen EF, Couch GS, Pearson ZJ, Morris JH, Ferrin TE. Protein Sci. 2023 Nov;32(11):e4792.
- H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne, The Protein Data Bank (2000) Nucleic Acids Research 28: 235-242 https://doi.org/10.1093/nar/28.1.235.

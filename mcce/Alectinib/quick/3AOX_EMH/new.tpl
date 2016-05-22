### This is a temporary parameter file made for residue EDO ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST EDO        EDOBK 

NATOM    EDOBK      4

IATOM    EDOBK  C1     0
IATOM    EDOBK  O1     1
IATOM    EDOBK  C2     2
IATOM    EDOBK  O2     3

ATOMNAME EDOBK    0  C1 
ATOMNAME EDOBK    1  O1 
ATOMNAME EDOBK    2  C2 
ATOMNAME EDOBK    3  O2 

CONNECT  EDOBK  C1  ion        0    O1   0    C2 
CONNECT  EDOBK  O1  ion        0    C1 
CONNECT  EDOBK  C2  ion        0    C1   0    O2 
CONNECT  EDOBK  O2  ion        0    C2 


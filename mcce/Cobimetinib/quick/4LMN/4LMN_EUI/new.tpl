### This is a temporary parameter file made for residue ANP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST ANP        ANPBK 

NATOM    ANPBK      31

IATOM    ANPBK  PG     0
IATOM    ANPBK  O1G    1
IATOM    ANPBK  O2G    2
IATOM    ANPBK  O3G    3
IATOM    ANPBK  PB     4
IATOM    ANPBK  O1B    5
IATOM    ANPBK  O2B    6
IATOM    ANPBK  N3B    7
IATOM    ANPBK  PA     8
IATOM    ANPBK  O1A    9
IATOM    ANPBK  O2A   10
IATOM    ANPBK  O3A   11
IATOM    ANPBK  O5'   12
IATOM    ANPBK  C5'   13
IATOM    ANPBK  C4'   14
IATOM    ANPBK  O4'   15
IATOM    ANPBK  C3'   16
IATOM    ANPBK  O3'   17
IATOM    ANPBK  C2'   18
IATOM    ANPBK  O2'   19
IATOM    ANPBK  C1'   20
IATOM    ANPBK  N9    21
IATOM    ANPBK  C8    22
IATOM    ANPBK  N7    23
IATOM    ANPBK  C5    24
IATOM    ANPBK  C6    25
IATOM    ANPBK  N6    26
IATOM    ANPBK  N1    27
IATOM    ANPBK  C2    28
IATOM    ANPBK  N3    29
IATOM    ANPBK  C4    30

ATOMNAME ANPBK    0  PG 
ATOMNAME ANPBK    1  O1G
ATOMNAME ANPBK    2  O2G
ATOMNAME ANPBK    3  O3G
ATOMNAME ANPBK    4  PB 
ATOMNAME ANPBK    5  O1B
ATOMNAME ANPBK    6  O2B
ATOMNAME ANPBK    7  N3B
ATOMNAME ANPBK    8  PA 
ATOMNAME ANPBK    9  O1A
ATOMNAME ANPBK   10  O2A
ATOMNAME ANPBK   11  O3A
ATOMNAME ANPBK   12  O5'
ATOMNAME ANPBK   13  C5'
ATOMNAME ANPBK   14  C4'
ATOMNAME ANPBK   15  O4'
ATOMNAME ANPBK   16  C3'
ATOMNAME ANPBK   17  O3'
ATOMNAME ANPBK   18  C2'
ATOMNAME ANPBK   19  O2'
ATOMNAME ANPBK   20  C1'
ATOMNAME ANPBK   21  N9 
ATOMNAME ANPBK   22  C8 
ATOMNAME ANPBK   23  N7 
ATOMNAME ANPBK   24  C5 
ATOMNAME ANPBK   25  C6 
ATOMNAME ANPBK   26  N6 
ATOMNAME ANPBK   27  N1 
ATOMNAME ANPBK   28  C2 
ATOMNAME ANPBK   29  N3 
ATOMNAME ANPBK   30  C4 

CONNECT  ANPBK  PG  ion        0    O1G  0    O2G  0    O3G  0    N3B
CONNECT  ANPBK  O1G ion        0    PG 
CONNECT  ANPBK  O2G ion        0    PG 
CONNECT  ANPBK  O3G ion        0    PG 
CONNECT  ANPBK  PB  ion        0    O1B  0    O2B  0    N3B  0    O3A
CONNECT  ANPBK  O1B ion        0    PB 
CONNECT  ANPBK  O2B ion        0    PB 
CONNECT  ANPBK  N3B ion        0    PG   0    PB 
CONNECT  ANPBK  PA  ion        0    O1A  0    O2A  0    O3A  0    O5'
CONNECT  ANPBK  O1A ion        0    PA 
CONNECT  ANPBK  O2A ion        0    PA 
CONNECT  ANPBK  O3A ion        0    PB   0    PA 
CONNECT  ANPBK  O5' ion        0    PA   0    C5'
CONNECT  ANPBK  C5' ion        0    O5'  0    C4'
CONNECT  ANPBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  ANPBK  O4' ion        0    C4'  0    C1'
CONNECT  ANPBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  ANPBK  O3' ion        0    C3'
CONNECT  ANPBK  C2' ion        0    C3'  0    O2'  0    C1'
CONNECT  ANPBK  O2' ion        0    C2'
CONNECT  ANPBK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  ANPBK  N9  ion        0    C1'  0    C8   0    N7   0    C4 
CONNECT  ANPBK  C8  ion        0    N9   0    N7   0    C5 
CONNECT  ANPBK  N7  ion        0    N9   0    C8   0    C5 
CONNECT  ANPBK  C5  ion        0    C8   0    N7   0    C6   0    C4 
CONNECT  ANPBK  C6  ion        0    C5   0    N6   0    N1 
CONNECT  ANPBK  N6  ion        0    C6 
CONNECT  ANPBK  N1  ion        0    C6   0    C2 
CONNECT  ANPBK  C2  ion        0    N1   0    N3 
CONNECT  ANPBK  N3  ion        0    C2   0    C4 
CONNECT  ANPBK  C4  ion        0    N9   0    C5   0    N3 

### This is a temporary parameter file made for residue MG_ ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MG_        MG_BK 

NATOM    MG_BK      1

IATOM    MG_BK  MG     0

ATOMNAME MG_BK    0  MG 

CONNECT  MG_BK  MG  ion      


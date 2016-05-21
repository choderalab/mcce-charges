### This is a temporary parameter file made for residue ACP ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST ACP        ACPBK 

NATOM    ACPBK      31

IATOM    ACPBK  PG     0
IATOM    ACPBK  O1G    1
IATOM    ACPBK  O2G    2
IATOM    ACPBK  O3G    3
IATOM    ACPBK  PB     4
IATOM    ACPBK  O1B    5
IATOM    ACPBK  O2B    6
IATOM    ACPBK  C3B    7
IATOM    ACPBK  PA     8
IATOM    ACPBK  O1A    9
IATOM    ACPBK  O2A   10
IATOM    ACPBK  O3A   11
IATOM    ACPBK  O5'   12
IATOM    ACPBK  C5'   13
IATOM    ACPBK  C4'   14
IATOM    ACPBK  O4'   15
IATOM    ACPBK  C3'   16
IATOM    ACPBK  O3'   17
IATOM    ACPBK  C2'   18
IATOM    ACPBK  O2'   19
IATOM    ACPBK  C1'   20
IATOM    ACPBK  N9    21
IATOM    ACPBK  C8    22
IATOM    ACPBK  N7    23
IATOM    ACPBK  C5    24
IATOM    ACPBK  C6    25
IATOM    ACPBK  N6    26
IATOM    ACPBK  N1    27
IATOM    ACPBK  C2    28
IATOM    ACPBK  N3    29
IATOM    ACPBK  C4    30

ATOMNAME ACPBK    0  PG 
ATOMNAME ACPBK    1  O1G
ATOMNAME ACPBK    2  O2G
ATOMNAME ACPBK    3  O3G
ATOMNAME ACPBK    4  PB 
ATOMNAME ACPBK    5  O1B
ATOMNAME ACPBK    6  O2B
ATOMNAME ACPBK    7  C3B
ATOMNAME ACPBK    8  PA 
ATOMNAME ACPBK    9  O1A
ATOMNAME ACPBK   10  O2A
ATOMNAME ACPBK   11  O3A
ATOMNAME ACPBK   12  O5'
ATOMNAME ACPBK   13  C5'
ATOMNAME ACPBK   14  C4'
ATOMNAME ACPBK   15  O4'
ATOMNAME ACPBK   16  C3'
ATOMNAME ACPBK   17  O3'
ATOMNAME ACPBK   18  C2'
ATOMNAME ACPBK   19  O2'
ATOMNAME ACPBK   20  C1'
ATOMNAME ACPBK   21  N9 
ATOMNAME ACPBK   22  C8 
ATOMNAME ACPBK   23  N7 
ATOMNAME ACPBK   24  C5 
ATOMNAME ACPBK   25  C6 
ATOMNAME ACPBK   26  N6 
ATOMNAME ACPBK   27  N1 
ATOMNAME ACPBK   28  C2 
ATOMNAME ACPBK   29  N3 
ATOMNAME ACPBK   30  C4 

CONNECT  ACPBK  PG  ion        0    O1G  0    O2G  0    O3G  0    C3B
CONNECT  ACPBK  O1G ion        0    PG 
CONNECT  ACPBK  O2G ion        0    PG 
CONNECT  ACPBK  O3G ion        0    PG 
CONNECT  ACPBK  PB  ion        0    O1B  0    O2B  0    C3B  0    O3A
CONNECT  ACPBK  O1B ion        0    PB 
CONNECT  ACPBK  O2B ion        0    PB 
CONNECT  ACPBK  C3B ion        0    PG   0    PB 
CONNECT  ACPBK  PA  ion        0    O1A  0    O2A  0    O3A  0    O5'
CONNECT  ACPBK  O1A ion        0    PA 
CONNECT  ACPBK  O2A ion        0    PA 
CONNECT  ACPBK  O3A ion        0    PB   0    PA 
CONNECT  ACPBK  O5' ion        0    PA   0    C5'
CONNECT  ACPBK  C5' ion        0    O5'  0    C4'
CONNECT  ACPBK  C4' ion        0    C5'  0    O4'  0    C3'
CONNECT  ACPBK  O4' ion        0    C4'  0    C1'
CONNECT  ACPBK  C3' ion        0    C4'  0    O3'  0    C2'
CONNECT  ACPBK  O3' ion        0    C3'
CONNECT  ACPBK  C2' ion        0    C3'  0    O2'  0    C1'
CONNECT  ACPBK  O2' ion        0    C2'
CONNECT  ACPBK  C1' ion        0    O4'  0    C2'  0    N9 
CONNECT  ACPBK  N9  ion        0    C1'  0    C8   0    C4 
CONNECT  ACPBK  C8  ion        0    N9   0    N7   0    C5   0    C4 
CONNECT  ACPBK  N7  ion        0    C8   0    C5 
CONNECT  ACPBK  C5  ion        0    C8   0    N7   0    C6   0    C4 
CONNECT  ACPBK  C6  ion        0    C5   0    N6   0    N1 
CONNECT  ACPBK  N6  ion        0    C6 
CONNECT  ACPBK  N1  ion        0    C6   0    C2 
CONNECT  ACPBK  C2  ion        0    N1   0    N3 
CONNECT  ACPBK  N3  ion        0    C2   0    C4 
CONNECT  ACPBK  C4  ion        0    N9   0    C8   0    C5   0    N3 

### This is a temporary parameter file made for residue MG_ ###
### Make sure that all the parameters are verified before using this file as a global parameter file ###

CONFLIST MG_        MG_BK 

NATOM    MG_BK      1

IATOM    MG_BK  MG     0

ATOMNAME MG_BK    0  MG 

CONNECT  MG_BK  MG  ion      


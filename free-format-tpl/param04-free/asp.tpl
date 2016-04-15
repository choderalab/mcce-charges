# MCCE free format tpl file
# Format rules:
# 1) One key value pair per line
# 2) key and value are separated by ":"
# 3) fields in kwy and value are separated by ","
# 4) quotation marks are optional
# 5) leading and ending spaces in a field must be quoted.

# Section 1. Atom charge and radius, connectivity
CONNECT, ASPBK, " N  " : sp2, " C  ", " CA ", " H  "
CONNECT, ASPBK, " H  " : s, " N  "
CONNECT, ASPBK, " CA " : sp3, " N  ", " C  ", " CB ", " HA "
CONNECT, ASPBK, " HA " : s, " CA "
CONNECT, ASPBK, " C  " : sp2, " CA ", " O  ", " N  "
CONNECT, ASPBK, " O  " : sp2, " C  "
CONNECT, ASP01, " CB " : sp3, " CA ", " CG ", "1HB ", "2HB "
CONNECT, ASP01, "1HB " : s, " CB "
CONNECT, ASP01, "2HB " : s, " CB "
CONNECT, ASP01, " CG " : sp2, " CB ", " OD1", " OD2"
CONNECT, ASP01, " OD1" : sp3, " CG ", " HD1"
CONNECT  ASP01  HD1 s         0     OD1
CONNECT  ASP01  OD2 sp2       0     CG
CONNECT  ASP02  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASP02 1HB  s         0     CB
CONNECT  ASP02 2HB  s         0     CB
CONNECT  ASP02  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASP02  OD1 sp2       0     CG
CONNECT  ASP02  OD2 sp3       0     CG  0     HD2
CONNECT  ASP02  HD2 s         0     OD2
CONNECT  ASP-1  CB  sp3       0     CA  0     CG  0    1HB  0    2HB
CONNECT  ASP-1 1HB  s         0     CB
CONNECT  ASP-1 2HB  s         0     CB
CONNECT  ASP-1  CG  sp2       0     CB  0     OD1 0     OD2
CONNECT  ASP-1  OD1 sp2       0     CG
CONNECT  ASP-1  OD2 sp2       0     CG

# Section 2, Basic Conformer Information: name, pka, em, rxn.
PKA, ASP01: 0.0
PKA, ASP02: 0.0
PKA, ASP-1: 4.75
EM       ASP01      0.0
EM       ASP02      0.0
EM       ASP-1      0.0
RXN      ASP01      -2.93
RXN      ASP02      -3.13
RXN      ASP-1      -20.2

# Section 3, Rotamers if possible

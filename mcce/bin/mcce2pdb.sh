#!/bin/bash
# Cat Chenal 2011-08-02
# To convert mcce file (eg step2) to pdb format 
# for proper display/selection in pymol

if [ $# -lt 1 ]; then
  echo "Usage: file_name required [saveas name]."
  exit 0
fi
file_in=$1
if [[ ! -f $file_in ]]; then
  echo $(basename $0)":: " $file_in": File not found."
  exit 0
fi
if [ $# -eq 2 ]; then
  saveas=$2
else
  saveas=${file_in%.*}".PDB"
fi

MEM=$(eval grep -q 'MEM' $1 &> /dev/null)
if [[ $MEM -eq 0 ]]; then     #found
  grep MEM $1 > mem.tmp
  sed -i '/MEM/d' $1
fi

#ATOM      1  N  aARG A  17     -21.893 -17.001  28.775  1.00  0.00           N              :: pdb
#123456789.123456789.123456789.
#1     2      3       4         5       6         7       8          9           10
#MODEL        1      
#ATOM  80913  OE1bGLU B0148_049 -20.029 -24.298  -5.781   1.400      -0.550      -1O000M000
#ATOM  80914  OE2bGLU B0148_049 -18.520 -22.720  -5.706   1.400      -0.550      -1O000M000
#ATOM  80942  CB cGLU B0148_053 -19.877 -22.567  -2.819   2.000       0.000      -1H005M000
#ENDMDL
#1     2      3  4    5         6       7        8        9           10         11

awk 'BEGIN{prt12="%-6s%5d %4s%s%3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n"}
     { if (NF<=2) { if (NF==1){print} else { printf ("%-13s%d\n", $1, $2) } }
     else
     {
        { if (NF=11) { if (length($3)<4) { atm=sprintf(" %-3s",$3) } else { atm=$3 };
                       { if (length($4)==4){ alt=substr($4,1,1); res=substr($4,2,3) }  else { alt=" "; res=$4 } } }
          else { len=length($3); atm=substr($3,1,len-4); alt=substr($3,len-4,1); res=substr($3,len-3,3);
                 if (length($3)<4) { atm=sprintf(" %-3s",$3) } else { atm=$3 } }
               }
         { if ($4=="_CL") {$4=" CL"; $11="Cl"} else { if ($4=="CYD") {$4="CYS"} } }
         { printf(prt12,$1,$2,atm,alt,res,substr($5,1,1),substr($5,2,4),$6,$7,$8,$9,$10) } } }' $file_in > $saveas

sed -i 's/\( [A|B]\)00/\1  /; s/\( [A|B]\)0/\1 /; s/_CL/ CL/; s/CYD/CYS/' $saveas

if [[ $MEM -eq 0 ]]; then
  sed -i '$r mem.tmp' $saveas
  /bin/rm *tmp
fi

 echo "Converted pdb was saved as "$saveas

#!/bin/bash
#  To get all conformers of given residue id (= nnnn) into pdb that
#  can be read by pymol or vmd (by removing the history column of step2_out.pdb)
#
if [[ $# -lt 4 ]]; then
  echo "Required: (1): residue identifier: nnnn; (2): titration point: as in fort.38 header; (3): occ_threshold%: nn; (4): run id"
  exit 0
fi
Res=$1
occ=$(printf "%6.2f" $(echo "scale=2; $3/100"|bc))
Run=$4

# get occ_col:
Col=$(head -1 fort.38 | awk -v col="$2" '{ for(n=2;n<=NF;n++) {if ( $n == col ) print n} }')
echo 'Retrieve confs of '$Res 'at titration point '$titra' (col '$Col'), with cutoff of '$occ

awk -v res="$Res" -v col="$Col" -v lim="$occ" '$1 ~ res"_" { if ($col >= lim) {print substr($1,6,9)} }' fort.38 > confs
echo "Total confs for cutoff: " $(wc -l confs)

skip_next=1
if [ $skip_next != 0 ]; then

for id in $(cat confs)
do
 #  if [ -f $Run"_"$id"_confs.PDB" ]; then
 #   /bin/rm $Run"_"$id"_confs.PDB"
 #  fi
#   touch $Run"_"$id"_confs.PDB"

   grepline="'"${id%_*}"_000|"$id"'"
   eval egrep $grepline step2_out.pdb | cut -c -74  > $Run"_"$id".PDB"
done
/bin/rm confs
fi

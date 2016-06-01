#!/bin/bash
# This function outputs a pdb file of conformers for the given residue id.
# It formats the output using the MODEL/ENDMDL keywords so that pymol will load them at once,
# each representing a different state. Thus their states can be run as a movie
# The second argument is the titration point in fort.38 at which the most occupied conformers will be retrieved.
# If the 3rd argument is 0; the pdb file format will be the same as step2, else it will be a standard pdb.
# If only 3 arguments are given, the conformers with an occupancy > 0.05 (default) are used;
# the default is overwritten when a fourth argument is given (can be 0 or any other cutoff).
# Example: get_conf_models.sh A0100 2 1
#
if [[ $# -lt 3 ]]; then
   echo "Required: ARG1:res_id (Znnn) ARG2:occ_col ARG3:convert_to_pdb (0/1) [ ARG4:occ_lim ] "
   echo $(basename $0) " must be called from a working directory (step2_out.pdb and fort.38 are required)"
   exit 0
fi
RES=$1
COL=$2
if [[ $COL -eq 1 ]]; then
   COL=2
   echo "Given column was 1; defaulted to 2"
fi
most_occ_file=$RES"_"$COL"_mostocc.csv"
res_confs_pdb=$RES"_"$COL"_S2.pdb"
if [[ -f $most_occ_file ]]; then
   /bin/rm $most_occ_file
fi
if [[ -f $res_confs_pdb ]]; then
   /bin/rm $res_confs_pdb
fi

occ_lim="0.05"
if [[ $# -eq 4 ]]; then
   occ_lim=$4
fi
if [[ $occ_lim == "0" ]]; then
   awk -v res="$RES" -v col="$COL" '{if ($1 ~ res) { print $1, $col } }' fort.38 |sort -k2nr > $most_occ_file
else
   awk -v res="$RES" -v col="$COL" -v lim="$occ_lim" '{if ($1 ~ res) { if ($col > lim){ print $1, $col}}}' fort.38 |sort -k2nr > $most_occ_file
fi
# GLU01A0148_008 0.272
#      678901234:9

touch $res_confs_pdb

conf_bkb=$RES"_000"

for conf in $(awk '{print substr($1,6,9)}' $most_occ_file)
do
   i=$(echo $conf| awk '{ num=substr($1,7,3)*1; printf ("%d", num) }' )
   echo "MODEL        "$i >> $res_confs_pdb
   egrep_line="'"$conf_bkb"|"$conf"'"
   eval egrep $egrep_line step2_out.pdb >> $res_confs_pdb
   echo "ENDMDL" >> $res_confs_pdb
done

if [[ $3 ]]; then

   mcce2pdb.sh $res_confs_pdb # saves the converted output with a ".PDB" extension
   /bin/rm $res_confs_pdb

fi
if [[ -f $res_confs_pdb ]]; then
   /bin/em $res_confs_pdb
fi


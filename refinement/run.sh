
$DIR_PATH ="/mnt/nas1/lanwei-125/FGF5/disulfide_peptide_cluster/"
$RESULT_PATH ="/mnt/nas1/lanwei-125/FGF5/disulfide_peptide_cluster/refinement/"
$recet_path ="/mnt/nas1/lanwei-125/FGF5/dock_prepare/FGF5.pdb"
$np_process_number ="12"
$rec_id ="A"
$rec_id ="E"

nohup python3 /home/weilan/script/refinement/rosetta_refinement.py -i $DIR_PATH -o $RESULT_PATH -rec $recet_path -rec_id $rec_id -lig_id $rec_id -n $np_process_number > /home/weilan/script/refinement/rosetta_refinement.log 2>&1 &
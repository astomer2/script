#!/bin/bash

$txt_path="/mnt/nas1/lanwei-125/CD44/Sequence/OPN.txt"
$config="/mnt/nas1/lanwei-125/CD44/docking/prepare_dock/CD44S.trg"
$modle_num=5
$num_steps=3500000
$core_num=12
$repeat_times=5
$work_path="/mnt/nas1/lanwei-125/CD44/docking/ADCP"
$log_path=$work_path/ADCP_log.txt"




python /mnt/sdc/lanwei/script-1/ADCP/new-adcp.py -t "$txt_path" -c "$config"  -m "$modle_num" -n "$num_steps" -p "$core_num" -r "$repeat_times" -w "$work_path" > "$log_path" 2>&1 

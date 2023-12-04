#!/bin/bash

# 定义参数
WORK_PATH="/mnt/nas1/lanwei-125/FGF5/FGF5-pos/new_pos/MD/pos/CYMAEYMYC"
log_path=$WORK_PATH"/log.txt"
CUDA_DEVICE_ID=3 
PARAMETER_FILE="/mnt/sdc/lanwei/script-1/MD/amber_paramter_files"
CMAP_PATH="/mnt/sdc/lanwei/script-1/MD/CMAP_files"
 


nohup python /mnt/sdc/lanwei/script-1/MD/run-simulation.py -i "$WORK_PATH" -g "$CUDA_DEVICE_ID" -CAMP 1 -cys 1 > "$log_path" 2>&1 & 
#!/bin/bash

# 定义参数
WORK_PATH="/mnt/nas1/lanwei-125/PRLR/MD/NEPLYHLVT/"
CUDA_DEVICE_ID=2
PARAMETER_FILE="/mnt/sdc/lanwei/script-1/MD/amber_paramter_files"
CMAP_PATH="/mnt/sdc/lanwei/script-1/MD/CMAP_files"
 
# 执行命令
nohup python /mnt/sdc/lanwei/script-1/MD/run-simulation.py -i "$WORK_PATH" -g "$CUDA_DEVICE_ID" -p "$PARAMETER_FILE" -c "$CMAP_PATH" &

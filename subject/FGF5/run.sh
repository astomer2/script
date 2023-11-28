#!/bin/bash

# 定义参数
WORK_PATH="/mnt/nas1/lanwei-125/FGF5/FGF5-pos/cyco/CTDKNCPLG/"
CUDA_DEVICE_ID=0
PARAMETER_FILE="/home/weilan/script/MD/amber_paramter_files"
CMAP_PATH="/home/weilan/script/MD/CMAP_files"
 
# 执行命令
nohup python /home/weilan/script/MD/run-simulation.py -i "$WORK_PATH" -g "$CUDA_DEVICE_ID" -p "$PARAMETER_FILE" -c "$CMAP_PATH" &

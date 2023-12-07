#!/bin/bash

# 定义参数
WORK_PATH="/mnt/nas1/lanwei-125/FGF5/Dis_to_cov/CC_to_SS/MD/SEMAEYMYS/"
log_path=$WORK_PATH"/log.txt"
CUDA_DEVICE_ID=1

nohup python /mnt/sdc/lanwei/script-1/MD/run-simulation.py -i "$WORK_PATH" -g "$CUDA_DEVICE_ID"  -cyc 1 > "$log_path" 2>&1 & 
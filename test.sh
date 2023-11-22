#!/bin/bash

# 获取当前文件夹中子文件夹数量
current_folder_count=$(ls -l |grep "^d" | wc -l)

# 暂存之前的子文件夹数量
previous_folder_count=$current_folder_count

while true; do
  sleep 1200 # 暂停30秒

  # 获取当前文件夹中子文件夹数量
  current_folder_count=$(ls -l |grep "^d" | wc -l)

  # 比较当前子文件夹数量与之前的数量是否相等
  if [ $current_folder_count -eq $previous_folder_count ]; then
    # 执行关闭并重新运行的命令
    pkill ADCP-mutidock.py
    python3 ADCP-mutidock.py

    # 更新之前的子文件夹数量为当前数量，以便下一次比较
    previous_folder_count=$current_folder_count
  else
    # 更新之前的子文件夹数量为当前数量，以便下一次比较
    previous_folder_count=$current_folder_count
  fi

done
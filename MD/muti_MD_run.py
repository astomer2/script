import os
import time
import tqdm
import subprocess
import logging
from pynvml.smi import nvidia_smi

# 常量


# 初始化日志
logging.basicConfig(level=logging.INFO)

def get_tasks(subject):

    tasks = [file for file in os.listdir(subject) if os.path.isdir(os.path.join(subject, file))]
    return tasks

def monitor_gpu_usage():
    nvsmi = nvidia_smi.getInstance()
    gpu_ids = list(range(nvsmi.DeviceQuery('count')['count']))
    gpu_utils = [gpu['utilization']['gpu_util'] for gpu in nvsmi.DeviceQuery('utilization.gpu')['gpu']]
    gpu_usage = list(zip(gpu_ids, gpu_utils))
    return gpu_usage

def instantly_available_gpu():

    enable_gpuid = [gpu_id for gpu_id, gpu_util in monitor_gpu_usage() if gpu_util < 5]
    return enable_gpuid

def determine_available_gpu():

    pre_gpu_ids = instantly_available_gpu()
    time.sleep(20)
    upgrade_gpu_ids = instantly_available_gpu()
    gpu_ids = list(set(pre_gpu_ids) & set(upgrade_gpu_ids))
    return gpu_ids

def run_command(script_path, work_path, gpu_id):

    command = [script_path, "-i", work_path, "-g", str(gpu_id), "-cyc", "1"]
    subprocess.Popen(command, stdout=open(work_path, 'w'), stderr=subprocess.STDOUT, preexec_fn=os.setpgrp)
    logging.info(f"运行命令：{command}")

def run_tasks(script_path, subject, tasks):

    gpu_ids = determine_available_gpu()
    for task, gpu_id in zip(tasks, gpu_ids):
        logging.info(f"运行任务：{task}，GPU：{gpu_id}")
        work_path = os.path.join(subject, task)
        run_command(script_path, work_path, gpu_id)

def show_progress(tasks):
    for _ in tqdm.tqdm(tasks, desc="运行任务中"):
        pass

def output_log(tasks):

    for task in tasks:
        log_path = f"{task}.log"
        try:
            with open(log_path, "r") as f:
                logging.info(f"{task} 的日志：")
                logging.info(f.read())
        except FileNotFoundError:
            logging.warning(f"找不到任务 {task} 的日志文件")

def main():
    tasks = get_tasks(SUBJECT_DIR)

    # 监控 CUDA 核心使用率
    monitor_gpu_usage()

    # 并行运行任务
    run_tasks(SCRIPT_PATH, SUBJECT_DIR, tasks)

    # 使用 tqdm 显示进度
    show_progress(tasks)

    # 输出每个任务的日志
    output_log(tasks)

if __name__ == "__main__":
    SUBJECT_DIR = "subject"
    SCRIPT_PATH = "/mnt/sdc/lanwei/script-1/MD/muti_MD_run.py"
    main()

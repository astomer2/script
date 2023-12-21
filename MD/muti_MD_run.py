import os
import time
from tqdm import tqdm
import subprocess
import collections
import logging
from pynvml.smi import nvidia_smi
from pathlib import Path
# 常量


# 初始化日志
logging.basicConfig(level=logging.INFO)

def get_tasks(subject):
    subject_path = Path(subject)
    tasks = [task for task in subject_path.iterdir() if task.is_dir()]
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
    logging.info(f"运行命令：{command}")

def output_log(tasks):

    for task in tasks:
        log_path = f"{task}.log"
        try:
            with open(log_path, "r") as f:
                logging.info(f"{task} 的日志：")
                logging.info(f.read())
        except FileNotFoundError:
            logging.warning(f"找不到任务 {task} 的日志文件")

def main(subject, script_path):
    gpu_ids = collections.deque(determine_available_gpu()) 
    logging.info(gpu_ids)
    ran_task = []
    tasks_list = get_tasks(subject)

    while len(ran_task) < len(tasks_list):
        # 使用 tqdm 创建进度条，总任务数为 len(tasks)
        tasks_queue = collections.deque(tasks_list)
        for _ in tqdm(range(len(tasks_list))):
            if not gpu_ids:
                logging.info("没有可用的 GPU，等待 5 秒后重试...")
                time.sleep(5)
                gpu_ids = collections.deque(determine_available_gpu())

            if gpu_ids:
                gpu_id = gpu_ids.popleft()
                work_path = tasks_queue.popleft()
                run_command(script_path, work_path, gpu_id)
                ran_task.append(work_path)

if __name__ == "__main__":
    SUBJECT_DIR = "subject"
    SCRIPT_PATH = "/mnt/sdc/lanwei/script-1/MD/muti_MD_run.py"
    main(subject, script_path)

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
    """
    Get a list of tasks from the specified subject directory.

    :param subject: The directory containing the tasks.
    :type subject: str
    :return: A list of tasks.
    :rtype: list
    """
    tasks = [file for file in os.listdir(subject) if os.path.isdir(os.path.join(subject, file))]
    return tasks

def monitor_gpu_usage():
    """
    Generate a function comment for the given function body in a markdown code block with the correct language syntax.

    :return: A string representing the function comment in Markdown format.
    """
    nvsmi = nvidia_smi.getInstance()
    gpu_ids = list(range(nvsmi.DeviceQuery('count')['count']))
    gpu_utils = [gpu['utilization']['gpu_util'] for gpu in nvsmi.DeviceQuery('utilization.gpu')['gpu']]
    gpu_usage = list(zip(gpu_ids, gpu_utils))
    return gpu_usage

def instantly_available_gpu():
    """
    Retrieves the IDs of all available GPUs with GPU utilization less than 5%.

    Returns:
        List[int]: A list of GPU IDs that have utilization less than 5%.
    """
    enable_gpuid = [gpu_id for gpu_id, gpu_util in monitor_gpu_usage() if gpu_util < 5]
    return enable_gpuid

def determine_available_gpu():
    """
    Determines the available GPUs by checking the instantly available GPUs before and after a 20-second delay.
    
    Returns:
        List[int]: A list of available GPU IDs.
    """
    pre_gpu_ids = instantly_available_gpu()
    time.sleep(20)
    upgrade_gpu_ids = instantly_available_gpu()
    gpu_ids = list(set(pre_gpu_ids) & set(upgrade_gpu_ids))
    return gpu_ids

def run_command(script_path, work_path, gpu_id):
    """
    Run a command using the specified script path, work path, and GPU ID.

    Args:
        script_path (str): The path to the script to be executed.
        work_path (str): The working directory for the command.
        gpu_id (int): The ID of the GPU to be used.

    Returns:
        None

    Raises:
        None
    """
    command = [script_path, "-i", work_path, "-g", str(gpu_id), "-cyc", "1"]
    subprocess.Popen(command, stdout=open(work_path, 'w'), stderr=subprocess.STDOUT, preexec_fn=os.setpgrp)
    logging.info(f"运行命令：{command}")

def run_tasks(script_path, subject, tasks):
    """
    Run multiple tasks using the given script path, subject, and task names.

    Args:
        script_path (str): The path to the script to be run.
        subject (str): The subject of the tasks.
        tasks (list): A list of task names.

    Returns:
        None
    """
    gpu_ids = determine_available_gpu()
    for task, gpu_id in zip(tasks, gpu_ids):
        logging.info(f"运行任务：{task}，GPU：{gpu_id}")
        work_path = os.path.join(subject, task)
        run_command(script_path, work_path, gpu_id)

def show_progress(tasks):
    for _ in tqdm.tqdm(tasks, desc="运行任务中"):
        pass

def output_log(tasks):
    """
    Outputs the log for each task in the given list of tasks.

    Parameters:
    - tasks (list): A list of tasks for which the logs need to be output.

    Returns:
    - None

    Description:
    - For each task in the list of tasks:
        - Constructs the path to the log file using the task name.
        - Tries to open the log file in read mode.
        - If the log file is found:
            - Logs an informational message indicating the task name.
            - Logs the contents of the log file.
        - If the log file is not found:
            - Logs a warning message indicating that the log file for the task could not be found.
    """
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

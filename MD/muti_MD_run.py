import os
import time
#from click import File
from tqdm import tqdm
import subprocess
import collections
import logging
from pynvml.smi import nvidia_smi
from pathlib import Path
import argparse
# 常量
now_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# 初始化日志
logging.basicConfig(level=logging.INFO)

def get_tasks(subject):
    """
    Retrieves a list of tasks within the specified subject.

    Parameters:
        subject (str): The path to the subject directory.

    Returns:
        List[Path]: A list of tasks (directories) within the subject directory.
    """
    subject_path = Path(subject)
    tasks = [task for task in subject_path.iterdir() if task.is_dir()]
    return tasks


def monitor_gpu_usage():
    """
    Monitors the GPU usage and returns a list of tuples containing the GPU ID and its corresponding utilization percentage.

    Returns:
        List[Tuple[int, float]]: A list of tuples, where each tuple contains the GPU ID (int) and its utilization percentage (float).
    """
    nvsmi = nvidia_smi.getInstance()
    gpu_ids = [item['uuid'] for item in nvsmi.DeviceQuery('uuid')['gpu']]
    gpu_utils = [gpu['utilization']['gpu_util'] for gpu in nvsmi.DeviceQuery('utilization.gpu')['gpu']]
    gpu_usage = list(zip(gpu_ids, gpu_utils))
    return gpu_usage

def instantly_available_gpu():

    enable_gpuid = [gpu_id for gpu_id, gpu_util in monitor_gpu_usage() if gpu_util < 5]
    return enable_gpuid

def determine_available_gpu():
    """
    Determines the available GPU by calling the `instantly_available_gpu` function four times with a 15-second delay between each call. It then finds the intersection of the GPU IDs returned by each call and returns the resulting list of available GPU IDs.

    Returns:
        List[int]: A list of available GPU IDs.

    """

    gpu_ids_0 = instantly_available_gpu()
    time.sleep(15)
    gpu_ids_1 = instantly_available_gpu()
    time.sleep(15)
    gpu_ids_2 = instantly_available_gpu()
    time.sleep(15)
    gpu_ids_3 = instantly_available_gpu()

    gpu_ids = list(set(gpu_ids_0) & set(gpu_ids_1) & set(gpu_ids_2) & set(gpu_ids_3))
    return gpu_ids

def run_command(script_path, work_path, extra_functions, gpu_id):
    """
    Runs a command using the given script path, work path, extra functions, and GPU ID.

    Args:
        script_path (str): The path to the script to be executed.
        work_path (str): The working directory where the command will be executed.
        extra_functions (list): A list of additional functions to be included in the command.
        gpu_id (str): The ID of the GPU to be used for execution.

    Returns:
        None
    """
    extra_functions = ' '.join(['-' + x for x in extra_functions])
    command = ("nohup python3 %s -i %s -g %s %s > %s.log 2>&1 &" % (script_path, work_path, gpu_id, extra_functions, work_path ))
    os.system(command)
    logging.info(f"{now_time} 运行命令：{command}")

def output_log(tasks):

    for task in tasks:
        log_path = f"{task}.log"
        try:
            with open(log_path, "r") as f:
                logging.info(f"{task} 的日志：")
                logging.info(f.read())
        except FileNotFoundError:
            logging.warning(f"找不到任务 {task} 的日志文件")

def run_muti_MD(subject: str, script_path: str,extra_functions: list):  
    """
    This function is the main function of the program. It takes in three parameters:
    
    - subject: The subject of the program.
    - script_path: The path to the script.
    - extra_functions: A list of extra functions.
    
    It initializes a deque object `gpu_ids` with the available GPUs. It then logs the `gpu_ids`.
    
    It initializes an empty list `ran_task` and gets the list of tasks using the `get_tasks` function.
    
    It enters a while loop that continues until the number of `ran_task` is equal to the number of tasks.
    
    Inside the while loop, it creates a list `tasks` that contains tasks that are not in `ran_task`.
    It then initializes a deque object `tasks_queue` with the tasks.
    
    If `gpu_ids` is empty, it logs a message indicating that there are no available GPUs and waits for 45 seconds before retrying.
    
    If `gpu_ids` is not empty, it pops the leftmost GPU id from `gpu_ids` and the leftmost task from `tasks_queue`.
    
    It checks if a file named "mdinfo" exists in the `work_path`. If it does, it appends `work_path` to `ran_task` and logs a message indicating that the task has already been run.
    
    If the file does not exist, it calls the `run_command` function with the `script_path`, `work_path`, `extra_functions`, and `gpu_id` as arguments. It then appends `work_path` to `ran_task`.
    
    After each iteration of the while loop, it calculates the progress as the ratio of `ran_task` to the total number of tasks and writes the progress to the console using the tqdm library.
    """
    gpu_ids = collections.deque(determine_available_gpu())
    logging.info(gpu_ids)
    ran_task = []
    tasks_list = get_tasks(subject)

    # 创建一个字典，用于存储已完成的任务
    finished_tasks = {}

    while len(ran_task) < len(tasks_list):
        # 使用 tqdm 创建进度条，总任务数为 len(tasks)
        tasks = [task for task in tasks_list if task not in ran_task]
        tasks_queue = collections.deque(tasks)

        # 检查是否有可用的 GPU
        if not gpu_ids:
            # 如果没有可用的 GPU，则等待 45 秒
            logging.info(f"{now_time} 没有可用的 GPU，等待 45 秒后重试...")
            time.sleep(45)
            gpu_ids = collections.deque(determine_available_gpu())

        # 如果有可用的 GPU，则分配一个 GPU 并运行任务
        if gpu_ids:
            gpu_id = gpu_ids.popleft()
            work_path = tasks_queue.popleft()

            # 检查任务是否已完成
            if work_path in finished_tasks:
                ran_task.append(work_path)
                logging.info(f"{now_time} {work_path} 已经运行，跳过")
                continue

            # 运行任务
            run_command(script_path, work_path, extra_functions, gpu_id)
            ran_task.append(work_path)
            finished_tasks[work_path] = True

        # 使用 tqdm 创建进度条，总任务数为 len(tasks)，len(ran_task) 为已经运行的任务数，使用这两个值作为tqdm进度
        progress = len(ran_task) / len(tasks_list)
        tqdm.write(f"{now_time} 进度：{progress:.2%}")


if __name__ == "__main__":
    default_script_path = '/mnt/nas1/lanwei-125/script/MD/run_simulation.py'
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--subject_dir", type=str, help="The directory of the subject")
    parser.add_argument("-s","--script_path", type=str,default=default_script_path, help=f"The default path is {default_script_path} ")
    parser.add_argument("-ex","--extra_functions", type=str, nargs = "+",default="", help="The keyword of the extra functions,more detials plesase see the python run_simulation.py -h script")
    
    args = parser.parse_args()
    subject_dir = args.subject_dir
    script_path = args.script_path
    extra_functions = args.extra_functions
    run_muti_MD(subject_dir, script_path, extra_functions)

import os
import time
import tqdm
from pynvml.smi import nvidia_smi
from pynvml.smi import 
def get_tasks(subject):

    tasks = []
    for file in os.listdir(subject):
        if os.path.isdir(os.path.join(subject, file)):
            tasks.append(file)
    return tasks

def monitor_gpu_usage():

    nvsmi = nvidia_smi.getInstance()
    gpu_ids = list(range(nvsmi.DeviceQuery('count')['count']))
    gpu_utils = [gpu['utilization']['gpu_util'] for gpu in nvsmi.DeviceQuery('utilization.gpu')['gpu']]
    gpu_usage = list(zip(gpu_ids, gpu_utils))

    return gpu_usage

def Instantly_available_GPU():

    enable_gpuid = []
    for gpu_id, gpu_util in monitor_gpu_usage():
        if gpu_util < 5:
            enable_gpuid.append(gpu_id)
    return enable_gpuid

def Determine_available_gpu():
    pre_gpu_ids = Instantly_available_GPU()
    time.sleep(20)
    Upgrade_gpu_ids = Instantly_available_GPU()
    gpu_ids = list(set(pre_gpu_ids) & set(Upgrade_gpu_ids))
    return gpu_ids

def run_command(script_path, work_path, gpu_id):
    
    os.system(f"nohup {script_path} -i {work_path} -g {gpu_id} -cyc 1 > {work_path} 2>&1 &")

def run_tasks(script_path, subject, tasks):
    
    gpu_ids = Determine_available_gpu()
    for task, gpu_id in zip(tasks, gpu_ids):
        print(f"正在运行任务：{task}，GPU：{gpu_id}")
        work_path = os.path.join(subject,task)
        run_command(script_path,work_path, gpu_id)

def show_progress(tasks):

    for task in tqdm.tqdm(tasks):
        pass

def output_log(tasks):

    for task in tasks:
        log_path = f"{task}.log"
        with open(log_path, "r") as f:
            print(f"{task}的log：")
            print(f.read())

def main():
    subject = "subject"
    tasks = get_tasks(subject)

    script_path = "/mnt/sdc/lanwei/script-1/MD/muti_MD_run.py"
    # 监控CUDA核心使用率
    monitor_gpu_usage()

    # 并行运行任务
    run_tasks(script_path, subject, tasks)

    # 添加tqdm进度条显示task进度
    show_progress(tasks)

    # 输出每个task的log
    output_log(tasks)

if __name__ == "__main__":
    main()
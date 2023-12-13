#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import argparse
import multiprocessing
from tqdm import tqdm  # Import tqdm for a progress bar

class ABian:
    def __init__(self, txt_path, config_file_path, model_num, num_steps, core_num, repeat_times, work_path, cyclic, cystein):
        print("*************程序开始执行*************")
        print("--start time:%s--\n" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
        self.txt_path = txt_path
        self.config_file_path = config_file_path
        self.model_num = model_num
        self.num_steps = num_steps
        self.core_num = core_num
        self.repeat_times = repeat_times
        self.work_path = work_path
        self.cyclic = cyclic
        self.cystein = cystein
    # 生成命令
    def create_executive(self, key):
        if self.cyclic:
            self.executive_txt = f"adcp -t {self.config_file_path} -s {key} -N {self.model_num} -n {self.num_steps} -c {self.core_num} -cyc -o dock > log.txt"
        elif self.cystein:
            self.executive_txt = f"adcp -t {self.config_file_path} -s {key} -N {self.model_num} -n {self.num_steps} -c {self.core_num} -cys -o dock > log.txt"
        else:
            self.executive_txt = f"adcp -t {self.config_file_path} -s {key} -N {self.model_num} -n {self.num_steps} -c {self.core_num} -o dock > log.txt"

    # 创建文件夹
    def make_dir(self, dir_name):
        max_files = self.repeat_times + 1
        self.make_dir_name = ["%s-%d" % (dir_name, i) for i in range(1, max_files)]
        self.make_dir_path = [os.path.join(self.work_path, item) for item in self.make_dir_name]
        
        for item in self.make_dir_path:
            if not os.path.exists(item):
                os.umask(0)
                os.makedirs(item)
            else:
                print("%s-->目录不为空，请检查！！！" % item)

    # 执行命令
    def executive(self, item):
        os.chdir(item)
        if os.path.exists(os.path.join(self.work_path, item, ".executed")):
            pass
        else:
            command = self.executive_txt
            subprocess.getstatusoutput(command)
            os.chdir(self.work_path)
            open(os.path.join(item, ".executed"), "w").close()


    def execute_in_parallel(self, items):
        with multiprocessing.Pool(processes=self.repeat_times) as pool:
            pool.map(self.executive, items)

    def run(self):
        
        if not os.path.exists(self.work_path):
            os.umask(0)
            os.makedirs(self.work_path)

        with open(self.txt_path, encoding='utf-8') as r_file:
            task_count = sum(1 for _ in r_file)  # Count the number of lines in the file

        with tqdm(total=task_count, desc="Processing", ncols=100) as pbar:
            with open(self.txt_path, encoding='utf-8') as r_file:

                for line in r_file.readlines():
                    file_and_key_name = line.strip('\n')  # Remove the newline character
                    self.make_dir(file_and_key_name)                                               
                    self.create_executive(file_and_key_name)
                    self.execute_in_parallel(self.make_dir_path)

                    for _ in self.make_dir_path:
                        pbar.update(1)
                    self.make_dir_name.clear()
                    self.make_dir_path.clear()


        print("\n***************程序运行完毕！！！***************")
        print("--结束时间：%s--\n" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

def main():
    parser = argparse.ArgumentParser(description="Run amber simulation with specified parameters")
    parser.add_argument("-t", "--txt_path", required=True, help="Path to the txt file")
    parser.add_argument("-c", "--config_file_path", required=True, help="Path to the config file")
    parser.add_argument("-m", "--model_num", type=int, required=True, help="Model number")
    parser.add_argument("-n", "--num_steps", type=int, required=True, help="Number of steps")
    parser.add_argument("-p", "--core_num", type=int, required=True, help="Number of cores")
    parser.add_argument("-r", "--repeat_times", type=int, required=True, help="Repeat times")
    parser.add_argument("-w", "--work_path", required=True, help="Path to the work directory")
    parser.add_argument(
            '-cyc', "--cyclic", dest="cyclic", action="store_true",
            default=False,
            help="option for cyclic peptide through backbone")
    parser.add_argument(
            '-cys', "--cystein", dest="cystein", action="store_true",
            default=False,
            help="option for cyclic peptide through CYS-S-S-CYS")


    args = parser.parse_args()

    abian = ABian(args.txt_path, args.config_file_path, args.model_num, args.num_steps, args.core_num, args.repeat_times, args.work_path, args.cyclic, args.cystein)
    abian.run()

if __name__ == '__main__':
    main()
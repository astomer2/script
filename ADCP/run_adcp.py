#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import argparse
import multiprocessing
from tqdm import tqdm 
import logging
import time
import pandas as pd
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
now_times = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

class ADCP:
    def __init__(self, input_file, config_file_path, output_model_num, run_step_num, core_num, repeat_times, work_path, cyclic, cystein):

        self.input_file = input_file
        self.config_file_path = config_file_path
        self.output_model_num = output_model_num
        self.run_step_num = run_step_num
        self.core_num = core_num
        self.repeat_times = repeat_times
        self.work_path = work_path
        self.cyclic = cyclic
        self.cystein = cystein

    def create_executive(self, key):
        if self.cyclic:
            self.executive_command = f"adcp -t {self.config_file_path} -s {key} -N {self.output_model_num} -n {self.run_step_num} -c {self.core_num} -cyc -o dock > log.txt"
        elif self.cystein:
            self.executive_command = f"adcp -t {self.config_file_path} -s {key} -N {self.output_model_num} -n {self.run_step_num} -c {self.core_num} -cys -o dock > log.txt"
        else:
            self.executive_command = f"adcp -t {self.config_file_path} -s {key} -N {self.output_model_num} -n {self.run_step_num} -c {self.core_num} -o dock > log.txt"

    def make_dir(self, dir_name):
        max_files = self.repeat_times + 1
        self.make_dir_name = ["%s-%d" % (dir_name, i) for i in range(1, max_files)]
        self.make_dir_path = [os.path.join(self.work_path, item) for item in self.make_dir_name]
        
        for item in self.make_dir_path:
            if not os.path.exists(item):
                os.umask(0)
                os.makedirs(item)
            else:
                logger.info("%s already exists" % item)

    def executive(self, item):
        os.chdir(item)
        if os.path.exists(os.path.join(self.work_path, item, ".executed")):
            pass
        else:
            command = self.executive_command
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

        with open(self.input_file, encoding='utf-8') as f:
            if self.input_file.endswith(".txt"):
                sequences = f.readlines()
            elif self.input_file.endswith(".csv"):
                df = pd.read_csv(self.input_file, encoding='utf-8')
                if 'Sequence' in df.columns:
                    sequences = df['Sequence'].tolist()
                else:
                    logger.error("The csv file does not contain the 'Sequence' column.")
                    exit()
            else:
                logger.error("The input file must be a txt or csv file.")
                exit()

        completed_tasks = 0

        for seq in sequences:
            file_and_key_name = seq.strip('\n')  # Remove the newline character
            self.make_dir(file_and_key_name)
            self.create_executive(file_and_key_name)
            self.execute_in_parallel(self.make_dir_path)
            completed_tasks += 1
            completion_ratio = completed_tasks / len(sequences)
            logger.info(f"Progress: {completion_ratio:.2%}")

            self.make_dir_name.clear()
            self.make_dir_path.clear()


def main():
    parser = argparse.ArgumentParser(description="Run amber simulation with specified parameters")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the txt file")
    parser.add_argument("-t", "--config_file_path", required=True, help="Path to the config file")
    parser.add_argument("-m", "--output_model_num", type=int, required=True, help="Model number")
    parser.add_argument("-n", "--run_step_num", type=int, required=True, help="Number of steps")
    parser.add_argument("-c", "--core_num", type=int, required=True, help="Number of cores")
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

    adcp = ADCP(args.input_file, args.config_file_path, args.output_model_num, args.run_step_num, args.core_num, args.repeat_times, args.work_path, args.cyclic, args.cystein)
    adcp.run()

if __name__ == '__main__':
    main()
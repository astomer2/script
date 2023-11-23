#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import time
import xml.etree.ElementTree as ET
import multiprocessing


class ABian:
    def __init__(self):
        print("*************程序开始执行*************")
        print("--start time:%s--\n" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
        self.i = 0
        self.executive_txt = None
        tree = ET.parse("/mnt/nas1/lanwei-125/IL8/v3/ADCP/monomer/application.xml")
        root = tree.getroot()
        for pathName in root.findall('pathName'):
            self.txt_path = pathName[0].text
            self.config_file_path = pathName[1].text
            self.work_path = pathName[2].text
        for param in root.findall('param'):
            self.model_num = int(param[0].text)
            self.num_steps = int(param[1].text)
            self.core_num = int(param[2].text )
            self.repeat_times = int(param[3].text)
        # 文件夹名称
        self.make_dir_name = []

        # 文件夹路径
        self.make_dir_path = []

    # 生成命令
    def create_executive(self, key):
        self.executive_txt = f"adcp -t {self.config_file_path} -s {key} -N {self.model_num} -n {self.num_steps} -c {self.core_num} -o dock > log.txt"

    # 生成名称
    def create_file_name(self, file_name):
        max_files = self.repeat_times + 1
        for i in range(1, max_files):
            self.make_dir_name.append("%s-%d" % (file_name, i))

    # 创建文件夹
    def make_dir(self, dir_name):
        self.create_file_name(dir_name)
        # 创建文件夹
        for item in self.make_dir_name:
            child_dir_name = os.path.join(self.work_path, item)
            if not os.path.exists(child_dir_name):
                os.makedirs(child_dir_name)
                self.make_dir_path.append(child_dir_name)
            else:
                print("%s-->目录不为空，请检查！！！" % child_dir_name)

    # 执行命令
    def executive(self):
            # 创建进程池，设置并行数量为4
        pool = multiprocessing.Pool(processes=self.repeat_times)
        for item in self.make_dir_path:
            # 判断是否已经执行过命令
            if not os.path.exists(os.path.join(item, ".executed")):
                # 生成命令
                command = self.executive_txt
                # 将任务加入进程池
                pool.apply_async(run_command, (command, item))
        # 关闭进程池
        pool.close()
        pool.join()


        if self.i % 20 == 0:
            print("\n------------------执行过程日志------------------")
            print("\n--执行时间：%s--" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))
            print("执行命令---%d----完毕\n" % self.i)

        self.i += 1

    def run(self):
        # 读取txt文件
        with open(self.txt_path, encoding='utf-8') as r_file:
            for line in r_file.readlines():
                file_and_key_name = line.strip('\n')  # 去掉列表中每一个元素的换行符
                if not os.path.exists(self.work_path):
                    os.makedirs(self.work_path)
                # 创建文件夹
                self.make_dir(file_and_key_name)
                # 在文件夹下面执行命令
                self.create_executive(file_and_key_name)
                self.executive()
                self.make_dir_name.clear()
                self.make_dir_path.clear()
                # 删除标记文件
                for item in self.make_dir_path:
                    if os.path.exists(os.path.join(item, ".executed")):
                        os.remove(os.path.join(item, ".executed"))

        print("\n***************程序运行完毕！！！***************")
        print("--结束时间：%s--\n" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))


# 记录文件夹中是否运行过程序
def run_command(command, path):
    # 进入子文件夹
    os.chdir(path)
    # 执行命令
    subprocess.getstatusoutput(command)
    # 创建标记文件
    open(".executed", "w").close()


if __name__ == '__main__':
    abian = ABian()
    abian.run()


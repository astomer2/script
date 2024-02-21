import re
from re import S
import time
from tqdm import tqdm
import pandas as pd
from logging import Logger
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.by import By
from datetime import datetime

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager



def read_sequences(file_path, cys):
    sequences = []

    # 读取文件内容
    if file_path.endswith('.txt'):
        with open(file_path, 'r') as f:
            lines = f.readlines()
    elif file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
        lines = df['Sequence']

    # 处理每一行序列
    for line in lines:
        seqs = line.strip()
        if cys:
            cys_pos = [substr.start() for substr in re.finditer('C', seqs)]
            if len(cys_pos) == 2:
                seqs = f"{seqs} -ss {cys_pos[0]+1} {cys_pos[1]+1}"
        sequences.append(seqs)

    seq_groups = [sequences[i:i+10] for i in range(0, len(sequences), 10)]
    return seq_groups

def submit_sequence(browser, seq_list):
    fasta_seq = browser.find_element(By.XPATH, '//*[@id="fastaseq2"]')
    fasta_seq.send_keys('\n'.join(seq_list))
    time.sleep(1.5)

def submit_receptor(browser):
    pdb_file = browser.find_element(By.XPATH,'//*[@id="pdbfile1"]')
    pdb_file.send_keys(Receptor_file)
    time.sleep(1.5)


def submit_optional(browser):
    sub_optional = browser.find_element(By.XPATH,'//*[@id="option1"]')
    sub_optional.click()
    
'''def submit_reference(browser):
    ref_pdb_file = browser.find_element(By.XPATH,'//*[@id="optionForm1"]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/UL[1]/LI[1]/INPUT[1]')
    ref_pdb_file.send_keys(reference_file)
    time.sleep(1.5)'''

    
def submit_task(browser):
    submit_button = browser.find_element(By.XPATH,'//*[@id="form1"]/table/tbody/tr[8]/td[1]/input')
    submit_button.click()
    time.sleep(1.5)

def save_result(browser, log_file):
    url = browser.current_url
    current_time = datetime.now().strftime("%H:%M:%S")
    
    with open(log_file, 'a') as f:
        f.write(url + ' ' + current_time + '\n')

def run(seq_groups):
    url = 'http://huanglab.phys.hust.edu.cn/hpepdock/'
    #browser = webdriver.Chrome()
    driver_path = ChromeDriverManager().install()
    service = Service(driver_path)
    # 创建 Chrome WebDriver，并指定服务
    browser = webdriver.Chrome(service=service)

    for i, seq_list in enumerate(seq_groups):
        
        browser.get(url)
        submit_receptor(browser)
        submit_sequence(browser, seq_list)
        submit_optional(browser)
        #submit_reference(browser)
        submit_task(browser)
        save_result(browser, log_file)
    browser.close()    
    browser.quit()
        
if __name__ == '__main__':
    cys = True
    Receptor_file = r"Y:\lanwei-125\MC5R\Sequence\MC5R.pdb"
    seq_file = r"Y:\lanwei-125\MC5R\dock\HPEP\cpep_sequence\Sequence.csv"
    #reference_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\Afamelanotide.pdb"
    log_file = r"Y:\lanwei-125\MC5R\dock\HPEP\cpep_sequence\log.txt"
    seq_groups = read_sequences(seq_file, cys)
    run(seq_groups)
        

    
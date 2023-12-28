from re import S
import time
from tqdm import tqdm
import argparse
from DrissionPage import ChromiumOptions
from DrissionPage import ChromiumPage
from DrissionPage.common import By
from logging import Logger
from pathlib import Path

def read_sequences(file_path):
    with open(file_path, 'r') as f:
        seqs = [line.strip() for line in f]

    seq_groups = [seqs[i:i+10] for i in range(0, len(seqs), 10)]
    return seq_groups

def submit_sequence(page, seq_list):
    fasta_seq = page.ele((By.XPATH, '//*[@id="fastaseq2"]'))
    fasta_seq.input('\n'.join(seq_list))
    time.sleep(1.5)
'''
def submit_receptor(page, Receptor_file):
    pdb_file = page.ele((By.XPATH,'//*[@id="pdbfile1"]'))
    pdb_file.input(Receptor_file)
    print("Receptor file uploaded")
    time.sleep(1.5)
'''

def submit_receptor(page, Receptor_file):
    page.set.upload_files(Receptor_file)
    page.ele((By.XPATH,'//*[@id="pdbfile1"]')).click()
    page.wait.upload_paths_inputted()
    time.sleep(1.5)


def submit_optional(page):
    sub_optional = page.ele((By.XPATH,'//*[@id="option1"]'))
    sub_optional.click()

def submit_reference(page, reference_file):
    ref_pdb_file = page.ele((By.XPATH,'//*[@id="optionForm1"]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/UL[1]/LI[1]/INPUT[1]'))
    ref_pdb_file.input(reference_file)
    time.sleep(1.5)

def submit_task(page):
    submit_button = page.ele((By.XPATH,'//*[@id="form1"]/table/tbody/tr[8]/td[1]/input'))
    submit_button.click()
    time.sleep(1.5)

def save_result(page, log_file):
    url = page.url
    current_time = time.strftime("%H:%M:%S")

    with open(log_file, 'a') as f:
        f.write(url + ' ' + current_time + '\n')

def run(seq_file, Receptor_file, log_file, reference_file):
    seq_groups = read_sequences(seq_file)
    url = 'http://huanglab.phys.hust.edu.cn/hpepdock/'
    driver = ChromiumOptions().set_paths(browser_path='/snap/bin/chromium')
    driver.set_argument('--headless')
    page = ChromiumPage(addr_driver_opts=driver)
    for i, seq_list in enumerate(tqdm(seq_groups)):
        page.get(url)
        submit_receptor(page, Receptor_file)
        submit_sequence(page, seq_list)
        submit_optional(page)
        if reference_file is not None:
            submit_reference(page, reference_file)
        else:
            continue
        submit_task(page)
        save_result(page, log_file)
    page.quit()
'''
if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--seq', type=str, help='sequence file')
    parser.add_argument('-r', '--receptor', type=str, help='receptor file')
    parser.add_argument('-l', '--log', type=str, help='log file')
    parser.add_argument('-ref', '--reference', type=str, help='reference file')
    args = parser.parse_args()
    seq_file = args.seq
    Receptor_file = args.receptor
    log_file = args.log
    reference_file = args.reference
    
    run(seq_file, Receptor_file, log_file, reference_file)
'''

if __name__ == '__main__':
    seq_file = "/mnt/nas1/lanwei-125/PRLR/GA-generator/HPEP/test.txt"
    Receptor_file = "/mnt/nas1/lanwei-125/PRLR/structure_prepare/PRLR.pdb"
    log_file =  "/mnt/nas1/lanwei-125/PRLR/GA-generator/HPEP/log.txt"
    reference_file = ""
    run(seq_file, Receptor_file, log_file, reference_file)
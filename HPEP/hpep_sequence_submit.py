from re import S
import time
from tqdm import tqdm
from DrissionPage import ChromiumOptions
from DrissionPage import ChromiumPage
from DrissionPage.common import By

def read_sequences(file_path):
    with open(file_path, 'r') as f:
        seqs = [line.strip() for line in f]

    seq_groups = [seqs[i:i+10] for i in range(0, len(seqs), 10)]
    return seq_groups

def submit_sequence(page, seq_list):
    fasta_seq = page.ele((By.XPATH, '//*[@id="fastaseq2"]'))
    fasta_seq.input('\n'.join(seq_list))
    time.sleep(1.5)

def submit_receptor(page):
    pdb_file = page.ele((By.XPATH,'//*[@id="pdbfile1"]'))
    pdb_file.input(Receptor_file)
    time.sleep(1.5)

def submit_optional(page):
    sub_optional = page.ele((By.XPATH,'//*[@id="option1"]'))
    sub_optional.click()

'''def submit_reference(page, reference_file):
    ref_pdb_file = page.ele((By.XPATH,'//*[@id="optionForm1"]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/UL[1]/LI[1]/INPUT[1]'))
    ref_pdb_file.input(reference_file)
    time.sleep(1.5)'''

def submit_task(page):
    submit_button = page.ele((By.XPATH,'//*[@id="form1"]/table/tbody/tr[8]/td[1]/input'))
    submit_button.click()
    time.sleep(1.5)

def save_result(page, log_file):
    url = page.url
    current_time = time.strftime("%H:%M:%S")

    with open(log_file, 'a') as f:
        f.write(url + ' ' + current_time + '\n')

def run(seq_groups):
    url = 'http://huanglab.phys.hust.edu.cn/hpepdock/'
    driver = ChromiumOptions().set_paths(browser_path='/snap/bin/chromium')
    driver.set_argument('--headless')
    page = ChromiumPage(addr_driver_opts=driver)
    for i, seq_list in enumerate(tqdm(seq_groups)):
        page.get(url)
        submit_receptor(page)
        submit_sequence(page, seq_list)
        submit_optional(page)
        # submit_reference(page, reference_file)
        submit_task(page)
        save_result(page, log_file)

    page.quit()

if __name__ == '__main__':
    Receptor_file = "/mnt/nas1/lanwei-125/CD44/docking/prepare_dock/CD44S-H.pdb"
    seq_file = "/mnt/nas1/lanwei-125/CD44/Sequence/OPN.txt"
    #reference_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\Afamelanotide.pdb"
    log_file = "/mnt/nas1/lanwei-125/CD44/docking/HPEP/OPN-log.txt"
    seq_groups = read_sequences(seq_file)
    run(seq_groups)

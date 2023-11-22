import requests
from lxml import etree
from tqdm import tqdm
import time
import pandas as pd
import os

def get_estimated_solubility(sequence, proxies=None):
    data = {
    'hideInputFields': 'no' ,
    'nTerm': '(NH2-)' ,
    'sequence': sequence ,
    'cTerm': '(-COOH)' ,
    'aaCode': '0' ,
    'disulphideBonds': '' 
    }

    resp = requests.post(url='https://pepcalc.com/ppc.php' , data=data , proxies=proxies)

    root = etree.HTML(resp.text)
    root.xpath('//div[@class="neat"]/div/text()')

    tr_list = root.xpath('//div[@class="neat"]/div/table//tr')
    for tr_i in tr_list[1:7]:
        if 'Estimated solubility' in tr_i.xpath('./td[1]/text()')[0] :
            return tr_i.xpath('./td[2]//text()')[0]

if __name__ == '__main__':

    seq_file_path = r"D:\jobwork\subject\TRPV1\to_channel_hole\dock_v3\describe.txt"
    proxies = { "http": "http://192.168.1.254:7890", "https": "http://192.168.1.254:7890", }


    with open(seq_file_path) as f:
        seq_list = f.readlines()

    seq_estimated_solubility_list = []
    for seq_i in tqdm(seq_list):
        seq_dict = {}
        seq_dict['Sequence'] = seq_i.strip()
        # seq_dict['estimated_solubility'] = get_estimated_solubility(seq_i.strip() , proxies)
        seq_dict['estimated_solubility'] = get_estimated_solubility(seq_i.strip())
        seq_estimated_solubility_list.append(seq_dict)
        time.sleep(0.1)

    df = pd.DataFrame(seq_estimated_solubility_list)

    dir_name = os.path.dirname(seq_file_path)  
    file_name = os.path.basename(seq_file_path)
    new_file_name = 'estimated_solubility_' + file_name
    new_seq_file_path = os.path.join(dir_name, new_file_name)

    df.to_csv(new_seq_file_path , index=False)
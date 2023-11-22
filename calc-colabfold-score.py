import re
from pprint import pprint
import os
import pandas as pd
import numpy as np

def calc_score(log_file, score_file):
    with open(log_file) as f:
        lines = f.readlines()

    scores = {} 

    for i, line in enumerate(lines):
        if 'Query ' in line:
            name = re.search(r'Query.*: (.*) \(length', line).group(1)
            #print(name)
        if 'rank_001' in line:
            parts = line.split()
            pLDDT = float(parts[-3].split('=')[-1])  
            pTM = float(parts[-2].split('=')[-1])
            ipTM = float(parts[-1].split('=')[-1])
            weighted_score = round(ipTM*0.8 + pTM*0.2, 2)         
            scores[name] = {'pLDDT': pLDDT, 'pTM': pTM, 'ipTM': ipTM, 'score': weighted_score}
            #pprint(scores)

    if score_file.endswith ('.txt'):
        with open(score_file,'a+') as f:
            f.write('Sequence\tscore\tpLDDT\tpTM\tipTM\n')
            for name, data in scores.items():
                weighted_score = data['score']
                pLDDT = data['pLDDT']
                pTM = data['pTM']
                ipTM = data['ipTM']

                f.write(f'{name}\t{weighted_score}\t{pLDDT}\t{pTM}\t{ipTM}\n')
    elif score_file.endswith ('.csv'):
        df = pd.DataFrame(scores)
        df.to_csv(score_file, index=False)
    elif score_file.endswith ('.xlsx'):
        df = pd.DataFrame(scores)
        df.to_excel(score_file, index=True)
    else:
        raise Exception('Unsupported file type')
    print('yadaze!')


if __name__ == '__main__':
    log_file = '/mnt/sdc/lanwei/ADIPOR/output/log.txt'
    score_file = '/mnt/sdc/lanwei/ADIPOR/output/colab-score.xlsx'
    calc_score(log_file, score_file)
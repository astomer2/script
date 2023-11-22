# 氨基酸序号重排
from pprint import pprint

pdb = r"C:\Users\123\Desktop\jobwork\subject\IL-8\IL8-R.pdb"

with open(pdb, 'r') as f:
    lines = f.readlines()

all_lines = []


prev_res_id = None
res_num = 0

all_lines = []
for line in lines:
    if  line.startswith('ATOM'):
        res_id = line[22:27].strip()
        if prev_res_id != res_id: 
            res_num += 1
        line = line[:22] + str(res_num).rjust(4) + line[26:]
        prev_res_id = res_id
        all_lines.append(line)
    elif line.startswith('TER'):
        all_lines.append(line)
    
    elif line.startswith('END'):
        all_lines.append(line)

pprint(all_lines)

with open(pdb, 'w') as f:
    f.writelines(all_lines)
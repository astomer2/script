import os

input_file = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v3\hpep"
result_dir = r'C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v3\hpep\hpep-result'

os.makedirs(result_dir, exist_ok=True)

peptide_scores = {}

for root, dirs, files in os.walk(input_file):
    if 'hpepdock_all.out' in files:
        with open(os.path.join(root, 'hpepdock_all.out')) as f:
            for line in f:
                data = line.split() 
                peptide = data[4]
                score = float(data[3])
                Cluster = data[0]
                nmuber = data[1]
                if peptide not in peptide_scores:
                    peptide_scores[peptide] = []
                peptide_scores[peptide].append((Cluster, nmuber, score))

    if 'hpepdock_all.pdb' in files:
        with open(os.path.join(root, 'hpepdock_all.pdb')) as f:
            lines = []
            peptide = ''
            for line in f:
                if line.startswith('REMARK Cluster:'):
                    cluster, nmuber = line.split()[2:]
                    if nmuber == '1':
                        peptide = os.path.basename(root)
                        lines.append(line)
                elif line.startswith('END'):
                    break
                else:
                    lines.append(line)
        
        if lines:
            with open(os.path.join(result_dir, f'{peptide}.pdb'), 'w') as f:
                f.write(''.join(lines))
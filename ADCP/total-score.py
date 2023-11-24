import os
import numpy as np

def read_score():
    score_file = open(path + '/score.txt', 'w')
    score_file.write('name\tscore\n')
    for root, dirs, files in os.walk(path):
        for dir in dirs:
            sub_path = os.path.join(root, dir)
            if 'log.txt' not in os.listdir(sub_path):
                continue
            log_file = open(sub_path + '/log.txt', 'r')
            if sum(1 for line in log_file) < 25:
                print(f'your docking job {dir} was failed,please cheak or rerun')
                continue
            log_file.seek(0)
            lines = log_file.readlines()
            for i in range(14, 25):
                if lines[i].split()[0] == '1':
                    score = lines[i].split()[1]
                    break
            score_line = dir + '\t' + score + '\n'
            score_file.write(score_line)
    score_file.close()

def write_result():
    with open(path + '/score.txt', "r") as f:
        lines = f.readlines()
    name_score_dict = {}
    for line in lines:
        data = line.strip().split()
        if len(data) != 2:
            continue
        name, score = data
        try:
            score = float(score)
        except ValueError:
            continue
        name = name.split("-")[0]
        if name not in name_score_dict:
            name_score_dict[name] = []
        name_score_dict[name].append(score)

    with open(path + "/result.txt", "w") as f:
        f.write("sequence\tmix\tmax\tavg\tmed\tvar\n")
        for name, scores in name_score_dict.items():
            min =  np.min(scores)
            max = np.max(scores)
            med = np.median(scores)
            avg = np.mean(scores)
            var = np.var(scores)
    #        d = b/a + 2.5 * (variance_b/a)**0.5
            f.write(f"{name.upper()}\t{min}\t{max}\t{avg:.2f}\t{med:.2f}\t{var:.2f}\n")

if __name__ == '__main__':
    #path = os.getcwd()
    path = '/mnt/nas1/lanwei-125/IL8/v4/ADCP/monomer/'
    read_score()
    write_result()
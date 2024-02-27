import os



def find_trajectory_flie(dir_path):
    for root, dirs, files in os.walk(dir_path):
        if 'complex_md1.mdcrd' and 'mdinfo' in files:
            mdinfo_path = os.path.join(root, 'mdinfo')
            rajectory_path = os.path.join(root, 'complex_md.mdcrd')
            if os.path.exists(rajectory_path) and os.path.exists(mdinfo_path):
                with open(mdinfo_path) as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'Total steps' in line:
                            percentage_str = line.split('(')[1].split('%')[0]
                            percentage = float(percentage_str)
                            if percentage > 99.0:
                                os.remove(rajectory_path)    

def check_trajectory_file(dir_path: str)->tuple:
    MD = 0
    MD_path = []
    for root, dirs, files in os.walk(dir_path):
        for file in files:
            if file=='mdinfo':
                MD += 1
                md_path = os.path.join(root, file)
                with open(md_path) as f:
                    for line in f:
                        if line.startswith('|'):
                            try:
                                if line.split()[1] == 'Total':
                                    times = int(line.split()[3])/500000
                                    a = (root ,times)
                                    MD_path.append(a)
                            except IndexError:
                                pass
    return MD, MD_path

if __name__ == '__main__':
    dir_path = '/mnt/nas1/lanwei-125/MC5R/MD/pos/'
    find_trajectory_flie(dir_path)
    MD , MD_path = check_trajectory_file(dir_path)
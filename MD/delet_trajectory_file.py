import os 


def find_trajectory_flie(dir_path):
    for root, dirs, files in os.walk(dir_path):
        if 'complex_md1.mdcrd' and 'mdinfo' in files:
            mdinfo_path = os.path.join(root, 'mdinfo')
            rajectory_path = os.path.join(root, 'complex_md1.mdcrd')
            with open(mdinfo_path) as f:
                lines = f.readlines()
                for line in lines:
                    if 'Total steps' in line:
                        percentage_str = line.split('(')[1].split('%')[0]
                        percentage = float(percentage_str)
                        if percentage > 99.0:
                            os.remove(rajectory_path)    

if __name__ == '__main__':
    dir_path = '/mnt/nas1/lanwei-125/IL8/v3/'
    find_trajectory_flie(dir_path)
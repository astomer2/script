import os 

workpath = '/mnt/nas1/lanwei-125/FGF5/FGF5-pos/cyco/'

for root, dirs, files_names in os.walk(workpath):
    # 这里只关注文件夹路径
    subdirs_paths = [os.path.join(root, subdir) for subdir in dirs]
    
    # 获取子文件夹中的文件路径
    files_paths = [os.path.join(subdir_path, file_name) for subdir_path in subdirs_paths for _, _, files in os.walk(subdir_path) for file_name in files]
    
    # 打印文件路径
    for file_path in files_paths:
        print(file_path)


for root, dirs, files_names in os.walk(workpath):
    # 这里只关注文件夹路径
    for subdir in dirs:
        subdir_path = os.path.join(root, subdir)
        
        # 接下来可以对子文件夹进行操作，比如获取其中的文件路径
        for subroot, subdirs, subfiles in os.walk(subdir_path):
            for file_name in subfiles:
                file_path = os.path.join(subroot, file_name)
                print(file_path)
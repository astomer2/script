import pandas as pd

# 读取 Excel 文件
filename = r"D:\jobwork\subject\TRPV1\天然多肽设计抑制TRPV1.xlsx"
xl = pd.ExcelFile(filename)

# 遍历每个工作表
for sheet_name in xl.sheet_names:
    # 读取工作表
    df = pd.read_excel(filename, sheet_name=sheet_name, header=0)

    # 遍历每一列
    for i, col_name in enumerate(df.columns):
        # 生成 CSV 文件名
        csv_name = f"{sheet_name}_{i+1}.csv"

        # 将当前列的数据保存为 CSV 文件
        sub_df = df.loc[:, col_name].dropna()
        sub_df.to_csv(csv_name, index=False, header=['Sequence'])

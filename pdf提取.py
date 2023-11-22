'''import pdfplumber
from pprint import pprint
with pdfplumber.open(r"D:\jobwork\subject\soluble\C6185QJTG0-batchExceptPlateMap-cd64c3\C6185QJTG0-1-Peptide Qualitative Solubility Test Report.pdf") as pdf:
    page01 = pdf.pages[0] #指定页码
    #table1 = page01.extract_table()#提取单个表格
    table2 = page01.extract_tables()#提取多个表格
    
    pprint(table2)'''

import os
import pdfplumber
import pandas as pd
from pprint import pprint
# 定义函数，提取表格并转换为 DataFrame
def extract_table_df(table):
    # 转换为 DataFrame
    df = pd.DataFrame(table[1:], columns=table[0])
    # 去除空值
    df.dropna(how='all', inplace=True)
    # 填充 NaN
    df.fillna('', inplace=True)
    return df

# 定义函数，处理 PDF 文件并提取表格
def process_pdf(file_path):
    tables = []
    with pdfplumber.open(file_path) as pdf:
        for page in pdf.pages:
            for table in page.extract_table():
                df = extract_table_df(table)
                tables.append(df)
    return tables

# 定义函数，处理文件夹中的所有 PDF 文件
def process_folder(input_folder, output_file):
    all_tables = []
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".pdf"):
            file_path = os.path.join(input_folder, file_name)
            tables = process_pdf(file_path)
            all_tables += tables
    all_tables_df = pd.concat(all_tables)
    all_tables_df.to_excel(output_file, index=False)
    #pprint(all_tables_df)

# 使用示例
input_folder = r"D:\jobwork\subject\soluble\C6185QJTG0-batchExceptPlateMap-cd64c3"
output_file = r"D:\jobwork\subject\soluble\output.xlsx"
process_folder(input_folder, output_file)
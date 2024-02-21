from pathlib import Path
import os 
import pandas as pd
import matplotlib.pyplot as plt


contact = Path("/mnt/nas1/lanwei-125/MC5R/MD/pos/ZMK2-110/contact_frac_byatom.dat")
df1= pd.read_csv(contact, sep="\s+", skiprows=2)

contact_res = df1['Contact']
S = contact_res.str.split(':', expand=True)
df1['res_id']=(S[1].str.split('@', expand=True))[0]
result = (df1.groupby('res_id')['Nframes'].sum()).sort_values(ascending=False)
I = pd.DataFrame(result)
I.columns = ['Nframes']
I['resid'] = I.index
I['resid'] = I['resid'].astype(int)
I = I.sort_values('resid', ascending=True)
I = I.reset_index(drop=True)

all_res_id = pd.Series(range(1, 273), name='resid')
merged_result = pd.merge(all_res_id, I, left_on='resid', right_on='resid', how='left').fillna(0)
plt.plot(merged_result['resid'], merged_result['Nframes'])
plt.xlabel('MC5R_res_id')
plt.ylabel('Contact_num')
plt.title('Contact_res')
plt.grid(True)
for i, row in I.iterrows():
    if row['Nframes'] > 10000:
        plt.annotate(f"{int(row['resid'])}", (row['resid'], row['Nframes']), xytext=(5, 5), textcoords='offset points', ha='left', arrowprops=dict(arrowstyle='->'))
plt.gcf().set_size_inches(20, 8)
plt.tight_layout() 
# plt.show()
plt.savefig(contact.parent/'contact_res.png', dpi=300)
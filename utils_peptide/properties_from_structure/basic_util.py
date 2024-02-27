import shutil
from pathlib import Path
from icecream import ic
from pandas import DataFrame

ic.configureOutput(includeContext=True, argToStringFunction=str)
ic.lineWrapWidth = 120


residues_dict = {
    "GLY": "G",
    "ALA": "A",
    "SER": "S",
    "THR": "T",
    "CYS": "C",
    "VAL": "V",
    "LEU": "L",
    "ILE": "I",
    "MET": "M",
    "PRO": "P",
    "PHE": "F",
    "TYR": "Y",
    "TRP": "W",
    "ASP": "D",
    "GLU": "E",
    "ASN": "N",
    "GLN": "Q",
    "LYS": "K",
    "ARG": "R",
    "HIS": "H",
    "ACE": "",
    "NH2": "",
}


def calc_disulfide_pep(pdb_files, out_file, pdb_out_dir, pep_squences=None):
    """Returns:
    df with 2 columns: pep_squences, is_disulfide, with is_disulfide always True.
    """
    pdb_out_dir = Path(pdb_out_dir)
    if pep_squences is None:
        pep_squences = [file.name for file in pdb_files]
    assert len(pep_squences) == len(pdb_files)
    # 遍历所有pdb文件
    results = []
    for i, file in enumerate(pdb_files):
        # 打开pdb文件
        with open(file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        # print(file)
        # 找到SG原子所在行
        sg_atoms = []
        for line in lines:
            # HETATM
            if line.startswith("ATOM") and line[12:16].strip() == "SG":
                sg_atoms.append(line)
        # ic(file.name, sg_atoms)
        if len(sg_atoms) == 2:
            # 提取SG坐标
            x1, y1, z1 = [
                float(sg_atoms[0][30:38]),
                float(sg_atoms[0][38:46]),
                float(sg_atoms[0][46:54]),
            ]
            x2, y2, z2 = [
                float(sg_atoms[1][30:38]),
                float(sg_atoms[1][38:46]),
                float(sg_atoms[1][46:54]),
            ]
            # 计算距离
            dist = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5

            if dist < 2.05:
                is_disulfide = True
            else:
                is_disulfide = False
            results.append((pep_squences[i], is_disulfide))
            if is_disulfide:
                shutil.copy(file, pdb_out_dir / f"{pep_squences[i]}.pdb")
    df = DataFrame(results, columns=["pep_squences", "is_disulfide"])
    disulfide_df = df[df["is_disulfide"]]
    disulfide_count = len(disulfide_df)
    disulfide_ratio = disulfide_count / len(df)
    ic(disulfide_count, len(df), disulfide_ratio)
    disulfide_df.to_csv(out_file, index=False, sep=",")
    return disulfide_df


if __name__ == "__main__":
    pass

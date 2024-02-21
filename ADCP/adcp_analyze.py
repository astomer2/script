from dataclasses import dataclass
from pathlib import Path
import os
import numpy as np
import pandas as pd

@dataclass
class DockingResult:
    name: str
    score: float

def read_score(log_path: Path):
    docking_results = []
    score_file_path = log_path / 'score.txt'

    with score_file_path.open('w') as score_file:
        score_file.write('name\tscore\n')

        for subdir in log_path.iterdir():
            if subdir.is_dir() and 'log.txt' in os.listdir(subdir):
                log_file_path = subdir / 'log.txt'

                if sum(1 for _ in log_file_path.open()) < 25:
                    print(f'Your docking job {subdir.name} failed, please check or rerun')
                    continue

                with log_file_path.open('r') as log_file:
                    lines = log_file.readlines()

                    for i in range(14, 25):
                        if lines[i].split()[0] == '1':
                            score = float(lines[i].split()[1])
                            docking_results.append(DockingResult(subdir.name, score))
                            break

        for docking_result in docking_results:
            score_file.write(f'{docking_result.name}\t{docking_result.score}\n')
    return docking_results

def write_result(result_path: Path, docking_results: list):
    name_score_dict = {}

    for docking_result in docking_results:
        name = docking_result.name.split("-")[0]
        if name not in name_score_dict:
            name_score_dict[name] = []
        name_score_dict[name].append(docking_result.score)

    result_df = pd.DataFrame(columns=["sequence", "mix", "max", "avg", "med", "var", "machine_score"])
    result_data = []
    for name, scores in name_score_dict.items():
        min_score = np.min(scores)
        max_score = np.max(scores)
        med_score = np.median(scores)
        avg_score = np.mean(scores)
        var_score = np.var(scores)
        machine_score = avg_score + 2.5*(float(var_score)/float(len(scores)))**0.5

        result_data.append({
            "Sequence": name.upper(),
            "mix": np.around(min_score,3),
            "max": np.around(max_score,3),
            "avg": np.around(avg_score,3),
            "med": np.around(med_score,3),
            "var": np.around(var_score,3),
            "machine_score": np.around(machine_score,3)
        })
        result_df = pd.DataFrame(result_data)

    result_df.to_csv(result_path / "result.csv", index=False)

if __name__ == '__main__':
    path = Path('/mnt/nas1/lanwei-125/MC5R/dock/ADCP/CPEP_motif_4/AGRP/')    
    docking_results = read_score(path)
    write_result(path, docking_results)

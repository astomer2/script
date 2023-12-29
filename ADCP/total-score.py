from dataclasses import dataclass
from pathlib import Path
import os
import numpy as np

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

    result_file_path = result_path / "result.txt"

    with result_file_path.open("w") as result_file:
        result_file.write("sequence\tmix\tmax\tavg\tmed\tvar\tmachine_score\n")

        for name, scores in name_score_dict.items():
            min_score = np.min(scores)
            max_score = np.max(scores)
            med_score = np.median(scores)
            avg_score = np.mean(scores)
            var_score = np.var(scores)
            machine_score = avg_score + 2.5*(float(var_score)/float(len(scores)))**0.5

            result_file.write(f"{name.upper()}\t{min_score}\t{max_score}\t{avg_score:.2f}\t{med_score:.2f}\t{var_score:.2f}\t{machine_score:.2f}\n")

if __name__ == '__main__':
    path = Path('/mnt/nas1/lanwei-125/FGF5/GA_generator/')    
    docking_results = read_score(path)
    write_result(path, docking_results)

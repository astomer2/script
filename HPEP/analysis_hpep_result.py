import os
import statistics
import csv


def analyze_scores():
    peptide_scores = {}
    result_path = os.path.join(input_file, "results.csv")
    for root, dirs, files in os.walk(input_file):
        if "hpepdock_all.out" in files:
            with open(os.path.join(root, "hpepdock_all.out")) as f:
                for line in f:
                    data = line.split()
                    peptide = data[4]
                    score = float(data[3])

                    if peptide not in peptide_scores:
                        peptide_scores[peptide] = []

                    peptide_scores[peptide].append(score)
  
        with open(result_path, "w+", newline='') as csvfile:
            fieldnames = ['sequence', 'mix', 'max', 'avg', 'med', 'var', 'machine_score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for peptide, scores in peptide_scores.items():
                min_score = min(scores)
                max_score = max(scores)
                avg_score = round(statistics.mean(scores), 3)
                median_score = round(statistics.median(scores), 3)
                variance = round(statistics.variance(scores), 3)
                machine_score = round(avg_score + 2.5 * (float(variance) / float(len(scores))) ** 0.5,3)

                writer.writerow({
                    'sequence': peptide,
                    'mix': min_score,
                    'max': max_score,
                    'avg': avg_score,
                    'med': median_score,
                    'var': variance,
                    'machine_score': machine_score
                })

if __name__ == "__main__":

    input_file = "/mnt/nas1/lanwei-125/test/"

    analyze_scores()
    print("输出成功")
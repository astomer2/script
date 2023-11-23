import os
import statistics

input_file = "/mnt/nas1/lanwei-125/PRLR/HPEP/"
result_path = "/mnt/nas1/lanwei-125/PRLR/HPEP/PRL-hpep-result.txt"


def analyze_scores():
    peptide_scores = {}

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

        with open(result_path, "w+") as f:
            f.write("sequence\tmix\tmax\tavg\tmed\tvar\n")

            for peptide, scores in peptide_scores.items():
                min_score = min(scores)
                max_score = max(scores)
                avg_score = round(statistics.mean(scores), 3)
                median_score = round(statistics.median(scores), 3)
                variance = round(statistics.variance(scores), 3)

                f.write(
                    f"{peptide}\t{min_score}\t{max_score}\t{avg_score}"
                    f"\t{median_score}\t{variance}\n"
                )


analyze_scores()
print("输出成功")
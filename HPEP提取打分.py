import os
import statistics

input_file = r"C:\Users\123\Downloads\655c772b0a311"
result_path = r'C:\Users\123\Downloads\655c772b0a311\hpep-result.txt'

def analyze_scores():
    

    peptide_scores = {}

    for root, dirs, files in os.walk(input_file):
        

        if 'hpepdock_all.out' in files:

            with open(os.path.join(root, 'hpepdock_all.out')) as f:
        
                for line in f:
          
                    data = line.split()  
                    peptide = data[4]
                    score = float(data[3])
          
                    if peptide not in peptide_scores:
                      peptide_scores[peptide] = []
            
                    peptide_scores[peptide].append(score)
          
        with open(result_path, 'w+') as f:
            f.write("sequence\tmix\tmax\tavg\tmed\tvar\n")
    
            for peptide, scores in peptide_scores.items():
    
                min_score = min(scores)
                max_score = max(scores)  
                avg_score = statistics.mean(scores)
                median_score = statistics.median(scores)
                variance = statistics.variance(scores)
      
                f.write(f'{peptide}\t{min_score}\t{max_score}\t{avg_score}'
                        f'\t{median_score}\t{variance}\n')
              
analyze_scores()
print('输出成功')
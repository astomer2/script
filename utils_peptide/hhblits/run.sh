# 
file=peptide_utils/hhblits/neff.py
nohup python $file \
    --a3mdir /mnt/nas/alphafold/af_out/tasks/golden_tests/Sanner2023paper/3r7g_A_B/msas/B \
    > $file.log 2>&1 &
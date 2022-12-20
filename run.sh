python3 build_rf_matrix.py  hapA_lst.txt hapB_lst.txt  dataset/simulated_data/trio_data/long-reads/kidONT.R94.30X.p1.fasta.gz p1
python3 build_rf_matrix.py  hapA_lst.txt hapB_lst.txt  dataset/simulated_data/trio_data/long-reads/kidONT.R94.30X2.p2.fasta.gz p2
cat p1 p2 > kidONT.R94.kmers.txt

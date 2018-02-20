tar -xzf genomes.tar.gz
../BRAG -t Neurospora.nwk -r Nc -x Sm -o best_Neuro -g alignments/Neurospora.*.os.txt -q alignments/Neurospora.*.seqs -C Neurospora-crassa_v12_fixed_centromeres.txt
python extra_tracks.py
../BRAG -t Neurospora.nwk -r Nc -x Sm -o best_Neuro -g alignments/Neurospora.*.os.txt -q alignments/Neurospora.*.seqs -C Neurospora-crassa_v12_fixed_centromeres.txt -T Ncra_extra_tracks.tsv



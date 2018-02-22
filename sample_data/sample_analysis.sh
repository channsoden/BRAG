# Download genomes from Open Science Framework
wget -O genomes.tar.gz https://osf.io/vmkc8/?action=download
tar -xzf genomes.tar.gz

PYTHONPATH=$PYTHONPATH:../hannsoden-bioinformatics:..

# Generate extra tracks of data
python2 extra_tracks.py

# Run BRAG
python2 ../BRAG -t Neurospora.nwk -r Nc -x Sm -o best_Neuro -g alignments/Neurospora.*.os.txt -q alignments/Neurospora.*.seqs -C Neurospora-crassa_v12_fixed_centromeres.txt -T Ncra_extra_tracks.tsv

# Generate figures from BRAG output
python2 paper_figures.py

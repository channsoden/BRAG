# Download genomes from Open Science Framework
#wget -O genomes.tar.gz https://osf.io/vmkc8/?action=download
#tar -xzf genomes.tar.gz

export PYTHONPATH=$PYTHONPATH:../hannsoden-bioinformatics:..

echo "Generating extra tracks of data"
python extra_tracks.py

echo "Running BRAG"
../BRAG -t brag_set.nwk -r crassa -x F_oxysporum Trichoderma_sp -o Ncrassa -g alignments/brag_set.*.os.txt -q alignments/brag_set.*.seqs -C Neurospora-crassa_v12_fixed_centromeres.txt -T Ncra_extra_tracks.tsv -W 20000 -S 4000

echo "Generating figures from BRAG output"
python paper_figures.py

# This will recreate the figures from the BRAG Paper.
# This analysis uses a few working-names for species that were changed by the time of publication.
# western_mountain = N. sierra
# alaska = N. glassae
# LNF6-4 = FGSC26633
# discreta_Olivares = N. discreta FGSC26659

# Download genomes from Open Science Framework
echo "Downloading genomes"
wget -O genomes.tar.gz https://osf.io/n3rcw/download
tar -xzf genomes.tar.gz
echo

export PYTHONPATH=$PYTHONPATH:../hannsoden-bioinformatics:..

echo "Generating extra tracks of data"
python extra_tracks.py
echo

echo "Running BRAG"
../BRAG -t Neurospora.nwk -r crassa -x S_macrospora -o Ncrassa -g alignments/brag_set.*.os.txt -q alignments/brag_set.*.seqs -C Neurospora-crassa_v12_fixed_centromeres.txt -T Ncra_extra_tracks.tsv -W 20000 -S 4000
echo

echo "Generating figures from BRAG output"
python paper_figures.py
echo

echo "Analyzing distribution of conserved genes"
python conserved_genes.py
echo

echo "Bundling results into ./results"
python bundle_results.py

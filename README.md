# Mapping-of-transcriptome-data-to-genome-scale-scale-model-reactions

####################################

src/GSMM_transcriptome_data_mapper/transcriptome_mapper.py contains the map_transcriptome_data function

####################################

Author: Soukaina Timouma
E-mail: soukaina.timouma@gmail.com


The purpose of this function is to enable the mapping of transcriptome data into a genome scale metabolic model reactions, enabling the imposition of restrictions that accurately represent the experimental conditions being investigated.


map_transcriptome_data() required 4 positional arguments: 'model', 'transcriptomeData', 'threshold_abundance', and 'max_bound'

- 'model': genome scale model
- 'transcriptomeData': dictionary containing the gene IDs (as written in the genome scale model) as keys and the transcription levels as values.
- 'threshold_abundance': threshold to consider that a gene is expressed rather than be noise. For example you can set the threshold to 10 as recommended in DESeq2 documentation.
- 'max_bound': the value of the upper bounds when there is no restriction. For example, in Yeast8 model, the lower and upper bounds range from -1000 to 1000. For Yeast8 model, the 'max_bound' is 1000.


Installation:
python3 -m pip install GSMM_transcriptome_data_mapper


Usage: on a python3 terminal

from GSMM_transcriptome_data_mapper import transcriptome_mapper
transcriptome_mapper.map_transcriptome_data(model,transcriptomeData,threshold_abundance,max_bound)



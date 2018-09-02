# Modisco vs. Homer Comparison

## Homer Results 

positives are created using 1k-surround the summits, negatives are subselected from 500k randomly chosen 1k regions on the genome with dinuc frequencies matching the positives, using Anna's gen_dinucleotide_freqs.py. Target neg:pos ratio set to 1.5, but after removing the duplicates, this ratio turns out to be lower.

[JUND_Homer](https://github.com/kundajelab/tfmodisco_bio_experiments/blob/master/results/nandi/JUND/JUND_GM12878_HOMER_18_08_31/homer/homer_out/homerResults.html)

[RELA_Homer](https://github.com/kundajelab/tfmodisco_bio_experiments/blob/master/results/nandi/RELA/RELA_GM12878_HOMER_18_08_31/homer/homer_out/homerResults.html)

[SPI1_Homer](https://github.com/kundajelab/tfmodisco_bio_experiments/blob/master/results/nandi/SPI1/SPI1_GM12878_HOMER_18_08_31/homer/homer_out/homerResults.html)

[STAT1_Homer](https://github.com/kundajelab/tfmodisco_bio_experiments/blob/master/results/nandi/STAT1/STAT1_GM12878_HOMER_18_08_31/homer/homer_out/homerResults.html)

## Modisco Results with stride 1
Standard pipeline with stride=1.

Hmm, using stride=1 did not deliver major improvements compare to stride=10.

[Summary single tasks training](single_tasks_18_08_31.tsv)


[~~JUND model validation~~](../../JUND/JUND_GM12878_18_08_31/JUND_GM12878_18_08_31.tsv), 
[~~JUND interpretation~~](../../JUND/JUND_GM12878_18_08_31/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb) two tasks merged as one (in progress)

[~~RELA model validation~~](../../RELA/RELA_GM12878_18_08_31/RELA_GM12878_18_08_31.tsv),
[~~RELA interpretation~~](../../RELA/RELA_GM12878_18_08_31/modisco.run1/tfmodisco-visualization-RELA-GM12878.ipynb) (in progress)

[SPI1 model validation](../../SPI1/SPI1_GM12878_18_08_31/SPI1_GM12878_18_08_31.tsv),
[SPI1 training details](../../SPI1/SPI1_GM12878_18_08_31/logs/analyze.txt),
[SPI1 interpretation](../../SPI1/SPI1_GM12878_18_08_31/modisco.run1/tfmodisco-visualization-SPI1-GM12878.ipynb)

[STAT1 model validation](../../STAT1/STAT1_GM12878_18_08_31/STAT1_GM12878_18_08_31.tsv),
[STAT1 training details](../../STAT1/STAT1_GM12878_18_08_31/logs/analyze.txt),
[STAT1 interpretation](../../STAT1/STAT1_GM12878_18_08_31/modisco.run1/tfmodisco-visualization-STAT1-GM12878.ipynb)


## Modisco Results when fine-tuning on same data that was fed to Homer (in progress)


## Older results with stride 10

[Summary single tasks training](../../bQTL/bQTL_18_08_24/single_tasks_18_08_25.tsv)


[JUND model validation](../../JUND/JUND_GM12878_18_08_25/JUND_GM12878_18_08_25.tsv), 
[JUND training details](../../JUND/JUND_GM12878_18_08_25/logs/analyze.txt),
[JUND interpretation](../../JUND/JUND_GM12878_18_08_25/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb) two tasks 
[JUND bQTL enrichment](../../JUND/JUND_GM12878_18_08_25/modisco.run1/SNP_enrichment.ipynb)

[RELA model validation](../../RELA/RELA_GM12878_18_08_25/RELA_GM12878_18_08_25.tsv), 
[RELA training details](../../RELA/RELA_GM12878_18_08_25/logs/analyze.txt),
[RELA interpretation](../../RELA/RELA_GM12878_18_08_25/modisco.run2/tfmodisco-visualization-RELA-GM12878.ipynb),
[RELA bQTL enrichment](../../RELA/RELA_GM12878_18_08_25/modisco.run2/SNP_enrichment.ipynb)

[SPI1 model validation](../../SPI1/SPI1_GM12878_18_08_25/SPI1_GM12878_18_08_25.tsv),
[SPI1 training details](../../SPI1/SPI1_GM12878_18_08_25/logs/analyze.txt),
[SPI1 interpretation](../../SPI1/SPI1_GM12878_18_08_25/modisco.run1/tfmodisco-visualization-SPI1-GM12878.ipynb),
[SPI1 bQTL enrichment](../../SPI1/SPI1_GM12878_18_08_25/modisco.run1/SNP_enrichment.ipynb)

[STAT1 model validation](../../STAT1/STAT1_GM12878_18_08_25/STAT1_GM12878_18_08_25.tsv),
[STAT1 training details](../../STAT1/STAT1_GM12878_18_08_25/logs/analyze.txt),
[STAT1 interpretation](../../STAT1/STAT1_GM12878_18_08_25/modisco.run2/tfmodisco-visualization-STAT1-GM12878.ipynb),
[STAT1 bQTL enrichment](../../STAT1/STAT1_GM12878_18_08_25/modisco.run2/SNP_enrichment.ipynb)



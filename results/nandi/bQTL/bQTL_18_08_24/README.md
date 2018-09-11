# bQTL multi-task and single-tasks models
## bQTL multi-task model
[Summary of bQTL multi-TF training](bQTL_18_08_24.tsv)

[interpretation1](modisco.run1/tfmodisco-visualization-bQTL-GM12878.ipynb) deeplift+modisco on train+valid set

[interpretation2](modisco.run2/tfmodisco-visualization-bQTL-GM12878.ipynb) deeplift+modisco on 1kb sequences centered at summit, excluding chr1 (in progress)

## bQTL single-task models

deeplift+modisco on 1k surround summit

[Summary of single task trainings](single_tasks_18_08_25.tsv)

[JUND model training](../../JUND/JUND_GM12878_18_08_25/JUND_GM12878_18_08_25.tsv), 
[JUND interpretation](../../JUND/JUND_GM12878_18_08_25/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb) two tasks together

[RELA model training](../../RELA/RELA_GM12878_18_08_25/RELA_GM12878_18_08_25.tsv), 
[RELA interpretation](../../RELA/RELA_GM12878_18_08_25/modisco.run2/tfmodisco-visualization-RELA-GM12878.ipynb)

[SPI1 model training](../../SPI1/SPI1_GM12878_18_08_25/SPI1_GM12878_18_08_25.tsv),
[SPI1 interpretation](../../SPI1/SPI1_GM12878_18_08_25/modisco.run1/tfmodisco-visualization-SPI1-GM12878.ipynb)

[STAT1 model training](../../STAT1/STAT1_GM12878_18_08_25/STAT1_GM12878_18_08_25.tsv),
[STAT1 interpretation](../../STAT1/STAT1_GM12878_18_08_25/modisco.run2/tfmodisco-visualization-STAT1-GM12878.ipynb)


## bQTL multi-task model refined on single-task data

deeplift+modisco on 1k surround summit

[Summary of single task refinement trainings](refine_tasks_18_09_04.tsv).

[JUND model training](../../JUND/JUND_GM12878_refine_18_09_04/JUND_GM12878_refine_18_09_04.tsv), 
[JUND interpretation](../../JUND/JUND_GM12878_refine_18_09_04/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb) two tasks together,
[JUND SNP enrichment](../../JUND/JUND_GM12878_refine_18_09_04/modisco.run1/SNP_enrichment.ipynb). 
Took [JUND multi-task model](../../JUND/JUND_18_09_03/finetune/model_files/), refined on [JUND GM12878 data](../../JUND/JUND_GM12878_18_08_31).

[RELA model training](../../RELA/RELA_GM12878_refine_18_09_04/RELA_GM12878_refine_18_09_04.tsv), 
[RELA interpretation](../../RELA/RELA_GM12878_refine_18_09_04/modisco.run2/tfmodisco-visualization-RELA-GM12878.ipynb), 
[RELA SNP enrichment](../../RELA/RELA_GM12878_refine_18_09_04/modisco.run1/SNP_enrichment.ipynb). 
Took [RELA multi-task model](../../RELA/RELA_18_09_03/finetune/model_files/), refined on [RELA GM12878 data](../../RELA/RELA_GM12878_18_08_31).

[SPI1 model training](../../SPI1/SPI1_GM12878_refine_18_09_04/SPI1_GM12878_refine_18_09_04.tsv),
[SPI1 interpretation](../../SPI1/SPI1_GM12878_refine_18_09_04/modisco.run1/tfmodisco-visualization-SPI1-GM12878.ipynba),
[SPI1 SNP enrichment](../../SPI1/SPI1_GM12878_refine_18_09_04/modisco.run1/SNP_enrichment.ipynb). 
Took [SPI1 multi-task model](../../SPI1/SPI1_18_09_03/finetune/model_files/), refined on [SPI1 GM12878 data](../../SPI1/SPI1_GM12878_18_08_31).

[STAT1 model training](../../STAT1/STAT1_GM12878_refine_18_09_04/STAT1_GM12878_refine_18_09_04.tsv),
[STAT1 interpretation](../../STAT1/STAT1_GM12878_refine_18_09_04/modisco.run2/tfmodisco-visualization-STAT1-GM12878.ipynb),
[STAT1 SNP enrichment](../../STAT1/STAT1_GM12878_refine_18_09_04/modisco.run1/SNP_enrichment.ipynb). 
Took [STAT1 multi-task model](../../STAT1/STAT1_18_09_03/finetune/model_files/), refined on [STAT1 GM12878 data](../../STAT1/STAT1_GM12878_18_08_31).

## bQTL multi-TF model, with single TF interpretation: run tfmodisco on only the deeplift scores for that TF/task:
[interpretation task1](modisco_task1/tfmodisco-visualization-bQTL-GM12878-JUND.ipynb) deeplift+modisco task1 (JUND)

[interpretation task2](modisco_task2/tfmodisco-visualization-bQTL-GM12878i-RELA.ipynb) deeplift+modisco task2 (RELA)

[interpretation task3](modisco_task3/tfmodisco-visualization-bQTL-GM12878.ipynb) deeplift+modisco task3 (SPI1)

[interpretation task4](modisco_task4/tfmodisco-visualization-bQTL-GM12878.ipynb) deeplift+modisco task4 (STAT1)




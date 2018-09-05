|  |Comparison of modisco approaches using JUND |
|---|--------------------------------------------|
|1| Single TF single task model (single cell-line) [JUND_GM12878_18_08_31 stride=1 modisco](../JUND_GM12878_18_08_31/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb), [JUND_GM12878_18_08_25 stride=10 modisco](../JUND_GM12878_18_08_25/modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb)|
|2| Single TF multi-task model (multiple cell-line)[JUND_18_09_03 stride=10](../JUND_18_09_03), [JUND_18_07_31 stride=20](../JUND_18_07_31)|
|3| Single TF multi-task model finetune on single task [JUND_GM12878_18_09_04 modisco](modisco.run1/tfmodisco-visualization-JUND-GM12878.ipynb) which takes [JUND_18_09_03](../JUND_18_09_03), and finetune on [JUND_GM12878_18_08_31](../JUND_GM12878_18_08_31) |
|4| Multi-TF multi-task model [bQTL_18_08_24 stride=10 modisco](../../bQTL/bQTL_18_08_24/modisco.run2/tfmodisco-visualization-bQTL-GM12878.ipynb)|

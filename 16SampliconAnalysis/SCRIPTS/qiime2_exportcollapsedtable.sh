
# Reference:
# https://forum.qiime2.org/t/lefse-after-qiime2/4496/8

module load QIIME2/2018.4 R/3.4.1-foss-2016b-X11-20160819
module list

# collapse the table.gza to the level 4 (Order)
qiime taxa collapse --i-table table-notrim.qza --o-collapsed-table collapse.table.qza --p-level 4 --i-taxonomy taxonomy-notrim.qza

# calculate relative-frequency for the collapsed table 
qiime feature-table relative-frequency --i-table collapse.table.qza --o-relative-frequency-table collapse.frequency.table.qza --output-dir collapse.frequency/

# export biom table
qiime tools export collapse.frequency.table.qza --output-dir collapse.frequency/

# convert biom to txt file
biom convert --input-fp collapse.frequency/feature-table.biom --output-fp collapse.frequency/collapsed-feature-table.txt --header-key “taxonomy” --to-tsv


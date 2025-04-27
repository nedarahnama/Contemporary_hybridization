# /home/alle/Neda/repeat_analysis/parents # working directory on local PC
- input: 
07_3_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci.vcf


conda install conda-forge::py-bgzip

bgzip file.vcf

bcftools index nem_sag_onlypolymorphic.vcf.gz
bcftools index 07_3_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci.vcf.gz

bcftools isec -n 2 -p shared_variants/ 07_3_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci.vcf.gz nem_sag_onlypolymorphic.vcf.gz


#bcftools merge ./nem_sag_onlypolymorphic.vcf.gz ./07_3_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci.vcf.gz -O z -o mergedF2_parents.vcf.gz

#bcftools isec -n=2 mergedF2_parents.vcf.gz -p shared_variants

#bcftools concat shared_variants/0000.vcf.gz shared_variants/0001.vcf.gz -Oz -o shared_variants.vcf.gz


- Rename vcf file of F2 that contains variants that are present in parents:

/home/alle/Neda/repeat_analysis/parents/shared_variants/0000.vcf -> 08_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci_sharedvariants.vcf

6014 sites 

- copy file to cluster:
scp 08_miss0.5_indels_biallelic_0.6indmiss_4meanDP40_5DP40_miss0.7_RADloci_sharedvariants.vcf nrahnama@cheops1.rrz.uni-koeln.de:/scratch/nrahnama1/bcf_new

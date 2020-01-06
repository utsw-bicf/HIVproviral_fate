#! /bin/bash
#### This script gets the number of HIV insertion in Loop vs Not loop
#### And it also gets HIV insertion in Tad vs not tad

##### Loop
# Count HIV in Loop; result:974
awk '{OFS="\t"}; {print $6,$7,$8}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/HIVexpression_scored_loop.bed | sort -k 1,1 -k 2,2n | uniq | wc -l
# Count HIV not in Loop; result:587
wc -l ../../../../Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed
  # 1561-974=587

# Sum up Loop; result:1896895500
awk '{OFS "\t"} {print $1,$2,$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.loop.2D_filtered.bed | bedops --merge - | awk '{sum+=($3-$2)} END {print sum}'
# Sum up not loop; result: 1122512763
awk '{sum+=($3-$2)} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_loop.bed

  ## Total sum:3019408263

##### Tad
# Count HIV in tad; result:731
awk '{OFS="\t"}; {print $6,$7,$8}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/HIVexpression_scored_tad.bed | sort -k 1,1 -k 2,2n | uniq | wc -l
# Count HIV not in Loop; result:830
wc -l ../../../../Bhive/bhive_singularity/BHIVE_for_single_provirus_transcriptomics/hiv_expression.bed
  # 1561-731=830

# Sum up Tad; result:1363364780
awk '{OFS "\t"} {print $1,$2,$3}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/lib.tad.2D_filtered.bed | bedops --merge - | awk '{sum+=($3-$2)} END {print sum}'
# Sum up not loop; result: 3007084417
awk '{sum+=($3-$2)} END {print sum}' /project/BICF/BICF_Core/shared/Projects/Dorso/Chromatin_organization/HOMER/std_HicHomerTagDir/merge_lib/not_tad.bed

  ## Total sum:4370449197

### Run hypergeometric test in R
### Did below in command
module load R

>library(broom)

>pt <- tidy(prop.test(x=c(974,1896895500),n=c(1561,3019408263), p=NULL, alternative="two.sided", conf.level=0.95, correct=F))
> pt
# A tibble: 1 x 9
  estimate1 estimate2 statistic p.value parameter conf.low conf.high method
      <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>     <dbl> <chr> 
1     0.624     0.628     0.122   0.727         1  -0.0283    0.0198 2-sam…
# ... with 1 more variable: alternative <chr>
> 


>pt2 <- tidy(prop.test(x=c(731,1363364780),n=c(1561,4370449197), p=NULL, alternative="two.sided", conf.level=0.95, correct=F))
> pt2
# A tibble: 1 x 9
  estimate1 estimate2 statistic  p.value parameter conf.low conf.high method
      <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>     <dbl> <chr> 
1     0.468     0.312      178. 1.50e-40         1    0.132     0.181 2-sam…
# ... with 1 more variable: alternative <chr>


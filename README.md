# selection_scan
Scripts for genome-wide calculation of population genetic metrics, including FST, DXY, and π, for population genomics analyses.
## Description
This repository contains scripts for population genomic analyses in rice, including:
- Genome-wide calculation of population genetic metrics (FST, DXY, π, θ, Tajima's D)
- Haplotype analysis for specific genes, such as MADS51

The scripts are a mix of Perl and Python programs and can be adapted to other species or datasets with phased genotype or allele frequency files.

## Repository Contents

### Perl Scripts
1. **1_cal.fr.pl**  
   - Calculates allele frequencies from input data.
2. **2_get.pos.idx.pl**  
   - Generates an index of SNP positions to speed up downstream calculations.
3. **3_cal.seg.swd.step.pop.pl**  
   - Performs sliding-window analysis for a given genomic segment and population data.
4. **4_cal.gene.pop.pl**  
   - Summarizes per-gene population genetic metrics across groups.

### Python Script
- **haplotype_analysis.py**  
  - Generates phased haplotypes from a VCF file.
  - Calculates haplotype frequencies per group and per sample.
  - Produces a summary output with haplotype distributions.

## Usage

perl 1_cal.fr.pl <input_file> <output_file>
perl 2_get.pos.idx.pl <input_file> 
perl 3_cal.seg.swd.step.pop.pl <segment_name> <chr> <start> <end> <window> <step> <output>
perl 4_cal.gene.pop.pl <up_length> <down_length> <output>
python3 haplotype_analysis.py


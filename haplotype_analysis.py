#!/usr/bin/env python3
"""
MADS51 Haplotype Analysis in Rice

This script parses a VCF file with phased genotypes and generates haplotypes
for each sample, calculates haplotype frequencies, and outputs distributions
across predefined rice groups.

Inputs:
- xx.vcf : VCF file containing phased genotypes (using |)
- xx.list : tab-delimited file mapping sample names to population groups

Output:
- xx.out : summary of haplotypes and their distribution
"""

# Predefined rice groups
gpn = ["Niv1a", "Niv1b", "Niv1c", "Niv2a", "Niv2b", "Niv2c",
       "Ruf1a", "Ruf1b", "Ruf2a", "Ruf2b", "Ruf2c", "Ruf2d"]

# Load sample-to-group mapping
with open("xx.list", 'r') as f:
    grp = {line.split('\t')[0]: line.split('\t')[1].strip() for line in f}

# Read VCF file
with open('xx.vcf', 'r') as f:
    vcf_lines = []
    header_line = ""
    for line in f:
        if line.startswith("#CHROM"):
            header_line = line.strip()
        if line.startswith("#"):
            continue
        vcf_lines.append(line.strip())

# Get sample names from VCF
vcf_cols = header_line.split('\t')
samples = vcf_cols[9:]
sam_pos = {s: i + 9 for i, s in enumerate(samples)}

# Filter samples to include only those in predefined groups
sam = [s for s in samples if grp.get(s) in gpn]

# Organize samples into groups
pgp = [[] for _ in range(len(gpn))]
for s in sam:
    for i, group in enumerate(gpn):
        if grp[s] == group:
            pgp[i].append(s)

# Initialize haplotype sequences
hap_phs = {s + '.1': "" for s in sam}
hap_phs.update({s + '.2': "" for s in sam})

# Process each variant in VCF
for var in vcf_lines:
    cols = var.split('\t')
    ref = cols[3]
    alt_list = cols[4].split(',')
    for sample in sam:
        genotype = cols[sam_pos[sample]].split('|')
        # phased haplotype sequence
        if len(alt_list) == 1:
            # Biallelic site
            if genotype == ['0', '0']:
                hap_phs[sample + '.1'] += ref
                hap_phs[sample + '.2'] += ref
            elif genotype == ['0', '1']:
                hap_phs[sample + '.1'] += ref
                hap_phs[sample + '.2'] += alt_list[0]
            elif genotype == ['1', '0']:
                hap_phs[sample + '.1'] += alt_list[0]
                hap_phs[sample + '.2'] += ref
            elif genotype == ['1', '1']:
                hap_phs[sample + '.1'] += alt_list[0]
                hap_phs[sample + '.2'] += alt_list[0]
        else:
            # Multiallelic site
            for h in [0, 1]:
                if genotype[h] == '0':
                    hap_phs[sample + f'.{h+1}'] += ref
                else:
                    hap_phs[sample + f'.{h+1}'] += alt_list[int(genotype[h]) - 1]

# Collect all haplotypes
hap_dat = [hap_phs[s + '.1'] for s in sam] + [hap_phs[s + '.2'] for s in sam]
hap = []
for h in hap_dat:
    if h not in hap:
        hap.append(h)

# Count haplotype frequencies
hap_total = len(hap_dat)
hap_num = {h: hap_dat.count(h) for h in hap}
nhp = [h for h in hap if hap_num[h] / hap_total > 0.0]  # can adjust frequency threshold
nhp_no = {h: f'H{i+1}' for i, h in enumerate(nhp)}

# Write output
with open('xx.out', 'w') as W:
    W.write("Variation of sites\n")
    for i, var in enumerate(vcf_lines):
        W.write(f"{i + 1}\t{var}\n")
    W.write("\n")

    W.write(f"A total of {len(nhp)} gene haplotypes in rice\n")
    for i, h in enumerate(nhp):
        W.write(f"H{i+1}\t{h}\n")
    W.write("\n")

    # Distribution of haplotypes
    W.write(f"The distribution of gene haplotype in rice groups\n")
    W.write("Haplotype\tAll")
    for group in gpn:
        W.write(f"\t{group}")
    W.write("\n")

    for i, h in enumerate(nhp):
        W.write(f"H{i+1}")
        # All samples
        htn = sum(hap_phs[s + '.1'] == h or hap_phs[s + '.2'] == h for s in sam)
        W.write(f"\t{htn}")
        # Per group
        for j, group in enumerate(gpn):
            htn = sum(
                hap_phs[s + '.1'] == h or hap_phs[s + '.2'] == h
                for s in pgp[j]
            )
            W.write(f"\t{htn}")
        W.write("\n")
    W.write("\n")

    # Haplotype assignment for each sample
    W.write("The haplotype of each sample in each group\n")
    for j, group in enumerate(gpn):
        W.write(f"{group}\thaplotype\n")
        for s in pgp[j]:
            W.write(f"{s}_1\t{nhp_no[hap_phs[s + '.1']]}\n")
            W.write(f"{s}_2\t{nhp_no[hap_phs[s + '.2']]}\n")

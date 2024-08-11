## How to Install

```
git clone https://github.com/MichaelSebero/Y-Haplogroup-Predictor

python3 /path/to/y-haplogroup-predictor.py
```
## Description

Y-Haplogroup-Predictor detects SNPs associated with haplogroups when provided raw DNA data which is formatted in rsID format. This script calculates which SNPS are the most prevalent in a sample and decides the haplogroup.

The intention of this software is to provide a FOSS alternative to detecting haplogroups. Databreaches occur often and institutions shouldn't be trusted with processing / storing client genetic data over a long period of time. Y-Haplogroup-Predictor is completely private and can be ran offline without sending PIA to a website.

## How to use Y-Haplogroup-Predictor

- Aquire a copy of your DNA data from a genetic testing platform.
- Open up the DNA file and make sure it's formatted in rsID format.

> The script looks for the rsid and genotype to determine results.

    rsid         chromosome  position   allele1	allele2

    rs17250535   24          22739367   T           T

- Input the raw DNA file pathway into the console.

### This script only detects SNPs associated with haplogroups R1a, R1b, I1 and I2, along with their subclades for now.

<p align="center">
  <img src="https://i.postimg.cc/yd0jYTVL/Europe.png"/>
</p>


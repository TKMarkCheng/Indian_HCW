In-house script for Kemp, S.A., Cheng, M.T.K., Hamilton, W.L. et al. Transmission of B.1.617.2 Delta variant between vaccinated healthcare workers. Sci Rep 12, 10492 (2022). https://doi.org/10.1038/s41598-022-14411-7


# Indian_HCW

## 8th Sept 2021

Uploaded test sequences and the pipeline script.

To run the pipeline, clone the repository, replace the sequences with the desired sequences and run the following code.
if an alternative reference sequence is needed, replace MN908947.3 with the alternative chromosome name in line 9 and 24 of snp_site_matrix_unix.sh accordingly.

```console
bash snp_site_matrix_unix.sh
```

The two python scripts can be ran independently from the script for their respective function.
1) Split.py
```console
python Split.py input.vcf output.vcf
```
Split multiallelic variants to multiple biallelic variant for vcfs generated by snp-sites.

2) matrix.py
```console
python matrix.py input.vcf output.xlsx
```
Takes in biallelic snp-site vcf and generates 4 pairwise-matrixes: N-included snp-count and exact snp-change, and N-excluded snp-count and exact snp-change.

## 9th Sept 2021

Matrix.py
### Corrected behaviour
Now considers multiallelic sites as single mutations rather than two mutations (reverse + forward) compounded.
### Improved functionality
now also outputs a snp-dist matrix for comparison against snp-dists matrix generated by snp-dists.

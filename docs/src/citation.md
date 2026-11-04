# Citation

If you use HapApp_utils in your work, please cite the paper that describes the pipeline and its databases.

## Citing HapApp_utils

```text
Zhao D, Lin M, Taniguti CH, Sandercock AM, et al. (2026). Standardized microhaplotype
databases and frameworks for assessing and mining crop genetic diversity.
Theoretical and Applied Genetics 139(9):244.
https://doi.org/10.1007/s00122-026-05340-4
```

BibTeX:

```bibtex
@article{zhao2026microhaplotype,
  author = {Zhao, Dongyan and Lin, Meng and Taniguti, Cristiane H. and Sandercock, Alexander M. and Chen, Shufen and Babiker, Ebrahiem and Bassil, Nahla V. and Brummer, E. Charles and Camacho, Jose Roberto and Chatwin, Warren and Chen, Shu-Yun and Clare, Shaun J. and da Silva Pereira, Guilherme and David, Maria and Fraher, Simon Phillip and Hardigan, Michael and Hilton, Angelyn and Hislop, Lillian M. and Irish, Brian M. and Kante, Moctar and Kim, Tae Hwa and Lee, Chong-Wei and Lindqvist-Kreuze, Hannele and Loarca, Jenyne and Lu, Po-Hsien and Medina Culma, Cesar Agusto and Jim\'{e}nez Morales, Jose Fabi\'{a}n and Polashock, James and Price, John H. and Riday, Heathcliffe and Samac, Deborah A. and Sandhu, Devinder and Ssali, Reuben and Castro V\'{a}squez, Ruth Mayela and Wadl, Phillip A. and Wang, Xinwang and Webster, Seymour A. and Xu, Zhanyou and Yencho, G. Craig and Beil, Craig T. and Sheehan, Moira J.},
  title = {Standardized microhaplotype databases and frameworks for assessing and mining crop genetic diversity},
  journal = {Theoretical and Applied Genetics},
  year = {2026},
  volume = {139},
  number = {9},
  pages = {244},
  doi = {10.1007/s00122-026-05340-4}
}
```

## AI Assistance

AI assistance applies to version 2 of HapApp_utils and its documentation, not to version 1. Dongyan Zhao wrote version 1,
the step scripts and bash workflows. For version 2, Claude Opus 5.5 (Anthropic) helped convert those workflows into the
Python command-line tool, fix bugs in the step scripts, and write the tests and documentation. The authors directed the
work, reviewed and tested the results, and are ultimately responsible for the application.

## Related Resources

- Breeding Insight, [MADC file format](https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/)
- Breeding Insight, [MADC data processing](https://breeding-insight.github.io/BIGr_documentation/posts/madc-data-processing/)
- Breeding Insight, [MADCATS](https://madcats.readthedocs.io), which generates MADC tables from FASTQ files

## Example Acknowledgment

In your methods section:

> We assigned fixed microhaplotype IDs with HapApp_utils v2.0.0 (Zhao et al., 2026), keeping alleles with at least
> 2 reads in at least 10 samples or 5% of samples, whichever was smaller, and requiring 90% BLAST coverage and identity
> against the panel's allele database.

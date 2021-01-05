# DeDup
A merged read deduplication tool capable to perform merged read deduplication on paired-end sequencing data on BAM files. 

[![Build Status](https://travis-ci.org/apeltzer/DeDup.svg?branch=master)](https://travis-ci.org/apeltzer/DeDup)
[![codecov](https://codecov.io/gh/apeltzer/DeDup/branch/master/graph/badge.svg)](https://codecov.io/gh/apeltzer/DeDup)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/dedup/README.html)

Author: Alexander Peltzer <alexander.peltzer@uni-tuebingen.de>

###### Method
This DeDup tool is a PCR duplicate removal tool of paired-end and merged sequenced data designed for ultra-short DNA (e.g. ancient DNA).

The basic concept is that when dealing with ultra-short DNA, with typical sequencing chemistry cycles, you will sequence the entire molecule and therefore find both ends. 

Compared to typical deduplication tools, that only look for reads with the same starting position at the 5' end of a read, DeDup will remove 'true' duplicates using BOTH (5' and 3') ends of the reads. This can help increase coverage in low-preservation samples such used in ancient DNA by being more exact as to what are duplicates or not.

While the DeDup tool is designed to work primarily with merged reads, it can also account for additional unmerged reads in the same library. DeDup expects the different kinds of reads to have read names that begin with one of the following prefixes:

- `F_`
- `R_`
- `M_`

To remove PCR duplicates we only retain a single read for a given genomic position and read direction. For `M_` (merged) reads we know both the start and end of the sequenced fragment. Two `M_` reads with the same start end and end position will be deduplicated (i.e. best quality read retained). 

For F_ and R_ reads we only know the start or end of the sequenced fragment because the read length is variable. In this case, positions of only one end of the read will be considered for deduplication (as with standard deduplication tools). In the case of a mixture of two reads types, e.g. `M_` and `F_`, only the start position will be considered and `M_` reads will be preferred (with the assumption of higher quality reads). Note this could sometimes lead to slightly over-zealous deduplication, as the unmerged read (`F_`, in this example) may not truly have the same end position as the `M_`.

:warning: **DeDup is NOT designed for solely single-end data**. If single-end data is supplied and the `--merged` flag is used, both start and end positions will _still_ be considered for deduplication. In this case, you will sub-optimally deduplicate your reads as single-end data will not necessarily always represent the true 'end' of the molecule and therefore leaving possible duplicates in your BAM file. 

A little documentation is available at http://dedup.readthedocs.io/en/latest/

###### Citation

DeDup was originally published in:

Peltzer, A., Jäger, G., Herbig, A., Seitz, A., Kniep, C., Krause, J., & Nieselt, K. (2016). EAGER: efficient ancient genome reconstruction. Genome Biology, 17(1), 1–14. [https://doi.org/10.1186/s13059-016-0918-z](https://doi.org/10.1186/s13059-016-0918-z)

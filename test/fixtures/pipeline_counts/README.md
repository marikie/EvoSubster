# Historical pipeline counting baseline

These eight TSVs were copied byte-for-byte from
`test_results/Oikdio1_Oikalb2_Oikvan3/20260420/statistics/` in the source checkout.
They are historical expectations, not outputs regenerated with current code.

The outgroup was GCA_907165135.1 (Oikopleura dioica); the ingroups were
GCA_004367875.1 (Oikopleura albicans) and GCA_004367855.1 (Oikopleura vanhoeffeni).
Both whole-genome and non-CDS counts are required for single-base and dinucleotide
substitutions. The original run used the outgroup's `genomic.gff` annotation.
Original tool versions and input checksums were not recorded with this baseline.

`scripts/codex_checks.py` compares all eight files exactly, mapping only the run
date in filenames. It never rewrites these expectations to accommodate a failure.
Genomes, alignments, and training caches are not stored here.

Original SHA-256 checksums:

| File (relative to this directory) | SHA-256 |
| --- | --- |
| Oikalb2/dinuc/GCA_004367875.1_Oikalb2_20260420_dinuc.tsv | 438d235a90d6d5e711e709cad6e4004f9bbb8338469590830df47702423c6a52 |
| Oikalb2/dinuc/GCA_004367875.1_Oikalb2_20260420_dinuc_ncds.tsv | 5c9aac982fe5f3f546d2b35979a2ccfa4bbfbdf6f9f9986ebba1ce35b1dfdad8 |
| Oikalb2/singlenuc/GCA_004367875.1_Oikalb2_20260420.tsv | b751feb7e9d6087642d03239bc4b03eba0d305f6aa6dff012b2438989b0580ac |
| Oikalb2/singlenuc/GCA_004367875.1_Oikalb2_20260420_ncds.tsv | 01a3b061825e9f9514225785221b23d43ba0d8ab103daeb0e3db33535891f522 |
| Oikvan3/dinuc/GCA_004367855.1_Oikvan3_20260420_dinuc.tsv | 2e54dcb7c9245f431303984fd498c5b6b293c21662842fda5740ade74ffb6fb8 |
| Oikvan3/dinuc/GCA_004367855.1_Oikvan3_20260420_dinuc_ncds.tsv | aa02b847ad6d1a19c8e0cf59386a1aea24a01bb2c3fc6abd50bc312511d7518c |
| Oikvan3/singlenuc/GCA_004367855.1_Oikvan3_20260420.tsv | 5f13e7bb371f6c8252b15805b0a698bc5080caf3a3410106f0a73008725a1ca3 |
| Oikvan3/singlenuc/GCA_004367855.1_Oikvan3_20260420_ncds.tsv | 997779763ba844d57349b959a27ec09d1faeabbf58313461bf05d16951a57aab |

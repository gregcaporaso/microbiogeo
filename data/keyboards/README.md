Keyboards
=========

QIIME commands
--------------

The following commands were run with QIIME 1.8.0:

```bash
# Original OTU table is table.biom.orig. Filtered out all samples that didn't
# come from the three primary subjects (M2, M3, M9). Also did this for the
# rarefied tables, which were created from the original, unfiltered OTU table.
filter_samples_from_otu_table.py -i table.biom.orig -o table.biom -m map.txt -s 'HOST_SUBJECT_ID:232:M2,232:M3,232:M9'

```

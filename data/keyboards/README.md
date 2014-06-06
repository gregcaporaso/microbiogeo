Keyboards
=========

QIIME commands
--------------

TODO
----

Once we have the new table and mapping file, need to do the following
filtering:

Original OTU table is otu_table.biom.orig
Filtered out all samples that didn't come from the three primary subjects
(M2, M3, M9).

    filter_samples_from_otu_table.py -i otu_table.biom.orig -o otu_table.biom -m map.txt -s 'HOST_SUBJECT_ID:M2,M3,M9'

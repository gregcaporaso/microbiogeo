Keyboards
=========

QIIME commands
--------------

The following commands were run with QIIME 1.8.0:

```bash

# Original OTU table is table.biom.orig. Filter out all samples that didn't
# come from the three primary subjects (M2, M3, M9). Also do this for the
# rarefied tables, which were created from the original, unfiltered OTU table.
filter_samples_from_otu_table.py -i table.biom.orig -o table.biom -m map.txt -s 'HOST_SUBJECT_ID:232:M2,232:M3,232:M9'

for DEPTH in 200 500 1000
do
  for REP in {0..9}
  do
    INFP=rarefied-orig/${DEPTH}/rarefaction_${DEPTH}_${REP}.biom
    OUTFP=rarefied/${DEPTH}/rarefaction_${DEPTH}_${REP}.biom
    filter_samples_from_otu_table.py -i ${INFP} -o ${OUTFP} -m map.txt -s 'HOST_SUBJECT_ID:232:M2,232:M3,232:M9'
  done
done

```

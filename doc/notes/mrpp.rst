===========================================================================
Multiple Response Permutation Procedure (MRPP) Statistical Method Reference
===========================================================================

Introduction
------------
The Multiple Response Permutation Procedure (MRPP) tests whether there is a
significant difference between two or more groups of samples. It is similar to
ANOSIM except that it uses the original distances instead of ranks,
so it may be more sensitive to outliers. It is also conceptually similar to
ANOVA in that it compares distances within and between groups. MRPP is
non-parametric. The original MRPP paper can be referenced here
[:ref:`2 <mrppref2>`].

MRPP calculates a delta statistic, which is the overall weighted mean of
within-group means of the pairwise distances among samples. The samples and
their pairwise distances are then permuted and delta is calculated for each
permutation of the data. The significance test is the fraction of permuted
deltas that are less than the observed delta. If two groups of samples are
really different, the average of the within-group distances should be less than
the average of the distances between two random collection of samples drawn from
the entire population [:ref:`1 <mrppref1>`].

MRPP tests whether there are differences between two or more groups of samples,
but it may not always find differences in groups due to differences in means.
Instead, MRPP may find differences in groups based on spread (differences in
within-group distance). For example, it might find that two groups are
significantly different because one group has greater within-group distances
than the other. Thus, it has been recommended to use adonis instead of MRPP when
possible because adonis doesn't have this problem [:ref:`1 <mrppref1>`].

Existing Implementations
------------------------
There are several existing implementations of MRPP in the following statistical
packages:

* vegan package for R

* MRPP macro for SPSS package [:ref:`3 <mrppref3>`]

* FORTRAN implementation from original author

* possibly others...

As the vegan package has an MRPP implementation and it is free and open-source,
this implementation will be tested. A simple R script wrapping vegan's mrpp
function has been checked into the Qiimeutils repository under
:file:`microbiogeo/r/examples/`. The following sections of the document will
explain how to run the R script.

System Setup and Required Dependencies
--------------------------------------
:note: The following instructions have been tested on 64-bit Linux Mint (essentially Debian). However, they `should` work across different Linux distros and on Macs, though some commands may need to be tweaked, or different package names might have to be used. The instructions assume you use bash as your shell.

The first step is to install R. The following command downloaded and installed R
(for me, it was R version 2.13.1): ::

    sudo apt-get install r-base

Next, you must install the vegan and optparse packages in R. Run the following
commands: ::

    sudo R
    install.packages("vegan")
    install.packages("optparse")
    q()

The install process for the packages will prompt you to choose a mirror to
download them from. Other than that, it is completely automated. On my system, I
ended up with vegan version 2.0-2 and optparse version 0.9.4.

Next, your system must have a version of QIIME installed (I used the latest
version of QIIME in SVN). The MRPP script uses some R utility functions in QIIME
to load data.

Next, you must define an environment variable to tell the MRPP script where to
look for the R utility functions in QIIME. Run the following command, changing
the path to point to the location of your QIIME install: ::

    export QIIME_DIR=/home/jrideout/qiime/trunk

If you don't want to have to perform this step each time you open a new
terminal, run the following command to add it to your .bashrc: ::

    echo "export QIIME_DIR=/home/jrideout/qiime/trunk" >> ~/.bashrc
    source ~/.bashrc

Next, run the following command to test if you can run the MRPP script: ::

    R --slave --args -h < r/examples/mrpp.r

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The MRPP script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix from the QIIME overview tutorial. You can get the distance
matrix :download:`here <../downloads/overview_unweighted_unifrac_dm.txt>` and
the mapping file :download:`here <../downloads/Fasting_Map.txt>`.

Next, run the following command to execute the MRPP script: ::

    R --slave --args -d overview_unweighted_unifrac_dm.txt -m Fasting_Map.txt -c Treatment < r/examples/mrpp.r

The -c option specifies which column in the mapping file will be used to group
the samples. The `Treatment` column has two values: `Control` and `Fast`. Thus,
MRPP will be used to calculate the dissimilarity between the control and fast
groups.

Output Files
------------
The command in the previous section creates a single output file in the current
directory named :file:`mrpp_results.txt`. The resulting file should look like
this: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 
    
    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          Control Fast  
    delta 0.6237  0.6243
    n     5       4     

    Chance corrected within-group agreement A: 0.07164 
    Based on observed delta 0.624 and expected delta 0.6721 

    Significance of delta: 0.008 
    Based on  999  permutations

The second from the last line contains the p-value of the observed delta
statistic, which is 0.008. This indicates that the differences between `Control`
and `Fast` sample groups is significant, based on 999 permutations.

Testing Results
---------------
This section will describe different tests that were run on the MRPP script.

:note: Many of these tests will use empirical data from one of the several datasets that the team has access to. These data files will not be included for download due to their (usually) large size, but it should be clear what inputs were used.

References
----------
.. _mrppref1:

[1] R help page for vegan function mrpp

.. _mrppref2:

[2] http://www.jstor.org/stable/1940409

.. _mrppref3:

[3] http://lcai.bol.ucla.edu/programs.html

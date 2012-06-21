=========================
Redundancy Analysis (RDA)
=========================

Synopsis
--------
DB-RDA is an ordination method that shows grouping/clustering of samples based
on a category in the metadata mapping file and a distance matrix. This category
is used to explain the variability between samples. Thus, RDA is similar to PCoA
except that it is constrained, while PCoA is unconstrained (you must specify
which category should be used to explain the variability in your data).

Introduction
------------
Redundancy Analysis (RDA) is an ordination method similar to Canonical
Correspondence Analysis (CCA) and Principal Components Analysis (PCA) that
allows one to reduce the dimentionality of the data to something more
manageable. RDA assumes that species are linearly-related to environmental
gradients [:ref:`3 <rdaref3>`].

RDA can accept one to three matrices as input. X is the community distance
matrix that is the required input to RDA. Y is the constraining matrix
of environmental variables and is optional. It can be used to place constraints
on how variation is explained by RDA. If Y is used, it is called `constrained`
redundancy analysis. Z is the condition matrix and its effects can be removed
from X before X is processed any further. Z is also optional, and if it is used,
it is called `partial` redundancy analysis [:ref:`4 <rdaref4>`].

[:ref:`2 <rdaref2>`] (may need to be on NAU VPN to access PDF) details a very
interesting approach in what they call 'distance-based RDA', or 'DB-RDA' for
short. This method basically takes a distance matrix of any type, performs
principal coordinates analysis (PCoA) on it, corrects for negative eigenvalues,
and then performs RDA on the result and another matrix of environmental dummy
variables to analyze their relationship. This application of RDA might be the
most relevant to QIIME/microbial ecology, so it may be worth looking into.

Another paper challenges the method used in the previously described paper
[:ref:`6 <rdaref6>`]. The authors argue that the negative eigenvalue correction
step is not necessary.

:note: After meeting with our clients, it was determined that I should focus on DB-RDA rather than the other two implementations of traditional RDA since it is designed to work on a distance matrix, not a community data matrix like the others. The other two methods assume Euclidean distances, which may not always be useful for our purposes.

More research into DB-RDA yielded some more options that we might consider
including in our own implementation of the method, or at least test them out
during the evaluation phase. The first option deals with how many environmental
variables (i.e. categories in the mapping file) should be used as constraints.
The author of capscale recommends to ``not`` add too many constraining variables
because the DB-RDA ends up becoming unconstrained ordination
(e.g. PCA, CA, etc.) and you cannot perform hypothesis testing
[:ref:`8 <rdaref8>`]. He also argues that many of the constraining variables
used may end up being insignificant in their contribution to explaining
variation in samples. So should users be able to specify all environmental
variables? A large number of them? Or just one?

The second issue that might be of interest is using a second distance matrix
containing spatial distances to "partial out" the effect of spatial distance
before running DB-RDA, or to use the spatial distances as a constraining
environment variable. A series of posts details how to accomplish this using
DB-RDA [:ref:`9 <rdaref9>`]. It basically consists of running a Principal
Coordinates of Neighborhood Matrix (PCNM) over the spatial distance matrix and
using the output of that as an environment variable input to DB-RDA (can either
be a constraining variable or a partialling-out variable).

Existing Implementations
------------------------
There are existing implementations of RDA in the following statistical packages:

* XLSTAT [:ref:`1 <rdaref1>`]

* vegan package for R

* calibrate package for R

* open source R implementation by researcher [:ref:`5 <rdaref5>`]

* possibly others...

XLSTAT must be purchased and is only available on Windows and Mac OSX. The
implementations in R seems to be our best bet because it is open source and
people are already familiar with using R. There are three implementations that
I've found in R so far: vegan::rda, vegan::capscale, and calibrate::rda.

vegan::rda allows you to do partial and/or constrained RDA, while calibrate::rda
forces you to do constrained RDA. vegan::capscale is an implementation of
DB-RDA, in that it can accept a community data matrix or a distance matrix. It
allows you to decide whether you want to correct for negative eigenvalues or not
(see the discussion in the introduction section for more details on this
dispute).

I wrote a quick R script to demo vegan's capscale (DB-RDA) function on a QIIME
distance matrix. I used a series of posts [:ref:`7 <rdaref7>`] to a mailing list
as a guide for how to use capscale.

The script accepts a mapping file and a single category from the mapping file
that will be used as the constraining environmental variable (this can be
continuous, discrete, or categorical). The category will be used by capscale to
determine how much of the variability can be attributed to it. The script has
been checked into the Qiimeutils repository under :file:`microbiogeo/r/rda.r`.
The following sections of the document will explain how to set up your system to
run the script.

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
version of QIIME in SVN). The RDA script uses some R utility functions in QIIME
to load data.

Next, you must define an environment variable to tell the RDA script where to
look for the R utility functions in QIIME. Run the following command, changing
the path to point to the location of your QIIME install: ::

    export QIIME_DIR=/home/jrideout/qiime/trunk

If you don't want to have to perform this step each time you open a new
terminal, run the following command to add it to your .bashrc: ::

    echo "export QIIME_DIR=/home/jrideout/qiime/trunk" >> ~/.bashrc
    source ~/.bashrc

Next, run the following command to test if you can run the RDA script: ::

    R --slave --args -h < rda.r

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The RDA script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix and mapping file from the QIIME overview tutorial. You can get
the distance matrix
:download:`here <../downloads/overview_unweighted_unifrac_dm.txt>` and the
mapping file :download:`here <../downloads/Fasting_Map.txt>`.

Next, run the following command to execute the RDA script: ::

    R --slave --args -d overview_unweighted_unifrac_dm.txt -m Fasting_Map.txt -c Treatment < r/rda.r

Output Files
------------
The command in the previous section creates two output files named
:file:`rda_plot.pdf` and :file:`rda_results.txt`. The first file contains a 2D
plot of each of the samples. It seems very similar to the clustering shown by a
PCoA plot. The factor "Fast" overlayed on the plot is accompanied with a vector
showing what constraining factor grouped the fasting samples together. The other
output file contains information about the DB-RDA results. Notice that the
"Treatment" category accounts for 24.7% of the variability in the samples (this
information is found in the "Constrained" row of the results table).

Testing Results
---------------
This section will describe different tests that were run on the RDA script.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `BODY_SITE` category as a positive control. We expect to see
grouping of samples based on body site in the resulting plot.

**Command:** ::

    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c BODY_SITE < r/rda.r

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         159.1762                
    Real Total    165.4413     1.0000     
    Constrained    46.0873     0.2786   19
    Unconstrained 119.3540     0.7214  371
    Imaginary      -6.2651             213
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
        CAP1     CAP2     CAP3     CAP4     CAP5     CAP6     CAP7     CAP8 
    14.72239 10.95891  8.89776  3.26489  2.89957  1.41151  0.87627  0.69475 
        CAP9    CAP10    CAP11    CAP12    CAP13    CAP14    CAP15    CAP16 
     0.40960  0.35446  0.29999  0.24395  0.20137  0.18342  0.17567  0.15110 
       CAP17    CAP18    CAP19 
     0.13347  0.11498  0.09327 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    12.480  5.688  4.495  3.722  3.331  2.814  2.279  2.153 
    (Showed only 8 of all 371 unconstrained eigenvalues)

.. image:: ../images/rda/whole_body_test_1.png
   :align: center

The plot shows clear grouping of fecal samples at the bottom right of the
plot. Grouping of tongue samples can also be seen at the top right, and there
is also noticable grouping of outer ear samples at the bottom left. The plot
also contains overlayed vectors indicating which body sites explain the
grouping (not sure how better to explain this).

The output text shows that the `BODY_SITE` constraining variable explains
27.86% of the variability in the samples. These results seem to fall in line
with previous results seen in PCoA plots.

Test 2
~~~~~~
**Description:**

This test uses the `SEX` category as a negative control. We don't expect to see
grouping of samples due to previous analysis done on the Whole Body dataset.

**Command:** ::

    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c SEX < r/rda.r

**Results:**

The following output file is created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                     Inertia Proportion Rank
    Total         159.176211                
    Real Total    165.441288   1.000000     
    Constrained     1.146286   0.006929    1
    Unconstrained 164.295002   0.993071  371
    Imaginary      -6.265078             213
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
     CAP1 
    1.146 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    22.935 16.207 12.165  6.875  4.970  4.167  2.915  2.809 
    (Showed only 8 of all 371 unconstrained eigenvalues)

.. image:: ../images/rda/whole_body_test_2.png
   :align: center

The plot doesn't really show grouping of samples based on sex. The output text
shows that the `SEX` constraining variable explains only 0.6929% of the
variability in the samples. These results are what we'd expect.

Test 3
~~~~~~
**Description:**

This test uses the `BODY_SITE` category with three shuffled distance matrices.
We do not expect to see grouping of samples based on body site in the resulting
plots.

**Command:** ::

    compare_categories.py --method rda -i datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m datasets/whole_body/map.txt -c BODY_SITE -o rda_output
    compare_categories.py --method rda -i datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m datasets/whole_body/map.txt -c BODY_SITE -o rda_output
    compare_categories.py --method rda -i datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m datasets/whole_body/map.txt -c BODY_SITE -o rda_output

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                    Inertia Proportion Rank
    Total         159.17621                
    Real Total    165.44129    1.00000     
    Constrained     6.01504    0.03636   19
    Unconstrained 159.42624    0.96364  371
    Imaginary      -6.26508             213
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    1.25438 0.75405 0.62839 0.43759 0.40933 0.30193 0.27260 0.22559 0.22264 0.20888 
      CAP11   CAP12   CAP13   CAP14   CAP15   CAP16   CAP17   CAP18   CAP19 
    0.18437 0.17722 0.16618 0.16073 0.14860 0.13303 0.12340 0.11306 0.09308 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    22.109 15.686 11.616  6.796  4.856  4.118  2.868  2.817 
    (Showed only 8 of all 371 unconstrained eigenvalues)

.. image:: ../images/rda/whole_body_test_3_1.png
   :align: center
   :height: 600px
   :width: 600px

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                    Inertia Proportion Rank
    Total         159.17621                
    Real Total    165.44129    1.00000     
    Constrained     5.38679    0.03256   19
    Unconstrained 160.05450    0.96744  371
    Imaginary      -6.26508             213
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    1.38271 0.54748 0.47638 0.36795 0.32360 0.27671 0.26245 0.22565 0.19745 0.18654 
      CAP11   CAP12   CAP13   CAP14   CAP15   CAP16   CAP17   CAP18   CAP19 
    0.17411 0.16812 0.14466 0.13237 0.12705 0.10959 0.10331 0.09799 0.08266 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    21.897 15.966 11.812  6.790  4.881  4.176  2.897  2.835 
    (Showed only 8 of all 371 unconstrained eigenvalues)

.. image:: ../images/rda/whole_body_test_3_2.png
   :align: center
   :height: 600px
   :width: 600px

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                    Inertia Proportion Rank
    Total         159.17621                
    Real Total    165.44129    1.00000     
    Constrained     4.91953    0.02974   19
    Unconstrained 160.52176    0.97026  371
    Imaginary      -6.26508             213
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    0.92824 0.71286 0.43950 0.30075 0.28313 0.23922 0.22787 0.20940 0.19005 0.18279 
      CAP11   CAP12   CAP13   CAP14   CAP15   CAP16   CAP17   CAP18   CAP19 
    0.17583 0.15601 0.15023 0.14634 0.13792 0.13397 0.11478 0.10314 0.08748 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    22.346 15.693 12.002  6.790  4.880  4.234  2.923  2.802 
    (Showed only 8 of all 371 unconstrained eigenvalues)

.. image:: ../images/rda/whole_body_test_3_3.png
   :align: center
   :height: 600px
   :width: 600px

There doesn't appear to be a large amount of variability explained by
`BODY_SITE` when shuffled distance matrices are used.

Keyboard
^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category as a positive control. We expect
to see grouping of samples based on individual in the resulting plot.

**Command:** ::

    compare_categories.py --method rda -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda_out

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         21.60003                
    Real Total    21.64008    1.00000     
    Constrained    7.18171    0.33187   10
    Unconstrained 14.45837    0.66813  104
    Imaginary     -0.04005               6
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    4.61852 1.12311 0.36352 0.24780 0.20336 0.17231 0.14435 0.11807 0.10825 0.08242 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    1.1336 0.7402 0.7189 0.5341 0.5118 0.4736 0.4140 0.3991 
    (Showed only 8 of all 104 unconstrained eigenvalues)

.. image:: ../images/rda/keyboard_test_1.png
   :align: center

The plot shows three clear groups of samples, where each group contains the
samples for an individual. The output text shows that the `HOST_SUBJECT_ID`
constraining variable explains 33.19% of the variability in the samples. These
results seem to fall in line with previous results seen in PCoA plots.

Test 2
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category with three shuffled distance
matrices. We do not expect to see grouping of samples based on individual in the
resulting plots.

**Command:** ::

    compare_categories.py --method rda -i datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda_out
    compare_categories.py --method rda -i datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda_out
    compare_categories.py --method rda -i datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda_out

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         21.60003                
    Real Total    21.64008    1.00000     
    Constrained    1.97257    0.09115   10
    Unconstrained 19.66751    0.90885  104
    Imaginary     -0.04005               6
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    0.73504 0.22784 0.19892 0.16400 0.14474 0.12638 0.11865 0.09814 0.09152 0.06736 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    4.3752 1.3406 1.0699 0.6430 0.5843 0.4855 0.4559 0.4495 
    (Showed only 8 of all 104 unconstrained eigenvalues)

.. image:: ../images/rda/keyboard_test_2_1.png
   :align: center

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         21.60003                
    Real Total    21.64008    1.00000     
    Constrained    1.96717    0.09090   10
    Unconstrained 19.67291    0.90910  104
    Imaginary     -0.04005               6
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    0.75724 0.23330 0.21369 0.15925 0.13431 0.11837 0.11181 0.09178 0.08875 0.05866 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    4.3288 1.3103 1.0572 0.6717 0.5954 0.5139 0.4897 0.4379 
    (Showed only 8 of all 104 unconstrained eigenvalues)

.. image:: ../images/rda/keyboard_test_2_2.png
   :align: center

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         21.60003                
    Real Total    21.64008    1.00000     
    Constrained    1.82623    0.08439   10
    Unconstrained 19.81385    0.91561  104
    Imaginary     -0.04005               6
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
       CAP1    CAP2    CAP3    CAP4    CAP5    CAP6    CAP7    CAP8    CAP9   CAP10 
    0.43621 0.33219 0.20952 0.16133 0.14827 0.13586 0.12058 0.10970 0.09097 0.08159 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    4.7816 1.1089 1.0661 0.6794 0.5883 0.5182 0.4941 0.4513 
    (Showed only 8 of all 104 unconstrained eigenvalues)

.. image:: ../images/rda/keyboard_test_2_3.png
   :align: center

There doesn't appear to be a large amount of variability explained by
`HOST_SUBJECT_ID` when shuffled distance matrices are used, which is what we
would expect.

Glen Canyon
^^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category as a positive control. We expect
to see grouping of samples based on whether or not they are wet in the resulting
plot.

**Command:** ::

    compare_categories.py --method rda -i datasets/glen_canyon/unweighted_unifrac_dm.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o rda_out

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                  Inertia Proportion Rank
    Total         15.1717     1.0000     
    Constrained    3.5323     0.2328    1
    Unconstrained 11.6395     0.7672   92
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
     CAP1 
    3.532 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    1.4239 0.8082 0.4769 0.4267 0.3770 0.2911 0.2474 0.2418 
    (Showed only 8 of all 92 unconstrained eigenvalues)

.. image:: ../images/rda/glen_canyon_test_1.png
   :align: center

The plot shows two clear groups of samples, where one group contains the
samples that are currently wet and the other contains samples that are not. The
output text shows that the `CurrentlyWet` constraining variable explains 23.28%
of the variability in the samples. These results seem to fall in line with
previous results seen in PCoA plots.

Test 2
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category with three shuffled distance
matrices. We do not expect to see grouping of samples based on whether or not
they are wet in the resulting plots.

**Command:** ::

    compare_categories.py --method rda -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o rda_out
    compare_categories.py --method rda -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o rda_out
    compare_categories.py --method rda -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o rda_out

**Results:**

The following output files are created: ::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         15.17174    1.00000     
    Constrained    0.17157    0.01131    1
    Unconstrained 15.00018    0.98869   92
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
      CAP1 
    0.1716 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    3.7101 1.3615 0.7404 0.4472 0.4263 0.3606 0.2817 0.2445 
    (Showed only 8 of all 92 unconstrained eigenvalues)

.. image:: ../images/rda/glen_canyon_test_2_1.png
   :align: center

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                   Inertia Proportion Rank
    Total         15.17174    1.00000     
    Constrained    0.15236    0.01004    1
    Unconstrained 15.01938    0.98996   92
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
      CAP1 
    0.1524 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    3.6732 1.4065 0.7409 0.4493 0.4248 0.3709 0.2811 0.2463 
    (Showed only 8 of all 92 unconstrained eigenvalues)

.. image:: ../images/rda/glen_canyon_test_2_2.png
   :align: center

::

    Call: capscale(formula = as.dist(qiime.data$distmat) ~ factor, data =
    factors.frame)

                    Inertia Proportion Rank
    Total         15.171743   1.000000     
    Constrained    0.112111   0.007389    1
    Unconstrained 15.059632   0.992611   92
    Inertia is squared Unknown distance 

    Eigenvalues for constrained axes:
      CAP1 
    0.1121 

    Eigenvalues for unconstrained axes:
      MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    3.7049 1.4177 0.7392 0.4479 0.4232 0.3736 0.2812 0.2440 
    (Showed only 8 of all 92 unconstrained eigenvalues)

.. image:: ../images/rda/glen_canyon_test_2_3.png
   :align: center

There doesn't appear to be a large amount of variability explained by
`CurrentlyWet` when shuffled distance matrices are used, which is what we would
expect. It is hard to tell from the plots themselves, but the textual output
verifies this.

References
----------
.. _rdaref1:

[1] http://www.xlstat.com/en/products-solutions/feature/redundancy-analysis-rda.html

.. _rdaref2:

[2] http://www.jstor.org/stable/pdfplus/2657192.pdf?acceptTC=true

.. _rdaref3:

[3] http://ordination.okstate.edu/glossary.htm#RDA

.. _rdaref4:

[4] R's help page for vegan::rda

.. _rdaref5:

[5] http://www.bio.umontreal.ca/legendre/indexEn.html#RFunctions

.. _rdaref6:

[6] http://www.esajournals.org/doi/abs/10.1890/0012-9658(2001)082%5B0290:FMMTCD%5D2.0.CO;2

.. _rdaref7:

[7] http://r.789695.n4.nabble.com/R-question-about-capscale-vegan-td812694.html

.. _rdaref8:

[8] http://cc.oulu.fi/~jarioksa/opetus/metodi/mmmbeam2.pdf

.. _rdaref9:

[9] http://r.789695.n4.nabble.com/partial-dbRDA-or-CCA-with-two-distance-objects-in-Vegan-td2548762.html

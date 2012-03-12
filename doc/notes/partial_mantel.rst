.. _partial_mantel:

======================================================
Partial Mantel Test
======================================================

Synopsis
------------
The partial Mantel statistic is used to estimate the correlation between two matrices, 
A and B, while controlling for the effect of C. For instance, you could compare the
weighted and unweighted unifrac for soil, while controling for a third variable matrix,
such as pH or some measure of euclidian distance.

Introduction
------------
The partial Mantel test is a first-order correlation analysis that utilizes 
three distance(dis-similarity) matrices. [:ref:`1 <partial_mantelref1>`] This
builds on the simple Mantel which is a procedure that tests the hypothesis
that distances between the objects within a given matrix are linearly independent
of the distances withing those same objects in a separate matrix. [:ref:`2 <partial_mantelref2>`] 
It builds  on the simple Mantel by adding a third "control" matrix.

The goal is to test the correlation between
matrices A and B while controlling the effect of a third matrix C, in order to remove spurious 
correlations

The partial Mantel test accepts three matrices as input, the rest of its parameters it shares 
with the simple Mantel test, including permutations, strata and method[:ref:`3 <partial_mantelref3>`].
Partial Mantel statistic uses partial correlation conditioned on the third matrix. Only the first matrix 
its permuted so that the correlation structure between second and first matrices is kept constant. 
The goal is to test the correlation between matrices A and B while controlling the effect of a third 
matrix C, in order to remove spurious correlations


Existing Implementations
------------------------
There were found two existing implementations of the partial mantel test:

* vegan package for R

* zt 

vegan::mantel allows you to find the the partial Mantel statistics (using mantel.partial())
as the partial matrix correlation between three distance(dissimilarity) matrices.

zt is an implementation of the partial Mantel in C, it is free to under the GPL
license and it is open source. [:ref:`4 <partial_mantelref4>`]

Because zt is implemented to accept matrix data formatted for just the lower-left reflection
of the distance matrix, it is not immediately compatible with QIIME distance matrices, so, 
for now at least, I think the R implementation will be the best way to test.

System Setup and Required Dependencies
--------------------------------------
:note: The instructions that follow have been tested only on Mac OS x 10.7.3, but should be backward compatible for most Intel based Macs.

The following used a local install of QIIME 1.4.0-dev. First, install R. A binary is available in a self install PKG `here <http://cran.r-project.org/bin/macosx/>`_.

The python version is included and it does run (directions below), but I haven't been able to compare it to
the R version(werid optparse issue).

Next, you will need to install the vegan package and the optparse package
simply execute the following commands : ::

    sudo R
    install.packages("vegan")
    install.packages("optparse")
    q()

The install will require you to choose a mirror location for the package; 
any mirror will do, though you should choose one nearest your location.
The process is automatic after this step.

On this system I'm running vegan 2.0-2 and optparse version 0.9.4

Next, run the following command in order to test if you can run the partial Mantel script: ::

    R --slave --args -h < pmantel.r

:note: R script is not yet complete.

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The Partial Mantel script requires three distance(dissimilarity) matrices as input, 
for instance a Unifrac distance matrix as output by beta_diversity.py, a Euclidean
distance matrix and a median distance matrix as a control.
:download:`here <../downloads/partial_mantel_ready_to_run.zip>`. To derive the matrices yourself
you will need the original data and utilities(not finished,) which you can :download:`here <../downloads/partial_mantel_generate_own_matrices.zip>`

Please skip to the R script at the bottom if you don't wish to generate the matrices.

With everything set up we need to generate our three matrices. We start with the unweighted 
unifrac and extract the keyboard data and the sample subjects M2, M3 and M9. This is our 
first matrix(derived from the unweighted uniifrac of the keyboard study): ::

  filter_distance_matrix.py -i unweighted_unifrac_dm.txt -o unweighted_unifrac_dm_keyboard_only_239.txt -m meta_analysis_keyboard_map.txt -s 'COMMON_NAME:keyboard;HOST_SUBJECT_ID:M2,M3,M9'

:note: The data used is from the Fierer et al Keyboard Study.

To get our second matrix we will use the the euclidean distance between the physical keys
on the keyboard. In order to get a proper matrix, I had to use the list of element identifiers
from the previous output matrix, and using data extracted from a keyboard image annotated with  
GraphClick [:ref:`6 <partial_mantelref6>`]. The utility(not finished) is run as so: ::

  ./get_euclidian_dist_matrix.py -i unweighted_unifrac_dm_keyboard_only_239.txt -o unweighted_euclidean_dm.txt

The third(control) matrix will be the median unifrac distance between the sample subjects. We
first run make_distance_boxplots.py provided by the QIIME package: ::

  make_distance_boxplots.py -m meta_analysis_keyboard_map.txt -d unweighted_unifrac_dm_keyboard_only_239.txt -o .  -f HOST_SUBJECT_ID --suppress_all_between --suppress_all_within --save_raw_data

The output includes a PDF of the box plots, but we are interested in the "HOST_SUBJECT_ID_Distances.xls"
file. This will be used to generate our final matrix. We next run a second utility(not finished) as so: ::

  ./get_median_dist_matrix.py -i HOST_SUBJECT_ID_Distances.xls -m unweighted_unifrac_dm_keyboard_only_239.txt -o unifrac_median_dm.txt

We now should have all three of the matrices we need: unweighted_unifrac_dm_keyboard_only_239.txt, 
unweighted_euclidean_dm.txt and unifrac_median_dm.txt.

Finally, run the following command and execute the partial Mantel script: ::

    R --slave --args -a unweighted_unifrac_dm_keyboard_only_239.txt -b unweighted_euclidean_dm.txt -c unifrac_median_dm.txt < r/pmantel.r

To run the python version (requires test.py and compare_distance_matrices.py, included in download): ::

  ./compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239.txt,unweighted_euclidean_dm.txt,unifrac_median_dm.txt -o mantel_out.txt -n 1000 -m partial_mantel

Output Files
------------

R Version: outputs to stdout

Python Version: mantel_out.txt - includes the statistic as computed from the three matrices.


Testing Results
----------------
The keyboard data (outlined above) was used to test the Python implementation of the partial Mantel test


Keyboard
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the Python version of the partial Mantel. The three distance matrices derived
previously will be used.

**Command:** ::

  ./compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239.txt,unweighted_euclidean_dm.txt,unifrac_median_dm.txt -o mantel_out.txt -n 9999 -m partial_mantel

**Results:**

The following output file is created: ::

  DM1	DM2	DM3	Number of entries	Partial Mantel p-value
  unweighted_unifrac_dm_keyboard_only_239.txt	unweighted_euclidean_dm.txt	unifrac_median_dm.txt	74	0.506

The results are a little difficult to explain. The r-value basically splits the difference at ~0.5, which isn't definitive 
in support or in contrast to the hypothesis. We'll see also that the R version does not agree with Python version. 
For this reason, subsequent testing will use only the R version as it is more extensively documented and utilized.

Test 2
~~~~~~
**Description:**

In this test we use the vegan implementation of the partial Mantel test. As before we are looking at the
three distance matrices for the keyboard data set as derived above.

**Command:** ::

  R --slave --args -a unweighted_unifrac_dm_keyboard_only_239.txt -b unweighted_euclidean_dm.txt -c unifrac_median_dm.txt < pmantel.r 

**Results:**

The following was output to stdout: ::

  Mantel statistic r: 0.05618 
        Significance: 0.0583 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0451 0.0590 0.0723 0.0878
  

Test 3
~~~~~~
**Description:**

Negative Control: shuffle unifrac distmat

In this test we use the vegan implementation of the partial Mantel test.

**Command:** ::

  R --slave --args -a unweighted_unifrac_dm_keyboard_only_239_shuffled_1.txt -b unweighted_euclidean_dm.txt -c unifrac_median_dm.txt < pmantel.r

**Results:**

The following was output to stdout: ::

  Mantel statistic r: 0.01063 
        Significance: 0.3711 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0455 0.0584 0.0696 0.0847 


Test 4
~~~~~~
**Description:**

Negative Control: second shuffle of unifrac distmat

In this test we use the vegan implementation of the partial Mantel test.

**Command:** ::

  R --slave --args -a unweighted_unifrac_dm_keyboard_only_239_shuffled_2.txt -b unweighted_euclidean_dm.txt -c unifrac_median_dm.txt < pmantel.r

**Results:**

The following was output to stdout: ::

  Mantel statistic r: -0.02556 
        Significance: 0.768 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0455 0.0596 0.0722 0.0861


Test 5
~~~~~~~
**Description:**

Negative Control: third shuffle of unifrac distmat

In this test we use the vegan implementation of the partial Mantel test.

**Command:** ::
  
  R --slave --args -a unweighted_unifrac_dm_keyboard_only_239_shuffled_3.txt -b unweighted_euclidean_dm.txt -c unifrac_median_dm.txt < pmantel.r

**Results:**

The following was output to stdout: ::

  Mantel statistic r: -0.01656 
        Significance: 0.6841 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0449 0.0579 0.0708 0.0864 


88 Soils
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test compares the weighted and unwighted unifrac distances matrices for the 88 Soils
dataset. It uses a pH difference matrix as a third control matrix.

**Command:** ::

  R --slave --args -a unweighted_unifrac_dm.txt -b weighted_unifrac_dm.txt -c PH_dm.txt < pmantel.r  

**Results:**

The following output file is created: ::

  Mantel statistic r: 0.6818 
        Significance: 1e-04 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0449 0.0591 0.0710 0.0858 

Test 1
~~~~~~
**Description:**

In this test we compare the same PH_dm and weighted unifrac against three
unweighted unifrac matrices whose labels have been shuffled

**Command:** ::

  R --slave --args -a unweighted_unifrac_dm_shuffled_1.txt -b weighted_unifrac_dm.txt -c PH_dm.txt < pmantel.r  
  R --slave --args -a unweighted_unifrac_dm_shuffled_2.txt -b weighted_unifrac_dm.txt -c PH_dm.txt < pmantel.r  
  R --slave --args -a unweighted_unifrac_dm_shuffled_3.txt -b weighted_unifrac_dm.txt -c PH_dm.txt < pmantel.r  

**Results:**

The following was output, respective of the calls just above: ::

  Mantel statistic r: 0.02612 
        Significance: 0.2247 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0452 0.0588 0.0709 0.0836 

  Mantel statistic r: 0.06232 
        Significance: 0.0415 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0455 0.0590 0.0712 0.0844 

  Mantel statistic r: -0.009716 
        Significance: 0.5971 

  Empirical upper confidence limits of r:
     90%    95%  97.5%    99% 
  0.0443 0.0576 0.0686 0.0833 

We can see, in all three cases the the significance has changed dramatically, 100 fold in the least 
severe case (shuffled_2). This is in line with the expectation for the shuffled matrices.

References
---------- .. _partial_mantelref1:

[1] http://www.jstor.org/stable/2413122

.. _partial_mantelref2:

[2] http://www.bio.umontreal.ca/legendre/reprints/Partial_Mantel_paper.pdf

.. _partial_mantelref3:

[3] http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/mantel.html

.. _partial_mantelref4:

[4] http://www.jstatsoft.org/v07/i10/

.. _partial_mantelref5:

[5] http://www.bio.umontreal.ca/legendre/indexEn.html#RFunctions

.. _partial_mantelref6:

[6] http://www.arizona-software.ch/graphclick/


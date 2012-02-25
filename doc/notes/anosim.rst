==============================================================
Analysis of Similarities (ANOSIM) Statistical Method Reference
==============================================================

Introduction
------------
The Analysis of Similarities (ANOSIM) method tests for clusters by computing an
R statistic, which is based on mean ranks of within group and between group
dissimilarities, scaled into the range −1...+1. It is a non-parametric variant
of Analysis of Variance (ANOVA). This method operates on a distance matrix using
a given grouping of samples.

:note: ANOSIM only works with a categorical variable that is used to do the grouping. Mantel is recommended for continuous variables.

The R statistic that is calculated by ANOSIM is determined by the following
formula:

R = (rb-rw)/(N(N-1)/4)

where rb is the mean rank of all distances between groups and rw is the mean
rank of all distances within groups. An R value near +1 means
that there is dissimilarity between the groups [:ref:`2 <anosimref2>`].

The original paper referencing ANOSIM is a paper by K.R. Clarke
[:ref:`3 <anosimref3>`]. Here is a helpful description taken from
[:ref:`1 <anosimref1>`]:

The ANOSIM (Analysis of similarities) permutation method allows for testing for
group structure in the observations. If the original data contains abundance of
M species measured at N sites, an N – by – N matrix of (dis-)similarities D is
calculated. Suppose that the observations are from four different transects. As
a result, the matrix D contains a block structure:

.. image:: ../images/anosim.jpg
   :align: center
 
The sub-matrices Dii represent the (dis-)similarities between observations of
the same transect, and Dij between observations of different transects. A
statistic based on both the between and within sub-matrices is used to test the
differences among the 4 groups. Further details can be found in Legendre and
Legendre (1998), or in Chapter 10 of Zuur et al. (2007). A p-value for the
statistic is obtained by permutation. ANOSIM can be applied on 1-way data,
2-way nested data, 2-way crossed data with replication and 2-way crossed data
with no replication.

Existing Implementations
------------------------
There are several existing implementations of ANOSIM in statistical packages
that include:

* vegan package for R

* PERMANOVA add-in for PRIMER

* Fathom Toolbox for Matlab

* Brodgar

* possibly others...

ANOSIM has already been implemented in Python by Andrew Cochran but has not yet
been added to QIIME. The implementation has been checked into the Qiimeutils
repository under :file:`microbiogeo/python/`. The following sections of the
document will explain how to run Andrew's implementation of ANOSIM.

System Setup and Required Dependencies
--------------------------------------
:note: The following instructions have been tested on 64-bit Linux Mint (essentially Debian) using Python 2.6.7. However, they `should` work across different linux distros and on Macs. The instructions assume you use bash as your shell.

First, your system must have a version of QIIME installed (I used the latest
version of QIIME in SVN). The code also depends on NumPy, though this should
already be installed if you have QIIME installed. Next, you must add the area
where Andrew's code resides to your PYTHONPATH environment variable, changing
the path to point to the location of the microbiogeo checkout on your machine. I
also added the scripts area to my PATH for convenience: ::

    export PYTHONPATH=/home/jrideout/qiime/qiimeutils/microbiogeo:$PYTHONPATH
    export PATH=/home/jrideout/qiime/qiimeutils/microbiogeo/python/scripts:$PATH

If you don't want to have to perform this step each time you open a new
terminal, run the following command to add the path to your .bashrc: ::

    echo "export PYTHONPATH=/home/jrideout/qiime/qiimeutils/microbiogeo:$PYTHONPATH" >> ~/.bashrc
    echo "export PATH=/home/jrideout/qiime/qiimeutils/microbiogeo/python/scripts:$PATH" >> ~/.bashrc
    source ~/.bashrc

Next, run the following command to test if you can run the ANOSIM script: ::

    anosim.py -h

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The ANOSIM script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix from the QIIME overview tutorial. You can get the distance
matrix :download:`here <../downloads/overview_unweighted_unifrac_dm.txt>` and
the mapping file :download:`here <../downloads/Fasting_Map.txt>`.

Next, run the following command to execute the ANOSIM script: ::

    anosim.py -i overview_unweighted_unifrac_dm.txt -m Fasting_Map.txt -c Treatment -o anosim_results.txt

The -c option specifies which column in the mapping file will be used to group
the samples. The `Treatment` column has two values: 'Control' and 'Fast'. Thus,
ANOSIM will be used to calculate the dissimilarity between the control and fast
groups. The -o option specifies the file that we want the results written to.

Output Files
------------
The command in the previous section creates a single output file named
:file:`anosim_results.txt`. The resulting file should look like this: ::

    Input_filepath  ANOSIM_R_value  p_value
    overview_unweighted_unifrac_dm.txt      0.8125  NA

The first field lists the distance matrix file that was used as input. The
second field lists the R statistic that was computed (remember that this is the
primary output of ANOSIM). The final field lists the p-value, which is NA
because we did not specify the optional -p parameter (by default, the number of
p-trials is 0).

The value of the R statistic can fall between -1 and +1, with a positive value
close to 1 indicating that the groups are highly dissimilar. Thus, in this
example, the control and fast groups are dissimilar.

Testing Results
---------------
This section will describe different tests that were run on the ANOSIM script.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `BODY_SITE` category as a positive control. We expect there
to be significant clustering due to previous analysis done on the Whole Body
dataset.

**Command:** ::

    anosim.py -i datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c BODY_SITE -o anosim_results.txt -p 999

**Results:**

The following output file is created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/whole_body/unweighted_unifrac_dm.txt	0.469648075442	0.001

The R value of 0.469648075442 indicates that body sites are significantly
different (i.e. there is clustering) due to its relatively "large" positive
value. This is a result that we would expect. The p-value of 0.001 indicates
that the result is significant.

Test 2
~~~~~~
**Description:**

This test uses the `SEX` category as a negative control. We don't expect to see
significant clustering due to previous analysis done on the Whole Body dataset.

**Command:** ::

    anosim.py -i datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c SEX -o anosim_results.txt -p 999

**Results:**

The following output file is created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/whole_body/unweighted_unifrac_dm.txt	0.0354433583741	0.002

The R value of 0.0354433583741 indicates that there isn't significant clustering
due to sex of the subjects because it is close to zero. This result is what we
would expect. The only confusing thing is the p-value of 0.002. This is a really
small p-value, so it **is** indicating that there are significant differences
between the groups.

Test 3
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `BODY_SITE` category to
perform three negative control tests. Since the labels of the distance matrices
are shuffled, we don't expect to see clustering any more on this category.

**Command:** ::

    anosim.py -i datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m datasets/whole_body/map.txt -c BODY_SITE -o anosim_results.txt -p 999
    anosim.py -i datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m datasets/whole_body/map.txt -c BODY_SITE -o anosim_results.txt -p 999
    anosim.py -i datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m datasets/whole_body/map.txt -c BODY_SITE -o anosim_results.txt -p 999

**Results:**

The following output files are created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt	-0.0085666370674	0.771

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt	-0.00260471465844	0.571

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt	-0.00382322857638	0.632

The R values of -0.0085666370674, -0.00260471465844, and -0.00382322857638
indicate that body sites are no longer significantly different once the distance
matrices are shuffled, which is what we would expect.

Keyboard
^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category as a positive control. We expect
there to be significant clustering on host subjects due to previous analysis
done on the keyboard study dataset.

**Command:** ::

    anosim.py -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim_results.txt -p 999

**Results:**

The following output file is created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/keyboard/unweighted_unifrac_dm.txt	0.794026410205	0.001

The R value of 0.794026410205 indicates that samples taken from different hosts
are significantly different (i.e. there is clustering) due to its "large"
positive value. This is a result that we would expect. The p-value of 0.001
indicates that the result is significant.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `HOST_SUBJECT_ID`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    anosim.py -i datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim_results.txt -p 999
    anosim.py -i datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim_results.txt -p 999
    anosim.py -i datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim_results.txt -p 999

**Results:**

The following output files are created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt	-0.00712796151372	0.6

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt	0.00843082850421	0.342

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt	-0.00611883437807	0.59

The R values of -0.00712796151372, 0.00843082850421, and -0.00611883437807
indicate that samples taken from different host subjects are no longer
significantly different once the distance matrices are shuffled, which is what
we would expect.

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category as a positive control. We expect
there to be significant clustering on this category due to previous analysis
done on the Glen Canyon dataset.

**Command:** ::

    anosim.py -i datasets/glen_canyon/unweighted_unifrac_dm.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o anosim_results.txt -p 999

**Results:**

The following output file is created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/glen_canyon/unweighted_unifrac_dm.txt	0.9984007035	0.001

The R value of 0.9984007035 indicates that samples taken from wet and dry
environments are significantly different (i.e. there is clustering) due to the
really "large" positive value that is close to 1. This is a result that we would
expect, as there is also clear clustering in the 3D PCoA plots.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `CurrentlyWet`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    anosim.py -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o anosim_results.txt -p 999
    anosim.py -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o anosim_results.txt -p 999
    anosim.py -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o anosim_results.txt -p 999

**Results:**

The following output files are created: ::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt	0.0876180335381	0.129

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt	0.0074529653733	0.415

::

    Input_filepath	ANOSIM_R_value	p_value
    datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt	-0.0507653473398	0.711

The R values of 0.0876180335381, 0.0074529653733, and -0.0507653473398 indicate
that samples taken from wet vs. dry environments are no longer significantly
different once the distance matrices are shuffled, which is what we would
expect.

References
----------
.. _anosimref1:

[1] http://www.brodgar.com/manual/Chapter6BMS.pdf

.. _anosimref2:

[2] http://folk.uio.no/ohammer/past/multivar.html

.. _anosimref3:

[3] Clarke, K.R. 1993. Non-parametric multivariate analysis of changes in community structure. Australian Journal of Ecology 18:117-143.

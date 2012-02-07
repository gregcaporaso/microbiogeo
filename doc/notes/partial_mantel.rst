.. _partial_mantel:

======================================================
Partial Mantel Test
======================================================

Introduction
------------
The partial Mantel test is a first-order correlation analysis that utilizes 
three distance(dis-similarity) matrices. [:ref:`1 <partial_mantelref1>`] This
builds on the simple Mantel which is a procedure that tests the hypothesis
that distances between the objects within a given matrix are linearly independent
of the distances withing those same objects in a separate matrix. [:ref:`2 <partial_mantelref2>`] 
It builds  on the simple Mantel by adding a third "control" matrix.

The partial Mantel test accepts three matrices as input, the rest of its parameters it shares 
with the simple Mantel test, including permutations, strata and method. [:ref:`3 <partial_mantelref3`] 

Existing Implementations
------------------------
There were found two existing implementations of the partial mantel test:

* vegan package for R

* zt 

vegan::mantel allows you to find the the partial Mantel statistics (using mantel.partial())
as the partial matrix correlation between three distance(dissimilarity) matrices.

zt is an implementation of the partial Mantel in C, it is free to under the GPL
license and it is open source. [:ref:`4 <partial_mantelref4`]

Because zt is implemented to accept matrix data formatted for just the lower-left reflection
of the distance matrix, it is not immediately compatible with QIIME distance matrices, so, 
for now at least, I think the R implementation will be the best way to test.

System Setup and Required Dependencies
--------------------------------------
:note: The instructions that follow have been tested only on Mac OS x 10.7.2, but should 
be backward compatible for most Intel based Macs.

First, install R. A binary is available in a self install PKG `here <http://cran.r-project.org/bin/macosx/>`

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

Next, run the following command to test if you can run the RDA script: ::

    R --slave --args -h < pmantel.r

:note: R scrpt is not yet complete.
This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The Partial Mantel script requires three distance(dissimilarity) matrices as input, 
for instance a Unifrac distance matrix as output by beta_diversity.py, a Euclidean
distance matrix and a mean distance matrix as a control.
:download:`here <../downloads/not_yet_available.txt>`.

Next, run the following command to execute the RDA script: ::

    R --slave --args -d1 matrix1.txt -d2 matrix2.txt -d3 control_matrix.txt < r/pmantel.r

Output Files
------------
TBD, have not yet run succesfully

References
----------
.. _partial_mantelref1:

[1] http://www.jstor.org/stable/2413122

.. _partial_mantelref2:

[2] http://www.bio.umontreal.ca/legendre/reprints/Partial_Mantel_paper.pdf

.. _partial_mantelref3:

[3] http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/mantel.html

.. _partial_mantelref4:

[4] http://www.jstatsoft.org/v07/i10/

.. _partial_mantelref5:

[5] http://www.bio.umontreal.ca/legendre/indexEn.html#RFunctions


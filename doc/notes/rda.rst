.. _rda:

======================================================
Redundancy Analysis (RDA) Statistical Method Reference
======================================================

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
interesting approach in what they call 'distance-based RDA', or 'db-RDA' for
short. This method basically takes a distance matrix of any type, performs
principal coordinates analysis (PCoA) on it, corrects for negative eigenvalues,
and then performs RDA on the result and another matrix of environmental dummy
variables to analyze their relationship. This application of RDA might be the
most relevant to QIIME/microbial ecology, so it may be worth looking into.

Another paper challenges the method used in the previously described paper
[:ref:`6 <rdaref6>`]. The authors argue that the negative eigenvalue correction
step is not necessary.

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
db-RDA, in that it can accept a community data matrix or a distance matrix. It
allows you to decide whether you want to correct for negative eigenvalues or not
(see the discussion in the introduction section for more details on this
dispute).

I wrote a quick R script to demo vegan's RDA on a QIIME distance matrix. It
currently does not accept input for Y and Z matrices. It has been checked into
the Qiimeutils repository under :file:`microbiogeo/r/rda.r`. The following
sections of the document will explain how to set up your system to run the demo.

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
beta_diversity.py). I used the unweighted Unifrac distance matrix from the QIIME
overview tutorial. You can get the distance matrix
:download:`here <../downloads/overview_unweighted_unifrac_dm.txt>`.

Next, run the following command to execute the RDA script: ::

    R --slave --args -d overview_unweighted_unifrac_dm.txt < r/rda.r

Output Files
------------
The command in the previous section creates a single output file named
:file:`rda_plot.pdf`. This file contains a 2D plot of each of the samples. It
seems that the sample IDs colored black represent a PCoA plot, while the red
colored sample IDs are the RDA plot. From comparison with the 2D PCoA plot
generatd by QIIME's make_2D_plots.py, it appears the same clustering exists for
RDA that we see in PCoA.

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

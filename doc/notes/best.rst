.. _primere_best:

======================================================
Primer-e BEST
======================================================

Introduction
------------

BEST [:ref:`1 <primere_bestref1>`] Selects environmental variables, or species  
"best explaining" community pattern, by maximising a rank correlation between 
their respective resemblance matrices. Two algorithms are available. In the 
BIOENV algorithm  all permutations of the trial variables are tried. In the BVSTEP 
algorithm a stepwise search over the trial variables is tried. Use BVSTEP if 
there is a large number of trial variables and BIOENV is too slow.

The backbone of the of the BEST analyses are the two closely related algorithms, 
BIOENV and the stepwise version of BIOENV, BVSTEP. 

BIOENV analysis allows for the comparison of distance/similarity matrices between 
two sets of data having either samples or variables in common. 

BVSTEP is intended for a faster exploration of the subset combinations. It looks 
specifically at the BIOBIO type of exploration and addresses the concept of structural 
redundancy in composition of communities through the identification of "response units" 
or Taxonomic/functional groupings of species that changed in abundance in the same way 
over time.

Although the R vegan package has an implementation of the BIOENV algorithm,  

Existing Implementations
------------------------
There were found three existing implementations of the "BEST" test:

* Primer-e 

* Vegan BIOENV 

* Unnamed R implementation of the BIOENV and BVSTEP [:ref:`2 <primere_bestref2>`] 

Primer-e is a proprietary implementation of an undocumented (under the name BEST)
statistic.


System Setup and Required Dependencies
--------------------------------------
Windows XP+ and the Primer-e software is required to run the Primer-e version of 
the software.

The R versions require R:

:note: The instructions that follow have been tested only on Mac OS x 10.7.3, but should be backward compatible for most Intel based Macs.

The following used a local install of QIIME 1.4.0-dev. First, install R. A binary is available in a self install PKG `here <http://cran.r-project.org/bin/macosx/>`.

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

To run the BIOENV/BVSTEP algorithms in R, the required scripts are included 
(bioenv.r and bvstep.r) To run them simply import the functions into R from 
the R shell: ::

  > source('r/bvstep.r')
  > source('r/bioenv.r')

To test that this worked correctly, pull some sample data from the Vegan package: ::
  
  > library('vegan')
  > data(varespec)
  > data(varechem)

Then, to run the R functions we just sources: ::
  
  > res <- bio.env(varespec, varechem,  fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 

And to run the Vegan BIOENV implementation: ::

  > res <- bioenv(comm = varespec, env = varechem)

And you can summarize the results: ::

  > summary(res)

:note:Method still needs to be tested with empirical data form a collection of datasets.


Input Files
-----------

Sample data worksheet with variables tested to explain the community (active worksheet).
Resemblance matrix for the community to be explained (selected from drop down list). 
The actual samples considered are as selected in the active data worksheet.  Label 
matching between the worksheets is automatic but the resemblance matrix selection must 
contain at least these sample names.  They can be in any order. Data worksheet must 
have no missing values.


Output Files
------------

Outputs directly to R's stdout.

References
----------
.. _primere_bestref1:

[1] http://www.primer-e.com/

.. _primere_bestref2:

[2] http://menugget.blogspot.com/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

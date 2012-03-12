.. _primere_best:

======================================================
BEST (BIOENV and BVSTEP, via Primer-E)
======================================================

Introduction
------------
the BIOENV and BVSTEP algorithms look at the numerical environmental variables relating samples 
in a distance matrix. For instance, the unifrac distance matrix and pH and latitude (or any other 
number of variables) in soil samples, and ranks them in order of which best explain patterns in 
the communities.


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
two sets of data having either samples or variables in common. BIOENV uses all the 
available environmental variables to find the combination that ‘best explains’ the 
patterns in the biological data. Basically, BIOENV seems to take a distance matrix,
i.e. unifrac distance matrix and then a columnization of variables (pH, location, etc.)


However, when large numbers (>15 or 16) of environmental variables are used the procedure 
can become impractically slow, and computation time may be excessive. In such cases the BVSTEP 
option can be employed to carry out a stepwise search of the variables, employing both 
forward selection and backward elimination. Starting with the variable showing the maximum 
matching coefficient, variables are successively added, the combinations tested and (at each stage) 
the variable contributing least eliminated. BVSTEP, could, potentially be used to compare 
two distance matrices, though, the value in this is unclear, since the idea behind the analysis
is to compare multiple variables 

Existing Implementations
------------------------
There were found three existing implementations of the "BEST" test:

* Primer-e 

* Vegan BIOENV 

* Unnamed R implementation of the BIOENV and BVSTEP [:ref:`2 <primere_bestref2>`] 

Primer-e is a proprietary implementation of a (mostly) undocumented (under the name BEST)
statistic. So, this implementation may be effective for testing purposes, though, its difficult to 
determine the appropriate parameterization; there is little to no documentation.

The nameless R implementation is not convincingly discussed and does not explain its 
methodology in a sufficient detail.

The Vegan version of BIOENV is fairly straightforward to work with and the results seem reasonable.
For the time being this method will be used for testing, however, it would be desirable to generate
comparable results in Primer-e as a verification of the validity of the output.


System Setup and Required Dependencies
--------------------------------------
Windows XP+ and the Primer-e software is required to run the Primer-e version of 
the software.

The R versions require R:

:note: The instructions that follow have been tested only on Mac OS x 10.7.3, but should be backward compatible for most Intel based Macs.

The following used a local install of QIIME 1.4.0-dev. First, install R. A binary is available in a self install PKG `here <http://cran.r-project.org/bin/macosx/>`_.

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

:note: Method still needs to be tested with empirical data from a collection of datasets.


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


Testing Results
----------------
The hope is to use the vegan implementation and compare against the Primer-e.
So far, only the vegan has be run satisfactorily.

88 Soils
^^^^^^^^^^^^^^^
For these initial tests of the BIOENV algorithm I chose the 88 Soils data set.
This particular dataset was chosen for its multi-variate data in the mapping file.
Additionally, the BEST analysis seems suited for environmental data.

Test 0 (Vegan BIOENV)
~~~~~~~~~~~~~~~~~~~~~~
**Description:**

The version of BIOENV built-in to the Vegan package is tested here. Ideally 
this will line up with BEST test case to follow. However, it seems skewed toward
the vegan enironment datasets (i.e., varespec and varechem.)

:note: The row count must match in order to perform the comparison. The 88 Soils data had additional sample rows in the mapping file which were not included in the unifrac distance matrix. These were removed and a new file generated *vars.txt*

The variables table can only include numerical data, so any non-numerical columns were removed. 

:note: In order to compare the Primer-e and Vegan, we had to use "euclidean" as the dissimilarity index, since this was the only one common to both that allows negative values in the variables. The Primer-e software also states the this distance metric is well suited to environmental data.

**Command:** ::

  R --slave --args -c unweighted_unifrac_dm.txt -e vars.txt < r/best.r > best_result.txt

**Results:**
:note: Output is likely invalid, I've just discovered that the vegan version runs a vegdist() function on the input matrix, I expect this is not leaving the data in a usable state. Attempting to re-write the method. However, this method won't be used further for now.

And the output is: ::

  Subset of environmental variables with best correlation to community data.

  Correlations:      spearman 
  Dissimilarities:   euclidean 

  Best model has 1 parameters (max. 11 allowed):
  PH
  with correlation  0.7764964 


Test 1 (Primer-e)
^^^^^^^^^^^^^^^^^^
**Description:**

The Primer-e version has significantly more(seemingly) configuration options.
This is a positive control, using the original, valid, distance matrix. And the
variables: TOT_ORG_CARB, SILT_CLAY, ELEVATION, SOIL_MOISTURE_DEFICIT, CARB_NITRO_RATIO, ANNUAL_SEASON_TEMP, ANNUAL_SEASON_PRECPT, PH, CMIN_RATE, LONGITUDE, LATITUDE


**Command:**

There is no command, per-se, all of the methods in Primer-e are run in the
Windows GUI which lays on top of the software. However, the following steps were 
take:

* Open the unweighted_unifrac_dm.txt

  * Use the open file and choose the .txt file
  * Select the "Resemblance matrix" option and click Next>
  * Uncheck the "Title check box", Select "Distance" and click Next>
  * Click "Finish"

* Open the vars.txt file

  * Choose "Sample data" and click Next>
  * Uncheck the "Title" checkbox
  * Click "Samples as rows" and click Next>
  * Click "Finish"

Now, with the vars.txt selected 

* Choose Analyse > BEST...
* Choose the BIOENV tab and set the value there to 15(11 is actually sufficient for this data)
* Choose the "General" tab
* Click "Resemblance..." and choose "Euclidean" then click "OK"
* Reopen the BEST analysis window and click the BVSTEP radio button. Click "OK"


**Results:**

At the bottom of the analysis window you should have for the BIOENV: ::

  Best results
  No.Vars    Corr. Selections
        1    0.738 8
        9    0.419 1,2,4-6,8-11
        8    0.419 1,2,4-6,9-11
        8    0.419 1,2,4-6,8,10,11
        7    0.419 1,2,4-6,10,11
        8    0.419 1,2,4,6,8-11
        7    0.419 1,2,4,6,9-11
        7    0.419 1,2,4,6,8,10,11
        6    0.419 1,2,4,6,10,11
        8    0.418 1,2,4,5,8-11
        7    0.418 1,2,4,5,9-11
        7    0.418 1,2,4,5,8,10,11
        6    0.418 1,2,4,5,10,11
        7    0.418 1,2,4,8-11
        6    0.418 1,2,4,9-11

And for the BVSTEP: ::

  Best results
  Multiple   No.Vars    Corr.    Selections
  1             1       0.738     8

Where the variables are numbered as such: ::

  1 TOT_ORG_CARB
  2 SILT_CLAY
  3 ELEVATION
  4 SOIL_MOISTURE_DEFICIT
  5 CARB_NITRO_RATIO
  6 ANNUAL_SEASON_TEMP
  7 ANNUAL_SEASON_PRECPT
  8 PH
  9 CMIN_RATE
  10 LONGITUDE
  11 LATITUDE

I believe that the discrepancy between these results and the Vegan
results are due to the fact that vegan forces vegdist() call on the input 
distance matrix. 

To my understanding a "high" correlation between the
dissimilarity matrix and the variable indicates that there is a good variance, for instance when
looking at the PH the statistic is high for the "dissimilarity."


Test 2
~~~~~~~

**Description:**

In this test we wanted to use three shuffled distance matrices (each shuffled matrix is derived 
from the original unweighted_unifrac_dm.txt"

**Command:**

The same procedures were followed as outlined in Test 1. Once for each shuffled matrix.
(unweighted_unifrac_dm_shuffled_1, unweighted_unifrac_dm_shuffled_2, unweighted_unifrac_dm_shuffled_3)

**Results:**

The result files were actually identical to Test1, I believe this is because the BEST analsis 
actually matches sample names in the distance matrix to sample names in the "mapping"/variables data.
I'm not exactly sure what would make a good negative control at this point.

References
----------
.. _primere_bestref1:

[1] http://www.primer-e.com/

.. _primere_bestref2:

[2] http://menugget.blogspot.com/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

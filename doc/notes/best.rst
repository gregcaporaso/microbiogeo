.. _primere_best:

======================================================
Primer-e BEST (BIOENV and BVSTEP)
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
two sets of data having either samples or variables in common. BIOENV uses all the 
available environmental variables to find the combination that ‘best explains’ the 
patterns in the biological data. Basically, BIOENV seems to take a distance matrix,
i.e. unifrac distance matrix and then a colomnization of variables (pH, location, etc.)
BIOENV 

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

Vegan BIOENV
^^^^^^^^^^^^^

To test the vegan implementation of the BIOENV algorithm I chose the 88 Soils data set.
This particular dataset was chosen for its multi-variate data in the mapping file.

:note: The row count must match in order to perform the comparison. The 88 Soils data had additional sample rows in the mapping file which were not included in the unifrac distance matrix. These were removed and a new file generated *vars.txt*

The variables table can only include numerical data, so any non-numerical columns were removed. 

to run the the vegan a small script has been written to assist, simply run: ::

  R --slave --args -c unweighted_unifrac_dm.txt -e vars.txt < r/best.r > best_result.txt

:note: Output is likely invalid, I've just discovered that the vegan version runs a vegdist() function on the input matrix, I expect this is not leaving the data in a usable state. Attempting to re-write the method.

And the output is: ::

  2047 possible subsets (this may take time...)

  Call:
  bioenv(comm = cdm, env = edm) 

  Subset of environmental variables with best correlation to community data.

  Correlations:      spearman 
  Dissimilarities:   bray 

  Best model has 1 parameters (max. 11 allowed):
  PH
  with correlation  0.7649623 

                                                                                                                                                  size
  PH                                                                                                                                                 1
  SOIL_MOISTURE_DEFICIT PH                                                                                                                           2
  SOIL_MOISTURE_DEFICIT ANNUAL_SEASON_PRECPT PH                                                                                                      3
  SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_PRECPT PH                                                                                     4
  SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH                                                                  5
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_PRECPT PH LONGITUDE                                                                 6
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH LONGITUDE                                              7
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE                                    8
  TOT_ORG_CARB SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE                       9
  TOT_ORG_CARB SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE LATITUDE             10
  TOT_ORG_CARB SILT_CLAY ELEVATION SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE LATITUDE   11
                                                                                                                                                  correlation
  PH                                                                                                                                                   0.7650
  SOIL_MOISTURE_DEFICIT PH                                                                                                                             0.7166
  SOIL_MOISTURE_DEFICIT ANNUAL_SEASON_PRECPT PH                                                                                                        0.6582
  SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_PRECPT PH                                                                                       0.5944
  SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH                                                                    0.5435
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_PRECPT PH LONGITUDE                                                                   0.5072
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH LONGITUDE                                                0.4773
  SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE                                      0.4468
  TOT_ORG_CARB SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE                         0.4135
  TOT_ORG_CARB SILT_CLAY SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE LATITUDE                0.3796
  TOT_ORG_CARB SILT_CLAY ELEVATION SOIL_MOISTURE_DEFICIT CARB_NITRO_RATIO ANNUAL_SEASON_TEMP ANNUAL_SEASON_PRECPT PH CMIN_RATE LONGITUDE LATITUDE      0.3479 

Primer-e BIOENV
^^^^^^^^^^^^^^^^

References
----------
.. _primere_bestref1:

[1] http://www.primer-e.com/

.. _primere_bestref2:

[2] http://menugget.blogspot.com/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

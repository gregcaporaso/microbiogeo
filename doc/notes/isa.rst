.. _isa:

======================================================
Indicator Species Analysis
======================================================

Synopsis
------------

The identification of representative(characteristic) or indicator species is a traditional 
problem in biogeography and mecroecology. Studies in the field which describe sites and/or habitats 
usually mention one or multiple species that are representative of each habitat, in other words, which 
characterize the habitat. Obviously then, there is a need for the identification of
characteristic/indicator species in other fields. Since represenative/indicator species add ecological 
meaning to groups of sites discovered through clustering, they provide criteria to

* compare typologies derived from data analysis,
* identify where to stop dividing clusters into subsets, and 
* point out the main levels in a hierarchical classification of sites.

Good indicator species should be found mostly in a single
group of a typology and be present at most of the sites belonging to that group.


Introduction
------------
Indicator species are species that can be used as ecological indicators of community/habitat types, environmental conditions, 
or environmental fluctuations and changes. In macroecology indicator species are plants and animals that, in their presence, 
abundance, or chemical makeup, demonstrate some representatively distinctive aspect of the character or quality of the 
environment(in which they are found.) As an example, in areas where metal-rich minerals can be found at the soil surface, 
indicator species of plants accumulate large concentrations of those minerals in their tissues. 

It is plausible, given that ISA has found application in macroecology, that it is extensible to microbial ecology as well. 
There are multitude of situations where the specie and characteristics of that specie can  be used to indicate the 
environment where it is found, perhaps certain drug-resistant pathogens or organisms which thrive in certain pH environments
could be used to characterize their environment instead of the environments characterizign the microbiota. These are not the best
examples, but it is indicative of the fact that the utility of ISA in microbial ecological identification is a very real possibility.

Existing Implementations
------------------------
There were found two existing implementations of the Indicator Species Analysis:

* The LabDSV package for R[:ref:`3 <isaref3>`]
* PC-ORD[:ref:`2 <isaref2>`]

The LabDSV package is an R package and the ISA functionality is found in the indval() function, formerly duleg().

PC-ORD is a Windows software package which offers multiple statistical tools, among them is [some] version of ISA,
but it is proprietary software and it is not freely available.

Because the PC-ORD software is not freely available and verifiable, we will use the R version in the LabDSV
R package for testing and analysis of the method's utility. Unfortunately the majority of the functionality 
is actually written in FORTRAN, however, it is functioning correctly.

System Setup and Required Dependencies
--------------------------------------
:note: The instructions that follow have been tested only on Mac OS x 10.7.3, but should be backward compatible for most Intel based Macs.

The following used a local install of QIIME 1.4.0-dev. 

First, install R. A binary is available in a self install PKG `here <http://cran.r-project.org/bin/macosx/>`_.

Next, you will need to install the LabDSV package and the optparse package
simply execute the following commands: ::

    sudo R
    install.packages("labdsv")
    install.packages("optparse")
    q()

The install will require you to choose a mirror location for the package; 
any mirror will do, though you should choose one nearest your location.
The process is automatic after this step.

On this system I'm running vegan 2.0-2 and optparse version 0.9.4

Next, run the following command in order to test if you can run the partial Mantel script: ::

    R --slave --args -h < isa.r

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------

Indicator Species Analysis takes a table/matrix where the columns are species and the samples are rows.
Additionally it takes a vector of numeric cluster memberships for samples.

Output Files
------------

The output is printed to stdout. It can be summarized, but generally outputs 


Testing Results
----------------
It is not immediately obvious what datasets will be most helpful. It seems the Glen Canyon and 88 Soils would
be good candidates since they would both 

Specifically, though, it would seem, that the OTU tables, and not distance matrices will be used here.
The only issue, there, is that the OTU tables are inherently clustered and so the "clustering" classes used by
the method would be redundant. Though, the OTU table is representative sequences, so it may actually be okay.


Example Run(Fake Data)
^^^^^^^^^^^^^^^^^^^^^^^
**Description:**

The following is just a sample run; it uses the R functionality directly, rather than through the isa.r script.

**Command:**

Setup: ::

  library('labdsv')
  m = matrix( c(4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,3,3,8,8,8,8,8,4,4,4,4,4,6,6,6,6,6,4, 4,2,0,0,0,0,0,0,0,18,18,18,18,18,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),25,3)
  vec = c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5))
  colnames(m) = c("species_1", "species_2", "species_3")

* The second line a 25x3 (row,col) matrix. 
* The third line creates a a vec, which will be the cluster memberships for the samples
* The last line is purely aesthetic and labels the columns for easy identification.

Run: ::

  indval(m, vec)

**Results:**

The following is ouput to stdout: ::

  $relfrq
            1 2 3   4 5
  species_1 1 1 1 1.0 1
  species_2 1 1 1 0.6 0
  species_3 1 1 0 0.0 0

  $relabu
              1    2    3    4    5
  species_1 0.2 0.25 0.25 0.15 0.15
  species_2 0.4 0.20 0.30 0.10 0.00
  species_3 0.9 0.10 0.00 0.00 0.00

  $indval
              1    2    3    4    5
  species_1 0.2 0.25 0.25 0.15 0.15
  species_2 0.4 0.20 0.30 0.06 0.00
  species_3 0.9 0.10 0.00 0.00 0.00

  $maxcls
  species_1 species_2 species_3 
          2         1         1 

  $indcls
  species_1 species_2 species_3 
       0.25      0.40      0.90 

  $pval
  species_1 species_2 species_3 
      0.027     0.001     0.001 
  

Dataset 1
^^^^^^^^^^
Test 1
~~~~~~
**Description:**


**Command:** ::

  Input

**Results:**

The following output file is created: ::

  Output

Test 2
~~~~~~~
**Description:**

**Command:** ::

  Input

**Results:**

The following was output to stdout: ::

  Output

Test 3
~~~~~~~
**Description:**

**Command:** ::

  Output

**Results:**

Test 4
~~~~~~
**Description:**

**Command:** ::

  Input

**Results:**

The following was output to stdout: ::

  Output


Test 5
~~~~~~~
**Description:**

**Command:** ::

  Input
  
**Results:**

The following was output to stdout: ::

  Output


Dataset 2
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

**Command:** ::

  Input

**Results:**

The following output file is created: ::

  Output

Test 1
~~~~~~
**Description:**

**Command:** ::

  Input

**Results:**

The following was output to stdout: ::

  Output


References
-----------
.. _isaref1:

[1] http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2121141/ 

.. _isaref2:

[2] http://people.oregonstate.edu/~mccuneb/pcord.htm

.. _isaref3:

[3] http://rss.acs.unt.edu/Rdoc/library/labdsv/html/duleg.html

.. _isaref4:

[4] http://ecology.msu.montana.edu/labdsv/R/labs/

.. _isaref5:

[5] http://www.wsl.ch/info/mitarbeitende/moretti/download/DeCarceres_et_al_OIKOS2010 

.. _isaref6:

[6]


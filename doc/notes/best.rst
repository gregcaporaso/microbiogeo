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

:note: Need to follow up with the BIOENV and BVSTEP discovery. I do know R has an implementation of BIOENV. 


Existing Implementations
------------------------
There was found one existing implementation of the BEST test:

* Primer-e 

Primer-e is a proprietary implementation of an undocumented (under the name BEST)
statistic.


System Setup and Required Dependencies
--------------------------------------

Windows XP+ and the Primer-e software


Input Files
-----------

Sample data worksheet with variables tested to explain the community (active worksheet).
Resemblance matrix for the community to be explained (selected from drop down list). The actual samples considered are as selected in the active data worksheet.  Label matching between the worksheets is automatic but the resemblance matrix selection must contain at least these sample names.  They can be in any order. Data worksheet must have no missing values.


Output Files
------------

Not sure yet. Need to pursue the BIOENV and BVSTEP leads.

References
----------
.. _primere_bestref1:

[1] http://www.primer-e.com/


==============================================================
Mantel Test
==============================================================

Introduction
------------

This statistical method will test the correlation between two matrices. The data often represents the "distance" between objects. For example, one matrix might contain estimates of the genetic distances (i.e., the amount of difference between two different genomes) between all possible pairs of species in the study, obtained by the methods of molecular systematics; while the other might contain estimates of the geographical distance between the ranges of each species and every other species.

Existing Implementations
------------------------

Mantel has been implemented in python and has already been integrated into the QIIME package. The script is called compare_distance_matrices.py and was created by Dr. Caporaso.

System Setup and Required Dependencies
--------------------------------------

:note: The following instructions have been tested on Ubuntu 10.04.

Install the Quantative Insight Into Microbial Ecology (QIIME) software package. QIIME contains the compare_distance_matrices.py script.

Input Files
-----------

Required parameters

* Two matrices that will be compared. The script accepts tab delimited (*.txt) or comma delimited (*.csv) txt files. These text files should consist of matrices that have the same rank. 

* The file where the output will be written.

Optional parameters

* The number of iterations to perform

* Map of original sample ids to new sample ids


Execute the following command to run the script without additional parameters: ::

    compare_distance_matrices.py -i distance_matrix_1.txt,distance_matrix_2.txt -o mantel_output.txt
	
Execute the following command to run the script with additional parameters: ::

    compare_distance_matrices.py -i distance_matrix_1.txt,distance_matrix_2.txt -o mantel_output.txt -n 1000



Output Files
------------

The output will be in a tab seperated txt file with the following information:

* References to the files containing the distance matrices

* The number of entries

* The mantel p-value

Testing Results
---------------
This section will describe different tests that were run on the mantel script.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

88 Soils
^^^^^^^^^^
Test 1
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 2
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 3
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Face Site
^^^^^^^^^^

Test 1
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 2
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 3
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 2
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 3
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Keyboard
^^^^^^^^

Test 1
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 2
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

Test 3
~~~~~~

**Description:**

**Command:** ::

	Command here

**Results:**

References
----------

[1]
http://qiime.org/scripts/compare_distance_matrices.html
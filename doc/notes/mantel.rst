===========
Mantel Test
===========

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

These tests utilize the Python version of Mantel that was implemented in QIIME. The distance matrices used for these tests were derived from the 88 Soils dataset.

Test 1
~~~~~~

**Description:**

This test compares the PH distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm.txt,PH_dm.txt -o test1.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm.txt	PH_dm.txt	77	0.000

The p-value shows there is a very strong correlation between the two matrices.
	
Test 2
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_1.txt,PH_dm.txt -o test2.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_1.txt	PH_dm.txt	77	0.263

The p-value shows there is not a correlation between the two matrices.
	
Test 3
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_2.txt,PH_dm.txt -o test3.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_2.txt	PH_dm.txt	77	0.241

The p-value shows there is not a correlation between the two matrices.
	
Test 4
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_3.txt,PH_dm.txt -o test4.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_3.txt	PH_dm.txt	77	0.339

The p-value shows there is not a correlation between the two matrices.
	
Glen Canyon
^^^^^^^^^^^

These tests utilize the Python version of Mantel that was implemented in QIIME. The distance matrices used for these tests were derived from the Glen Canyon dataset.

Test 1
~~~~~~

**Description:**

This test compares the estimated years since submerged for plotting distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm.txt,estimated_years_since_submerged_for_plotting_dm.txt -o test1.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm.txt	estimated_years_since_submerged_for_plotting_dm.txt	94	0.000
	
The p-value shows there is a very strong correlation between the two matrices.

Test 2
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_1.txt,estimated_years_since_submerged_for_plotting_dm.txt -o test2.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_1.txt	estimated_years_since_submerged_for_plotting_dm.txt	94	0.442

The p-value shows there is not a correlation between the two matrices.
	
Test 3
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_2.txt,estimated_years_since_submerged_for_plotting_dm.txt -o test3.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_2.txt	estimated_years_since_submerged_for_plotting_dm.txt	94	0.762

The p-value shows there is not a correlation between the two matrices.
	
Test 4
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_shuffled_3.txt,estimated_years_since_submerged_for_plotting_dm.txt -o test4.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_shuffled_3.txt	estimated_years_since_submerged_for_plotting_dm.txt	94	0.539

The p-value shows there is not a correlation between the two matrices.
	
Keyboard
^^^^^^^^

These tests utilize the Python version of Mantel that was implemented in QIIME. The distance matrices used for these tests were derived from the Keyboard dataset.

Test 1
~~~~~~

**Description:**

This test compares the unweighted unifrac keyboard only 239 distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239.txt,unweighted_euclidean_dm.txt -o test1.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_keyboard_only_239.txt	unweighted_euclidean_dm.txt	74	0.197

The p-value shows there is not a correlation between the two matrices.
	
Test 2
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239_shuffled_1.txt,unweighted_euclidean_dm.txt -o test2.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_keyboard_only_239_shuffled_1.txt	unweighted_euclidean_dm.txt	74	0.363
	
The p-value shows there is not a correlation between the two matrices.

Test 3
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239_shuffled_2.txt,unweighted_euclidean_dm.txt -o test3.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_keyboard_only_239_shuffled_2.txt	unweighted_euclidean_dm.txt	74	0.426

The p-value shows there is not a correlation between the two matrices.
	
Test 4
~~~~~~

**Description:**

This test compares a shuffled unweighted unifrac distance matrix with the unweighted unifrac distance matrix using mantel.

**Command:** ::

	compare_distance_matrices.py -i unweighted_unifrac_dm_keyboard_only_239_shuffled_3.txt,unweighted_euclidean_dm.txt -o test4.txt -n 9999

**Results:** 

The following output file is created: ::

	# Number of entries refers to the number of rows (or cols) 
	# retained in each distance matrix after filtering the distance matrices 
	# to include only those samples that were in both distance matrices. 
	# p-value contains the correct number of significant digits.
	DM1	DM2	Number of entries	Mantel p-value
	unweighted_unifrac_dm_keyboard_only_239_shuffled_3.txt	unweighted_euclidean_dm.txt	74	0.683

	
The p-value shows there is not a correlation between the two matrices.
	
References
----------

[1]
http://qiime.org/scripts/compare_distance_matrices.html

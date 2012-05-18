===========================================================
Permutational Multivariate Analysis of Variance (PERMANOVA)
===========================================================


Introduction
------------

The permutational multivariate analysis of variance (PERMANOVA) method is a computer program for testing 
the simultaneous response of one or more variables to one or more factors in an ANOVA experimental design
on the basis of any distance metric. It returns a F value and a P value.

The first thing it does is calculate the distances between each pair of sampled units to obtain a distance matrix.
It then calculates the test-statistics from this according to the relevant experimental design.
See Anderson (2001a) and McArdle and Anderson (2001).


Existing Implementations
------------------------



System Setup and Required Dependencies
--------------------------------------

:note: The following instructions have been tested on 64-bit Fedora (essentially Redhat) using Python 2.7.2. However, they `should` work across different linux distros and on Macs. The instructions assume you use bash as your shell.

First, your system must have a version of QIIME installed (I used the latest
version of QIIME in SVN). The code also depends on NumPy, though this should
already be installed if you have QIIME installed. Next, you must add the area
where Andrew's code resides to your PYTHONPATH environment variable, changing
the path to point to the location of the microbiogeo checkout on your machine. I
also added the scripts area to my PATH for convenience: ::

    export PYTHONPATH=/home/aragorn/workspace/Bio:$PYTHONPATH
    export PATH=/home/aragorn/workspace/Bio:$PATH

If you don't want to have to perform this step each time you open a new
terminal, run the following command to add the path to your .bashrc: ::

    echo "export PYTHONPATH=/home/aragorn/workspace/Bio:$PYTHONPATH" >> ~/.bashrc
    echo "export PATH=/home/aragorn/workspace/Bio:$PATH" >> ~/.bashrc
    source ~/.bashrc

Next install lapack-devel and it's dependencies: ::

	yum install lapack-devel 

Libgc2.so needs to be in /usr/lib: ::

	yum install compat-gcc-34-g77
	ln -s /usr/lib64/libg2c.so.0 libg2c.so

Install pysparse 1.1-1 as found at https://launchpad.net/ubuntu/+source/pysparse/1.1-1: ::

	wget https://launchpad.net/ubuntu/+archive/primary/+files/pysparse_1.1.orig.tar.gz
	tar -xzvf pysparse_1.1.orig.tar.gz 
	cd pysparse-1.1/
	python setup.py install
	

Run the following command to test if you can run the Permanova script: ::

    ./permanova.py -h

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The Permanova script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix from the QIIME overview tutorial. You can get the distance
matrix :download:`here <../../doc/downloads/overview_unweighted_unifrac_dm.txt>` and
the mapping file :download:`here <../../doc/downloads/Fasting_Map.txt>`.

Next, run the following command to execute the Permanova script: ::

	./permanova.py -i unweighted_unifrac_dm.txt -m map.txt -c BODY_SITE -o permanova_results.txt -p 999

The -c option specifies which column in the mapping file will be used to group
the samples. The `Treatment` column has two values: 'Control' and 'Fast'. Thus,
Permanova will be used to calculate the dissimilarity between the control and fast
groups. The -o option specifies the file that we want the results written to.

Output Files
------------
The command in the previous section creates a single output file named
:file:`permanova_results.txt`. The resulting file should look like this: ::

	Input_filepath						permanova_F_value	p_value
	unweighted_unifrac_dm.txt           2.29665065171		0.009
	
The first field lists the distance matrix file that was used as input. The
second field lists the F statistic that was computed, that this is the
primary output of Permanova. The final field lists the p-value, which is .009 
signifying reasonable accuracy.

Testing Results
---------------
This section will describe different tests that were run on the PERMANOVA script.
These tests will use empirical data from one of the several datasets that the team 
has access to. These data files will not be included for download due to their
(usually) large size. Unless otherwise noted, the data files that were used can be 
found under the datasets directory.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `BODY_SITE` category as a positive control. We expect there
to be significant differences between body sites due to previous analysis done
on the Whole Body dataset.

**Command:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

	Input_filepath						permanova_F_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	13.2670596158		0.001
	
The F value of 13.2670596158 indicates that the body sites are significantly different due to its 
relatively “large” ratio. This is a result that we would expect. The p-value of .001 indicates that the 
result is significant.

Test 2
~~~~~~
**Description:**

This test uses the `SEX` category as a negative control.
We don't expect there to be significant differences based on sex due to previous
analysis done on the Whole Body dataset.

**Command:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c SEX -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

	Input_filepath						permanova_F_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	21.0188242485		0.001
	
The F value of 21.0188242485 indicates that the sexes are significantly different due to its 
relatively “large” ratio. This is a result that we would not expect because it is supposed to be a negative control. The p-value of .001 indicates that the 
result is significant.

Test 3
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the BODY_SITE category to perform three negative control 
tests. Since the labels of the distance matrices are shuffled, we don’t expect to see significant differences any more on 
this category.

**Commands:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999
	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999
	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999

**Results:**

The following output files were created: ::

	Input_filepath							permanova_F_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt	1.98060081904		0.031

::
	
	Input_filepath							permanova_F_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt	1.81015551855		0.623

::
		
	Input_filepath							permanova_F_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt	1.73759470202		0.929
	
The F values of 1.98060081904, 1.81015551855,and 1.73759470202 indicates that the body sites are not significantly different due to its 
relatively “small” ratio. This is a result that we would expect because the matricies are pre-shuffled. The p-values of 0.031, 0.623, and 0.929 indicates that the 
results are insignificant.

Keyboard
^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category as a positive control. We expect
there to be significant differences based on host subjects due to previous analysis
done on the keyboard study dataset.

**Command:** ::

    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

        Input_filepath                                          permanova_F_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm.txt       5.17880475397           0.001
	
The F value of 5.17880475397 indicates that the host id's are significantly different due to its 
relatively “large” ratio. This is a result that we would expect. The p-value of 0.001 indicates that the 
result is significant.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `HOST_SUBJECT_ID`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see significant differences any more on
this category.

**Command:** ::

    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_1.txt -p 999
    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_2.txt -p 999
    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_3.txt -p 999

**Results:**

The following output files are created: ::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt    1.04303546137           0.31

::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt    1.03699740907           0.317

::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt    0.959082333436          0.648
	
The F values of 1.04303546137, 1.03699740907, and 0.959082333436 indicates that the shuffled host id's are not significantly different due to its 
relatively “small” ratio. This is a result that we would expect. The p-value of 0.31, 0.317, 0.648 and indicates that the 
result are insignificant.

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category as a positive control. We expect
there to be significant differences on this category due to previous analysis
done on the Glen Canyon dataset.

**Command:** ::

    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

        Input_filepath                                          permanova_F_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm.txt    29.2130439798           0.001
	
The F value of 29.2130439798 indicates that there is significant differences between groups of samples based on whether they are currently wet, we can tell due to its 
relatively “large” ratio. This is a result that we would expect from previous experments. The p-value of 0.001 indicates that the 
result is significant.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `CurrentlyWet`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see significant differences any more on
this category.

**Command:** ::

    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999
    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999
    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999

**Results:**

The following output files are created: ::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt 2.03268471405           0.332

::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt 1.91859583429           0.448

::

        Input_filepath                                                  permanova_F_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt 1.68346774545           0.902
	
The F value of 2.03268471405, 1.91859583429, and 1.68346774545 indicates that the results when CurrentlyWet is shuffled are not significantly different due to its 
relatively “small” ratio. This is a result that we would expect. The p-values of 0.332, 0.448, and 0.902 indicates that the 
results are insignificant.

References
----------
.. _permanovaref1:

[1] www.stat.auckland.ac.nz/~mja/prog/PERMANOVA_UserNotes.pdf


========================================================================================
Permutational Multivariate Analysis of Variance (PERMANOVA) Statistical Method Reference
========================================================================================


Introduction
------------

The permutational multivariate analysis of variance (PERMANOVA) method is a computer program for testing 
the simultaneous response of one or more variables to one or more factors in an ANOVA experimental design
on the basis of any distance metric. It returns a R value and a P value.

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

    ./permanova.py -i ../../doc/downloads/overview_unweighted_unifrac_dm.txt -m ../../doc/downloads/Fasting_Map.txt -c Treatment -o permanova_results.txt

The -c option specifies which column in the mapping file will be used to group
the samples. The `Treatment` column has two values: 'Control' and 'Fast'. Thus,
Permanova will be used to calculate the dissimilarity between the control and fast
groups. The -o option specifies the file that we want the results written to.

Output Files
------------
The command in the previous section creates a single output file named
:file:`permanova_results.txt`. The resulting file should look like this: ::

	Input_filepath				permanova_R_value	p_value
	overview_unweighted_unifrac_dm.txt	2.29665065171		NA

The first field lists the distance matrix file that was used as input. The
second field lists the R statistic that was computed (remember that this is the
primary output of Permanova). The final field lists the p-value, which is NA
because we did not specify the optional -p parameter (by default, the number of
p-trials is 0).

Investigate why not like Jai's: ::

	The value of the R statistic can fall between -1 and +1, with a positive value
	close to 1 indicating that the groups are highly dissimilar. Thus, in this
	example, the control and fast groups are dissimilar. 


Testing Results
---------------
This section will describe different tests that were run on the PERMANOVA script.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `BODY_SITE` category as a positive control.
We expect there to be significant clustering due to previous analysis done on
the Whole Body dataset.

**Command:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	13.2670596158		0.001
	

Test 2
~~~~~~
**Description:**

This test uses the `SEX` category as a negative control.
We expect there to be less clustering due to previous analysis done on
the Whole Body dataset.

**Command:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c SEX -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	21.0188242485		0.001


Test 3
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the BODY_SITE category to perform three negative control 
tests. Since the labels of the distance matrices are shuffled, we don’t expect to see clustering any more on 
this category.

**Commands:** ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999
	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999
	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/whole_body/map.txt -c BODY_SITE -o permanova_results.txt -p 999

**Results:**

The following output files were created: ::

	Input_filepath							permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt	1.98060081904		0.031

::
	
	Input_filepath							permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt	1.81015551855		0.623

::
		
	Input_filepath							permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt	1.73759470202		0.929

Keyboard
^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category as a positive control. We expect
there to be significant clustering on host subjects due to previous analysis
done on the keyboard study dataset.

**Command:** ::

    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

        Input_filepath                                          permanova_R_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm.txt       5.17880475397           0.001

The R value of 0.794026410205 indicates that samples taken from different hosts
are significantly different (i.e. there is clustering) due to its "large"
positive value. This is a result that we would expect. The p-value of 0.001
indicates that the result is significant.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `HOST_SUBJECT_ID`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_1.txt -p 999
    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_2.txt -p 999
    ./permanova.py -i ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova_results_3.txt -p 999

**Results:**

The following output files are created: ::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt    1.04303546137           0.31

::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt    1.03699740907           0.317

::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt    0.959082333436          0.648


Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category as a positive control. We expect
there to be significant clustering on this category due to previous analysis
done on the Glen Canyon dataset.

**Command:** ::

    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999

**Results:**

The following output file is created: ::

        Input_filepath                                          permanova_R_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm.txt    29.2130439798           0.001

The R value of 0.9984007035 indicates that samples taken from wet and dry
environments are significantly different (i.e. there is clustering) due to the
really "large" positive value that is close to 1. This is a result that we would
expect, as there is also clear clustering in the 3D PCoA plots.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `CurrentlyWet`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999
    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999
    ./permanova.py -i ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m ../../datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o permanova_results.txt -p 999

**Results:**

The following output files are created: ::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt 2.03268471405           0.332

::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt 1.91859583429           0.448

::

        Input_filepath                                                  permanova_R_value       p_value
        ../../datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt 1.68346774545           0.902

References
----------

Jai's permanova.rst

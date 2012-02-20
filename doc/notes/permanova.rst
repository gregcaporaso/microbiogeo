==============================================================
Permutational Multivariate Analysis of Variance (PERMANOVA) Statistical Method Reference
==============================================================

Introduction
------------


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

	Input_filepath	permanova_R_value	p_value
	overview_unweighted_unifrac_dm.txt	2.29665065171	NA

The first field lists the distance matrix file that was used as input. The
second field lists the R statistic that was computed (remember that this is the
primary output of Permanova). The final field lists the p-value, which is NA
because we did not specify the optional -p parameter (by default, the number of
p-trials is 0).

Investigate why not like Jai's: ::

	The value of the R statistic can fall between -1 and +1, with a positive value
	close to 1 indicating that the groups are highly dissimilar. Thus, in this
	example, the control and fast groups are dissimilar. 


Whole Body
----------

Evaluated by ORIGINAL_SAMPLE_SITE: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c ORIGINAL_SAMPLE_SITE -o permanova_results.txt

Returns: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	9.05965249449		NA
	
	
Evaluated by TARGET_SUBFRAGMENT: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c TARGET_SUBFRAGMENT -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars
	
Evaluated by COMMON_SAMPLE_SITE: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c COMMON_SAMPLE_SITE -o permanova_results.txt

Returns: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	9.05965249449		NA
	
Evaluated by COMMON_NAME: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c COMMON_NAME -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars
	
Evaluated by BIOLOGICAL_SPECIMEN: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c BIOLOGICAL_SPECIMEN -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars
	
Evaluated by SAMP_SIZE: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c SAMP_SIZE -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars
	
Evaluated by BODY_HABITAT: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c BODY_HABITAT -o permanova_results.txt

Returns: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	30.1660567245		NA
	
Evaluated by SEQUENCING_METH: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c SEQUENCING_METH -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars
	
Evaluated by SAMP_COLLECT_DEVICE: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c SAMP_COLLECT_DEVICE -o permanova_results.txt

Returns: ::

	FloatingPointError: divide by zero encountered in double_scalars

Evaluated by ENV_MATTER: ::

	./permanova.py -i ../../datasets/whole_body/unweighted_unifrac_dm.txt -m ../../datasets/whole_body/map.txt -c ENV_MATTER -o permanova_results.txt

Returns: ::

	Input_filepath						permanova_R_value	p_value
	../../datasets/whole_body/unweighted_unifrac_dm.txt	28.067284235		NA
	

References
----------

Jai's anosim.rst

======================================================
Spacial Autocorrelation (Morans I)
======================================================
Introduction
------------
Morans I is also known as Morans autocorrelation coefficient. It is an extension of Pearsons product moment correlation coefficient (Pearsons r), which measures the strength of linear dependence between two variables. Morans I returns a value, sometimes referred to as the I value, that shows the relationship between samples in different areas.

The Moran's I equation is:

.. image:: ../images/morans_i/MoransIEquation.PNG
  :align: center

... where wij is the weight between observation i and j, and S0 is the sum of all wij 's. Which can be seen here:

.. image:: ../images/morans_i/MoransISubVar.PNG
  :align: center

This equation returns a value between -1 and 1. The value indicates what type of correlation to expect for the data analyzed. In the case of -1 it represents perfect dispersion, where there is an equal amount of variation throughout. A positive 1 is called perfect correlation, where two points always occur together. Lastly, when you get a value closer to 0 that indicates that there is random patterning being detected.

What the null hypothesis tests for when using Morans I, is if there's no spacial correlation among communities.

Selected Implementation
-----------------------

The implementation being used is from the R software package.

The documentation for the Morans I method in R can be found :ref:`here <moranref1>`

As well as the Gittleman and Kot paper referenced from the documentation can be found :ref:`here <moranref3>`

Alternative Implementations
---------------------------

There is a matlab version of Morans I floating around. It is written by a Felix Hebeler. The last update was in October 2008 :ref:`[4] <moranref4>`.

Package sp has several functions, including moran.test, that are more specifically targeted to the analysis of spatial data. Package spatial has the function correlogram that computes and plots spatial correlograms.:ref:[`5 <moranref5>`]


Input
-----
In order to use this method you will need the distance matrices(weighted or unweighted) and the mapping file that the information was created from. Lastly you need to select a non numerical type of data to test against. One example is date of birth.

From the command line type: 

``R --slave --args -d unifrac.txt -m Fasting_Map.txt -c DOB -o morans_i < morans_i.r``

REQUIRED options:
The following options must be provided under all circumstances.

``--slave``
    Make R run as quietly as possible. This option is intended to support programs which use R to compute results for them. It implies --quiet and --no-save. 

``--args``
    This flag does nothing except cause the rest of the command line to be skipped: this can be useful to retrieve values from it with commandArgs(TRUE).

``-i OR --input_path = INPUT_PATH``
	path to the input distance matrix file(s) (i.e., the output from beta_diversity.py).

``-o OR --output_path = OUTPUT_PATH``
	output path to the name of a single file

``-m OR --map_path = MAP_PATH``
	path to the location of the mapping file

``-c CATEGORY, --category=CATEGORY``
	String which coresponds to the column name containing grouping info

Output
------
The output of Morans I is a file that is placed in a directory specified by the -o argument. The file will be a text file with 4 values: observed, expected, sd, and p.value.

The observed value is Morans I index of x. This is computed based on the values passed in to be compared with the weights.

The expected value is the value of I under the null hypothesis.

The sd is the standard deviation of I under the null hypothesis.

P Value is the p-value of the test of the null hypothesis against the alternative hypothesis specified in alternative

Each of these values, except for the p-value, should be between -1 and 1. 

Testing Results
---------------
Testing needs to be performed further to understand the results and what I should be expecting as output. In the mean time, using the `QIIME Overview <http://qiime.org/tutorials/tutorial.html>`_ data I ran Morans I and received the following back:

===========  ===========  ===========  ===========
observed     expected     sd           p.value
===========  ===========  ===========  ===========
-0.06005486  -0.125       0.01590547   4.442088e-05
===========  ===========  ===========  ===========

I'm not especially sure how accurate this is, but it's useful to see what was received.

Testing With Whole Body Datasets
--------------------------------

I can't get my testing information to be consistent, as such I need to continue working on it.


System Setup and Required Dependencies
--------------------------------------
Step 1:
The first step is to install R. The following command downloads and installs R:

    sudo apt-get install r-base

Step 2:
Identify the qiime location for where it is installed. In the case of the AWS, using AMI:QIIME 1.4.0 EBS East (ami-458d5b2c). 

	QIIME location is: /software/qiime-1.4.0-release

Step 3:
You need to define an environment variable to tell the script where to look for the r utility functions in qiime. Run the following command, changing the path to point to the location of your qiime install:

    export qiime_dir=/home/<username>/qiime/trunk

If you dont want to have to perform this step each time you open a new terminal, run the following command to add it to your .bashrc:

    echo "export qiime_dir=/home/<username>/qiime/trunk" >> ~/.bashrc
    source ~/.bashrc

OR

Go into /etc/, and open the file /etc/environment. In this file youll want to put the line:

	QIIME_DIR="/software/qiime-1.4.0-release" 

The full information is:

	Directory: /etc/
	File: environment
	Full file path: /etc/environment
	String to add at bottom: QIIME_DIR="/software/qiime-1.4.0-release" 

Make sure to include the quotes. Once you do that you need to save and  restart. 

After all of this you can now type "echo $QIIME_DIR" in the terminal and it should print out the set path that was used above..

Step 4:
Youll need to install some R packages. If you can use the R console from the command line simply type R to get to it.

To get the packages type:
	install.packages(optparse)
	install.packages(ape)

If youre concerned about updating packages type "update.packages()" in the R console, excluding the quotes.

References
----------
.. _moranref1:

[1]R Documentation for Morans I

http://svitsrv25.epfl.ch/R-doc/library/ape/html/MoranI.html

.. _moranref2:

[2]How to Work with Morans I in R

http://www.ats.ucla.edu/stat/r/faq/morans_i.htm

.. _moranref3:

[3]Gittleman and Kot paper

http://www.jstor.org/pss/2992183

.. _moranref4:

[4]Hebeler Morans I version

http://www.mathworks.com/matlabcentral/fileexchange/13663-morans-i/content/morans_I.m

.. _moranref5:

[5]Morans I Paper by Emmanuel Paradis

http://cran.r-project.org/web/packages/ape/vignettes/MoranI.pdf


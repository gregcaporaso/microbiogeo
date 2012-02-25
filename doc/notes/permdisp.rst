============================================================================================================================
Permutational Analysis of Multivariate Dispersions (PERMDISP) Statistical Method Reference
============================================================================================================================

Introduction
------------

This program compares within-group multivariate dispersions (spread or variability) among groups on the basis of any distance or dissimilarity measure using permutations.

This statistical method has two steps 

* calculation of the distances from observations to their centroids 

* comparison of the average of these distances among groups, using ANOVA

A P-value is then obtained using permutation of the observations. The approach is a multivariate analogue to Levene test (Levene 1960). The analysis of these distances to centroids can be done for any two-factor design, just as in PERMANOVA. It is wellknown that the test for differences in location among groups in multivariate space (such as PERMANOVA) is sensitive to differences in dispersion among the groups. Thus, rejection of the null hypothesis for PERMANOVA suggests that groups may differ because of their location, their relative dispersion, or both. 
The program is designed to be used as a companion to PERMANOVA, to assist in unraveling the possible reasons for rejection of the null hypothesis. It is also useful in its own right, however, to investigate differences in dispersions alone when hypotheses of this nature arise.

Existing Implementations
------------------------

The only implementation of PERMDISP was a FORTRAN program that is offered by:

<http://www.stat.auckland.ac.nz/~mja/Programs.htm>

System Setup and Required Dependencies
--------------------------------------

:note: The following instructions have been tested on 32-bit Windows Vista, 64-bit Windows 7, and Ubuntu 10.04. There is also a distribution for Macs, but it has not been tested.

Download the PERMDISP that relates to your local environment (PERMDIPS.exe for Windows or Linux / PERMDISP.mac.sit for Mac). 

<http://www.stat.auckland.ac.nz/~mja/Programs.htm>

Once downloaded run the executible to start the program.

For Linux users:

Install wine or other software to run exe files in a Linux enviroment.

To do this I ran the following command in Ubuntu: ::

	sudo apt-get install wine1.2
	
Once software is installed to open exe files run the exe file.

Input Files
-----------

Once the exe file is ran a prompt will come up requesting information.

Input a file containing raw data or a file containing a symmetric matrix of distances or 
dissimilarities. In each case, the file must be saved in tab delimited (ASCII text, *.txt) or comma delimited 
(*.csv) format.

The program will prompt the user to name the file where the results will be written. Though it is not necessary it is encouraged to name the output file with the extension .txt.

There are many different specifications that the user can control. To cover all of these would take a considerable amount of time. Because of this one possible path is displayed below.

Nature of the data in the input file: 
1) raw data (n x p) 
2) distance matrix (n x n)
2

ANOVA Experimental design: 
1) One-way 
2) Two-way nested 
3) Two-way crossed (i.e. factorial or orthogonal) 
3

Experimental design of two-way crossed analysis: 
1) Fixed effects - both factors are fixed 
2) Random effects - both factors are random 
3) Mixed model - factor 1 is fixed, 2 is random 
4) Mixed model - factor 1 is random, 2 is fixed 
1 

What is the name of factor 1? 
distance

Type the number of levels for factor 1 
1

What is the name of factor 2? 
difference 

Type the number of levels for factor 2 
1
 
What is the number of replicates? 
4
 
How many permutations do you want for the tests? (i.e. 99, 499, 999, 4999, etc.) 
99
 
Type an integer to be used as the seed 
for the random permutations 
9


Output Files
------------

The output file will contain the "Experimental Design" that contains the information that the user input when the program was ran. It will also contain the result of the tests for heterogeneity in the average dissimilarities of points from the central location of their group, which will be displayed in an ASCII table.

References
----------

[1]
http://www.stat.auckland.ac.nz/~mja/Programs.htm

[2]
http://www.stat.auckland.ac.nz/~mja/prog/PERMDISP_UserNotes.pdf
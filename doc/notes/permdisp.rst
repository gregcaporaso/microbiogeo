==============================================================
Permutational Analysis of Multivariate Dispersions (PERMDISP) Statistical Method Reference
==============================================================

Introduction
------------

This program compares within-group multivariate dispersions (spread or variability) among groups on the basis of any distance or dissimilarity measure using permutations.

This statistical method has two steps 

* calculation of the distances from observations to their centroids 

* comparison of the average of these distances among groups, using ANOVA.

A P-value is then obtained using permutation of the observations. The approach is a multivariate analogue to Levene’s test (Levene 1960). The analysis of these distances to centroids can be done for any two-factor design, just as in PERMANOVA. It is wellknown that the test for differences in location among groups in multivariate space (such as PERMANOVA) is sensitive to differences in dispersion among the groups. Thus, rejection of the null hypothesis for PERMANOVA suggests that groups may differ because of their location, their relative dispersion, or both. 
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

The program accepts tab delimited (*.txt) or comma delimited (*.csv) formats.

Output Files
------------

Though it is not necessary it is encouraged to name the output file with the extension *.txt.

References
----------

[1]
http://www.stat.auckland.ac.nz/~mja/Programs.htm

[2]
http://www.stat.auckland.ac.nz/~mja/prog/PERMDISP_UserNotes.pdf
=============================================================
Permutational Analysis of Multivariate Dispersions (PERMDISP)
=============================================================

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

Selected Implementation
^^^^^^^^^^^^^^^^^^^^^^^

BETADISP (PERMDISP): an R program that is in the vegan package:

<http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/betadisper.html>

Alternative Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PERMDISP: FORTRAN program that is offered by:

<http://www.stat.auckland.ac.nz/~mja/Programs.htm>

System Setup and Required Dependencies
--------------------------------------

BETADISPER (R)
^^^^^^^^^^^^^^

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
	install.packages(vegan)

If youre concerned about updating packages type "update.packages()" in the R console, excluding the quotes.

PERMDISP (FORTRAN)
^^^^^^^^^^^^^^^^^^

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

BETADISPER (R)
^^^^^^^^^^^^^^

The method call from the R software package that is used with this script is:

``betadisper(formula, data, permutations = 999, method = "bray", strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", ...)``

Formula - This represents what the data is going to be fit to
Data - This is the data being used for this method
Permutations - This is the number of replications used for hypothesis tests.
Method - The specified manner in which pair wise distances are calculated
Strata - This groups the permutations based on the specified strata

From the command line type: ::

  R --slave --args -d distanceMatrix.txt -m Fasting_Map.txt -c Treatment -o betadisper < betadisper.r

REQUIRED script options:
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

PERMDISP (FORTRAN)
^^^^^^^^^^^^^^^^^^

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

BETADISPER (R)
^^^^^^^^^^^^^^

The output for permdisp is in a directory specified by the -o parameter. The results should be labeled "betadisper_results.txt".

PERMDISP (FORTRAN)
^^^^^^^^^^^^^^^^^^

The output file will contain the "Experimental Design" that contains the information that the user input when the program was ran. It will also contain the result of the tests for heterogeneity in the average dissimilarities of points from the central location of their group, which will be displayed in an ASCII table.

Testing Results
---------------
This section will describe different tests that were run on the BETADISPER implementation.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

All the tests below were done using **BETADISPER (R)** because PERMDISP (FORTRAN) did not allow the use of a mapping file.

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~

**Description:**

This test uses the original unweighted unifrac distance matrix and the CurrentlyWet category as a positive control.

**Command:** ::

	R --slave --args -d Glen\ Canyon/unweighted_unifrac_dm.txt -m Glen\ Canyon/map_25Jan2012.txt -c CurrentlyWet -o betadisper_positive < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			  Df    Sum Sq   Mean Sq F value   Pr(>F)    
	Groups     1 0.0060076 0.0060076  26.742 1.35e-06 ***
	Residuals 92 0.0206680 0.0002247                     
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			  Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
	Groups     1 0.0060076 0.0060076 26.742    999  0.001 ***
	Residuals 92 0.0206680 0.0002247                         
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
				No   Yes
	No             0.001
	Yes 1.3501e-06      
	
The p-value indicates that the results are significant.
	
Test 2
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	R --slave --args -d Glen\ Canyon/unweighted_unifrac_dm_shuffled_1.txt -m Glen\ Canyon/map_25Jan2012.txt -c CurrentlyWet -o betadisper_negative_1 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			  Df   Sum Sq    Mean Sq F value Pr(>F)
	Groups     1 0.000878 0.00087764  0.3079 0.5803
	Residuals 92 0.262210 0.00285011               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			  Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.000878 0.00087764 0.3079    999  0.589
	Residuals 92 0.262210 0.00285011                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			No   Yes
	No         0.586
	Yes 0.5803      

The p-value indicates that the results are insignificant.
	
Test 3
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	R --slave --args -d Glen\ Canyon/unweighted_unifrac_dm_shuffled_2.txt -m Glen\ Canyon/map_25Jan2012.txt -c CurrentlyWet -o betadisper_negative_2 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			  Df   Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.002333 0.0023333  0.8033 0.3725
	Residuals 92 0.267228 0.0029046               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			  Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.002333 0.0023333 0.8033    999  0.381
	Residuals 92 0.267228 0.0029046                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			 No   Yes
	No          0.387
	Yes 0.37245      

The p-value indicates that the results insignificant.
	
Test 4
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	R --slave --args -d Glen\ Canyon/unweighted_unifrac_dm_shuffled_3.txt -m Glen\ Canyon/map_25Jan2012.txt -c CurrentlyWet -o betadisper_negative_3 < betadisper.r
	
**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			  Df   Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.001018 0.0010178  0.3552 0.5526
	Residuals 92 0.263611 0.0028653               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			  Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.001018 0.0010178 0.3552    999  0.526
	Residuals 92 0.263611 0.0028653                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			 No   Yes
	No          0.526
	Yes 0.55264      
	
The p-value indicates that the results are insignificant.
	
Keyboard
^^^^^^^^

Test 1
~~~~~~

**Description:**

This test uses the original unweighted unifrac distance matrix and the HOST_SUBJECT_ID category as a positive control.

**Command:** ::

	R --slave --args -d Keyboard/unweighted_unifrac_dm.txt -m Keyboard/map.txt -c HOST_SUBJECT_ID -o betadisper_positive < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df    Sum Sq    Mean Sq F value    Pr(>F)    
	Groups     10 0.0089615 0.00089615  4.4126 3.586e-05 ***
	Residuals 104 0.0211214 0.00020309                      
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)   
	Groups     10 0.0089615 0.00089615 4.4126    999  0.002 **
	Residuals 104 0.0211214 0.00020309                        
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			 F1 L1 L3       M1       M2       M3       M9 R1 U1 U2 U3
	F1                0.736000 0.403000 0.241000 0.159000            
	L1                                                               
	L3                                                               
	M1 0.718987                0.243000 0.098000 0.127000            
	M2 0.388167       0.298617          0.691000 0.065000            
	M3 0.225243       0.110657 0.693715          0.016000            
	M9 0.158122       0.147042 0.074368 0.018900                     
	R1                                                               
	U1                                                               
	U2                                                               
	U3                                                  

The p-value indicates that the results are significant.	

Test 2
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the HOST_SUBJECT_ID category to perform a negative control test.

**Command:** ::

	R --slave --args -d Keyboard/unweighted_unifrac_dm_shuffled_1.txt -m Keyboard/map.txt -c HOST_SUBJECT_ID -o betadisper_negative_1 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value    Pr(>F)    
	Groups     10 0.024535 0.00245353  4.4774 2.968e-05 ***
	Residuals 104 0.056990 0.00054798                      
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)    
	Groups     10 0.024535 0.00245353 4.4774    999  0.001 ***
	Residuals 104 0.056990 0.00054798                         
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			  F1 L1 L3        M1        M2        M3        M9 R1 U1 U2 U3
	F1                 0.3710000 0.2000000 0.0680000 0.3010000            
	L1                                                                    
	L3                                                                    
	M1 0.3554280                 0.0310000 0.0010000 0.0550000            
	M2 0.1838147       0.0306839           0.8610000 0.7230000            
	M3 0.0588469       0.0014799 0.8623054           0.7980000            
	M9 0.2662109       0.0551588 0.7268720 0.8137601                      
	R1                                                                    
	U1                                                                    
	U2                                                                    
	U3                                                                    

The p-value indicates that the results are significant.

Test 3
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the HOST_SUBJECT_ID category to perform a negative control test.

**Command:** ::

	R --slave --args -d Keyboard/unweighted_unifrac_dm_shuffled_2.txt -m Keyboard/map.txt -c HOST_SUBJECT_ID -o betadisper_negative_2 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value    Pr(>F)    
	Groups     10 0.024182 0.00241823  3.5249 0.0004881 ***
	Residuals 104 0.071348 0.00068604                      
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)    
	Groups     10 0.024182 0.00241823 3.5249    999  0.001 ***
	Residuals 104 0.071348 0.00068604                         
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			 F1 L1 L3       M1       M2       M3       M9 R1 U1 U2 U3
	F1                0.353000 0.969000 0.776000 0.951000            
	L1                                                               
	L3                                                               
	M1 0.346654                0.025000 0.006000 0.060000            
	M2 0.960319       0.016261          0.485000 0.726000            
	M3 0.757185       0.004107 0.488713          0.351000            
	M9 0.943114       0.050123 0.737238 0.351259                     
	R1                                                               
	U1                                                               
	U2                                                               
	U3                                                               

The p-value indicates that the results are significant.

Test 4
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the HOST_SUBJECT_ID category to perform a negative control test.

**Command:** ::

	R --slave --args -d Keyboard/unweighted_unifrac_dm_shuffled_3.txt -m Keyboard/map.txt -c HOST_SUBJECT_ID -o betadisper_negative_3 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value    Pr(>F)    
	Groups     10 0.024199 0.00241989  3.7129 0.0002801 ***
	Residuals 104 0.067783 0.00065176                      
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)    
	Groups     10 0.024199 0.00241989 3.7129    999  0.001 ***
	Residuals 104 0.067783 0.00065176                         
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
			F1 L1 L3      M1      M2      M3      M9 R1 U1 U2 U3
	F1               0.96500 0.36700 0.35400 0.13000            
	L1                                                          
	L3                                                          
	M1 0.96735               0.47100 0.47400 0.20500            
	M2 0.36282       0.47355         0.76700 0.27500            
	M3 0.33691       0.44413 0.78665         0.16500            
	M9 0.12918       0.21136 0.28776 0.15521                    
	R1                                                          
	U1                                                          
	U2                                                          
	U3                               

The p-value indicates that the results are significant.	

Whole Body
^^^^^^^^^^

Test 1
~~~~~~

**Description:**

This test uses the original unweighted unifrac distance matrix and the BODY_SITE category as a positive control.

**Command:** ::

	R --slave --args -d Whole\ Body/unweighted_unifrac_dm.txt -m Whole\ Body/map.txt -c BODY_SITE -o betadisper_positive < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq   Mean Sq F value    Pr(>F)    
	Groups     19 0.092975 0.0048934  22.892 < 2.2e-16 ***
	Residuals 565 0.120776 0.0002138                      
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
	Groups     19 0.092975 0.0048934 22.892    999  0.001 ***
	Residuals 565 0.120776 0.0002138                         
	---
	Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
									 UBERON:ear canal UBERON:feces
	UBERON:ear canal                                    1.0000e-03
	UBERON:feces                           7.0921e-11             
	UBERON:glans penis                     5.2830e-01   1.0172e-12
	UBERON:hair                            9.8369e-01   9.8136e-17
	UBERON:labia minora                    7.0630e-01   1.2191e-10
	UBERON:mouth                           3.4896e-04   2.9259e-02
	UBERON:nose                            9.6989e-01   2.3641e-15
	UBERON:nostril                         5.4606e-01   4.5769e-16
	UBERON:nostrils                        6.4419e-01   7.9374e-16
	UBERON:skin of arm                     8.9785e-01   9.4252e-22
	UBERON:skin of finger                  4.3878e-03   3.9370e-12
	UBERON:skin of forearm                 8.3154e-04   4.2747e-08
	UBERON:tongue                          1.9110e-08   8.1147e-01
	UBERON:urine                           3.3623e-01   3.6922e-14
	UBERON:zone of skin of abdomen         4.7309e-01   3.1765e-22
	UBERON:zone of skin of foot            4.9330e-01   1.3516e-13
	UBERON:zone of skin of hand            9.6965e-05   9.1352e-13
	UBERON:zone of skin of head            7.9705e-01   1.9681e-20
	UBERON:zone of skin of knee            1.9000e-06   1.4137e-12
	UBERON:zone of skin of outer ear       7.2754e-02   6.8148e-28
									 UBERON:glans penis UBERON:hair
	UBERON:ear canal                         5.2400e-01  9.8600e-01
	UBERON:feces                             1.0000e-03  1.0000e-03
	UBERON:glans penis                                   3.4600e-01
	UBERON:hair                              3.6136e-01            
	UBERON:labia minora                      8.2297e-01  5.9163e-01
	UBERON:mouth                             8.0518e-05  1.6235e-06
	UBERON:nose                              4.2323e-01  9.4077e-01
	UBERON:nostril                           6.9718e-01  4.4414e-01
	UBERON:nostrils                          6.4347e-01  5.3491e-01
	UBERON:skin of arm                       3.1778e-01  8.4270e-01
	UBERON:skin of finger                    4.4064e-02  2.3261e-04
	UBERON:skin of forearm                   5.7989e-03  1.6493e-05
	UBERON:tongue                            1.0720e-09  5.7338e-13
	UBERON:urine                             1.7974e-01  3.0509e-01
	UBERON:zone of skin of abdomen           8.1279e-01  2.7965e-01
	UBERON:zone of skin of foot              7.8801e-01  4.3252e-01
	UBERON:zone of skin of hand              1.5732e-02  2.8212e-06
	UBERON:zone of skin of head              2.2295e-01  7.8098e-01
	UBERON:zone of skin of knee              6.7749e-06  5.0286e-10
	UBERON:zone of skin of outer ear         1.5410e-02  4.6130e-02
									 UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                          6.9500e-01   1.0000e-03  9.6700e-01
	UBERON:feces                              1.0000e-03   2.1000e-02  1.0000e-03
	UBERON:glans penis                        8.4400e-01   1.0000e-03  4.4200e-01
	UBERON:hair                               6.0800e-01   1.0000e-03  9.3700e-01
	UBERON:labia minora                                    3.0000e-03  6.4400e-01
	UBERON:mouth                              5.0377e-04               1.0000e-03
	UBERON:nose                               6.4790e-01   6.2234e-06            
	UBERON:nostril                            9.4961e-01   3.5844e-06  5.1009e-01
	UBERON:nostrils                           8.9607e-01   4.1740e-06  6.0458e-01
	UBERON:skin of arm                        5.8305e-01   8.5397e-09  9.2451e-01
	UBERON:skin of finger                     3.7129e-02   3.8930e-04  5.4673e-04
	UBERON:skin of forearm                    7.1244e-03   1.3987e-02  4.6927e-05
	UBERON:tongue                             3.8368e-08   8.6365e-02  6.7536e-12
	UBERON:urine                              2.9047e-01   1.4300e-05  2.8182e-01
	UBERON:zone of skin of abdomen            8.7543e-01   3.7728e-09  3.5476e-01
	UBERON:zone of skin of foot               9.9897e-01   5.2351e-05  4.8594e-01
	UBERON:zone of skin of hand               1.0398e-02   5.0891e-04  7.9937e-06
	UBERON:zone of skin of head               4.2111e-01   3.4007e-08  7.1837e-01
	UBERON:zone of skin of knee               2.2705e-05   1.5893e-03  4.3363e-09
	UBERON:zone of skin of outer ear          5.2309e-02   1.1184e-11  3.9893e-02
									 UBERON:nostril UBERON:nostrils
	UBERON:ear canal                     5.4800e-01      6.4800e-01
	UBERON:feces                         1.0000e-03      1.0000e-03
	UBERON:glans penis                   7.1300e-01      6.4700e-01
	UBERON:hair                          4.8200e-01      5.3700e-01
	UBERON:labia minora                  9.5500e-01      8.9000e-01
	UBERON:mouth                         1.0000e-03      1.0000e-03
	UBERON:nose                          5.4100e-01      6.2100e-01
	UBERON:nostril                                       9.0600e-01
	UBERON:nostrils                      9.0932e-01                
	UBERON:skin of arm                   4.3891e-01      5.5419e-01
	UBERON:skin of finger                1.1299e-03      1.3746e-03
	UBERON:skin of forearm               5.3088e-05      7.9106e-05
	UBERON:tongue                        1.3483e-12      2.5948e-12
	UBERON:urine                         3.9401e-02      1.0177e-01
	UBERON:zone of skin of abdomen       7.6539e-01      6.7780e-01
	UBERON:zone of skin of foot          9.1018e-01      8.4083e-01
	UBERON:zone of skin of hand          8.4170e-06      2.9934e-05
	UBERON:zone of skin of head          1.9027e-01      2.9318e-01
	UBERON:zone of skin of knee          8.3639e-09      8.9881e-09
	UBERON:zone of skin of outer ear     3.9525e-04      3.3584e-03
									 UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                         8.8000e-01            3.0000e-03
	UBERON:feces                             1.0000e-03            1.0000e-03
	UBERON:glans penis                       3.2600e-01            3.6000e-02
	UBERON:hair                              8.6100e-01            1.0000e-03
	UBERON:labia minora                      6.0800e-01            4.2000e-02
	UBERON:mouth                             1.0000e-03            1.0000e-03
	UBERON:nose                              9.1800e-01            1.0000e-03
	UBERON:nostril                           4.3500e-01            2.0000e-03
	UBERON:nostrils                          5.6600e-01            2.0000e-03
	UBERON:skin of arm                                             1.0000e-03
	UBERON:skin of finger                    1.0833e-05                      
	UBERON:skin of forearm                   2.3364e-07            1.5539e-01
	UBERON:tongue                            3.8901e-17            2.2204e-09
	UBERON:urine                             1.2405e-01            3.2015e-05
	UBERON:zone of skin of abdomen           2.6311e-01            1.8168e-03
	UBERON:zone of skin of foot              4.0209e-01            3.9845e-03
	UBERON:zone of skin of hand              1.3756e-08            6.9669e-01
	UBERON:zone of skin of head              5.4531e-01            4.2034e-06
	UBERON:zone of skin of knee              3.7791e-13            5.5803e-03
	UBERON:zone of skin of outer ear         4.0678e-03            4.9576e-11
									 UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                             2.0000e-03    1.0000e-03
	UBERON:feces                                 1.0000e-03    8.0800e-01
	UBERON:glans penis                           1.0000e-02    1.0000e-03
	UBERON:hair                                  1.0000e-03    1.0000e-03
	UBERON:labia minora                          1.1000e-02    1.0000e-03
	UBERON:mouth                                 1.7000e-02    8.1000e-02
	UBERON:nose                                  1.0000e-03    1.0000e-03
	UBERON:nostril                               1.0000e-03    1.0000e-03
	UBERON:nostrils                              1.0000e-03    1.0000e-03
	UBERON:skin of arm                           1.0000e-03    1.0000e-03
	UBERON:skin of finger                        1.3700e-01    1.0000e-03
	UBERON:skin of forearm                                     1.0000e-03
	UBERON:tongue                                3.1911e-06              
	UBERON:urine                                 5.9563e-06    3.5196e-11
	UBERON:zone of skin of abdomen               4.2037e-05    4.4352e-17
	UBERON:zone of skin of foot                  3.0044e-04    1.1156e-10
	UBERON:zone of skin of hand                  2.0933e-01    6.9328e-10
	UBERON:zone of skin of head                  1.5240e-07    3.7477e-16
	UBERON:zone of skin of knee                  4.2781e-01    1.4055e-09
	UBERON:zone of skin of outer ear             1.9562e-12    1.9919e-22
									 UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                   3.3100e-01                     4.4600e-01
	UBERON:feces                       1.0000e-03                     1.0000e-03
	UBERON:glans penis                 1.7200e-01                     8.2900e-01
	UBERON:hair                        3.3700e-01                     2.7400e-01
	UBERON:labia minora                2.8100e-01                     8.7000e-01
	UBERON:mouth                       1.0000e-03                     1.0000e-03
	UBERON:nose                        2.6300e-01                     3.4700e-01
	UBERON:nostril                     4.5000e-02                     7.9200e-01
	UBERON:nostrils                    9.8000e-02                     6.8700e-01
	UBERON:skin of arm                 1.3500e-01                     2.6000e-01
	UBERON:skin of finger              1.0000e-03                     3.0000e-03
	UBERON:skin of forearm             1.0000e-03                     1.0000e-03
	UBERON:tongue                      1.0000e-03                     1.0000e-03
	UBERON:urine                                                      9.4000e-02
	UBERON:zone of skin of abdomen     1.0238e-01                               
	UBERON:zone of skin of foot        6.6989e-03                     8.7032e-01
	UBERON:zone of skin of hand        1.2493e-09                     3.5528e-04
	UBERON:zone of skin of head        2.2343e-01                     1.5466e-01
	UBERON:zone of skin of knee        1.8788e-09                     4.8468e-11
	UBERON:zone of skin of outer ear   6.3640e-01                     2.8249e-03
									 UBERON:zone of skin of foot
	UBERON:ear canal                                  5.0100e-01
	UBERON:feces                                      1.0000e-03
	UBERON:glans penis                                7.9400e-01
	UBERON:hair                                       4.6000e-01
	UBERON:labia minora                               1.0000e+00
	UBERON:mouth                                      2.0000e-03
	UBERON:nose                                       5.1000e-01
	UBERON:nostril                                    8.9500e-01
	UBERON:nostrils                                   8.2700e-01
	UBERON:skin of arm                                4.1100e-01
	UBERON:skin of finger                             4.0000e-03
	UBERON:skin of forearm                            1.0000e-03
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      9.0000e-03
	UBERON:zone of skin of abdomen                    8.6900e-01
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                       1.4322e-05
	UBERON:zone of skin of head                       1.4724e-01
	UBERON:zone of skin of knee                       4.3105e-07
	UBERON:zone of skin of outer ear                  2.9003e-05
									 UBERON:zone of skin of hand
	UBERON:ear canal                                  1.0000e-03
	UBERON:feces                                      1.0000e-03
	UBERON:glans penis                                1.3000e-02
	UBERON:hair                                       1.0000e-03
	UBERON:labia minora                               1.3000e-02
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       1.0000e-03
	UBERON:nostril                                    1.0000e-03
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                1.0000e-03
	UBERON:skin of finger                             7.2100e-01
	UBERON:skin of forearm                            2.0000e-01
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      1.0000e-03
	UBERON:zone of skin of abdomen                    2.0000e-03
	UBERON:zone of skin of foot                       1.0000e-03
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                       5.6656e-10
	UBERON:zone of skin of knee                       1.1382e-02
	UBERON:zone of skin of outer ear                  2.2930e-19
									 UBERON:zone of skin of head
	UBERON:ear canal                                  7.9800e-01
	UBERON:feces                                      1.0000e-03
	UBERON:glans penis                                2.1200e-01
	UBERON:hair                                       8.1600e-01
	UBERON:labia minora                               4.2300e-01
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       7.4400e-01
	UBERON:nostril                                    2.3000e-01
	UBERON:nostrils                                   2.9300e-01
	UBERON:skin of arm                                5.5600e-01
	UBERON:skin of finger                             1.0000e-03
	UBERON:skin of forearm                            1.0000e-03
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      2.3600e-01
	UBERON:zone of skin of abdomen                    1.5400e-01
	UBERON:zone of skin of foot                       1.5600e-01
	UBERON:zone of skin of hand                       1.0000e-03
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                       7.4644e-13
	UBERON:zone of skin of outer ear                  1.8548e-02
									 UBERON:zone of skin of knee
	UBERON:ear canal                                  1.0000e-03
	UBERON:feces                                      1.0000e-03
	UBERON:glans penis                                1.0000e-03
	UBERON:hair                                       1.0000e-03
	UBERON:labia minora                               1.0000e-03
	UBERON:mouth                                      3.0000e-03
	UBERON:nose                                       1.0000e-03
	UBERON:nostril                                    1.0000e-03
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                1.0000e-03
	UBERON:skin of finger                             3.0000e-03
	UBERON:skin of forearm                            3.9700e-01
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      1.0000e-03
	UBERON:zone of skin of abdomen                    1.0000e-03
	UBERON:zone of skin of foot                       1.0000e-03
	UBERON:zone of skin of hand                       6.0000e-03
	UBERON:zone of skin of head                       1.0000e-03
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                  8.6734e-20
									 UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.065
	UBERON:feces                                                0.001
	UBERON:glans penis                                          0.014
	UBERON:hair                                                 0.055
	UBERON:labia minora                                         0.046
	UBERON:mouth                                                0.001
	UBERON:nose                                                 0.041
	UBERON:nostril                                              0.002
	UBERON:nostrils                                             0.007
	UBERON:skin of arm                                          0.007
	UBERON:skin of finger                                       0.001
	UBERON:skin of forearm                                      0.001
	UBERON:tongue                                               0.001
	UBERON:urine                                                0.666
	UBERON:zone of skin of abdomen                              0.008
	UBERON:zone of skin of foot                                 0.001
	UBERON:zone of skin of hand                                 0.001
	UBERON:zone of skin of head                                 0.022
	UBERON:zone of skin of knee                                 0.001
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are significant.

Test 2
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	R --slave --args -d Whole\ Body/unweighted_unifrac_dm_shuffled_1.txt -m Whole\ Body/map.txt -c BODY_SITE -o betadisper_negative_1 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value Pr(>F)
	Groups     19 0.009105 0.00047923  0.9237 0.5532
	Residuals 565 0.293125 0.00051880               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
	Groups     19 0.009105 0.00047923 0.9237    999  0.572
	Residuals 565 0.293125 0.00051880                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
									 UBERON:ear canal UBERON:feces
	UBERON:ear canal                                      0.475000
	UBERON:feces                             0.464375             
	UBERON:glans penis                       0.479368     0.106057
	UBERON:hair                              0.608997     0.862392
	UBERON:labia minora                      0.732379     0.938807
	UBERON:mouth                             0.804742     0.610133
	UBERON:nose                              0.885097     0.277488
	UBERON:nostril                           0.297785     0.593561
	UBERON:nostrils                          0.456901     0.858758
	UBERON:skin of arm                       0.306992     0.636067
	UBERON:skin of finger                    0.502670     0.947465
	UBERON:skin of forearm                   0.572271     0.736379
	UBERON:tongue                            0.902904     0.202145
	UBERON:urine                             0.261478     0.559708
	UBERON:zone of skin of abdomen           0.431273     0.730046
	UBERON:zone of skin of foot              0.182610     0.342667
	UBERON:zone of skin of hand              0.176138     0.324192
	UBERON:zone of skin of head              0.623155     0.796835
	UBERON:zone of skin of knee              0.294462     0.557057
	UBERON:zone of skin of outer ear         0.624390     0.616991
									 UBERON:glans penis UBERON:hair
	UBERON:ear canal                           0.469000    0.611000
	UBERON:feces                               0.104000    0.862000
	UBERON:glans penis                                     0.035000
	UBERON:hair                                0.037969            
	UBERON:labia minora                        0.223167    0.967774
	UBERON:mouth                               0.143752    0.696830
	UBERON:nose                                0.245540    0.256096
	UBERON:nostril                             0.044311    0.521140
	UBERON:nostrils                            0.069082    0.723694
	UBERON:skin of arm                         0.034900    0.531579
	UBERON:skin of finger                      0.084757    0.897576
	UBERON:skin of forearm                     0.032060    0.899276
	UBERON:tongue                              0.377261    0.382896
	UBERON:urine                               0.046247    0.531347
	UBERON:zone of skin of abdomen             0.050148    0.574915
	UBERON:zone of skin of foot                0.039632    0.402329
	UBERON:zone of skin of hand                0.038464    0.388619
	UBERON:zone of skin of head                0.185124    0.965359
	UBERON:zone of skin of knee                0.071940    0.550065
	UBERON:zone of skin of outer ear           0.144469    0.854787
									 UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                            0.730000     0.785000    0.892000
	UBERON:feces                                0.946000     0.614000    0.287000
	UBERON:glans penis                          0.238000     0.156000    0.253000
	UBERON:hair                                 0.963000     0.722000    0.273000
	UBERON:labia minora                                      0.805000    0.479000
	UBERON:mouth                                0.787023                 0.562000
	UBERON:nose                                 0.474318     0.537941            
	UBERON:nostril                              0.710934     0.339893    0.117387
	UBERON:nostrils                             0.853577     0.508698    0.204524
	UBERON:skin of arm                          0.722026     0.341720    0.108680
	UBERON:skin of finger                       0.964381     0.629843    0.268176
	UBERON:skin of forearm                      0.901162     0.734133    0.246657
	UBERON:tongue                               0.555273     0.616823    0.961025
	UBERON:urine                                0.714368     0.336816    0.117546
	UBERON:zone of skin of abdomen              0.760404     0.414183    0.152136
	UBERON:zone of skin of foot                 0.611759     0.246003    0.083308
	UBERON:zone of skin of hand                 0.600483     0.236734    0.079593
	UBERON:zone of skin of head                 0.953164     0.789362    0.431541
	UBERON:zone of skin of knee                 0.724716     0.369625    0.151272
	UBERON:zone of skin of outer ear            0.877153     0.850946    0.426549
									 UBERON:nostril UBERON:nostrils
	UBERON:ear canal                       0.308000        0.490000
	UBERON:feces                           0.617000        0.879000
	UBERON:glans penis                     0.042000        0.078000
	UBERON:hair                            0.508000        0.722000
	UBERON:labia minora                    0.718000        0.841000
	UBERON:mouth                           0.347000        0.518000
	UBERON:nose                            0.118000        0.221000
	UBERON:nostril                                         0.793000
	UBERON:nostrils                        0.778314                
	UBERON:skin of arm                     0.952145        0.812703
	UBERON:skin of finger                  0.559947        0.812016
	UBERON:skin of forearm                 0.358656        0.589274
	UBERON:tongue                          0.090955        0.220261
	UBERON:urine                           0.982987        0.780808
	UBERON:zone of skin of abdomen         0.957168        0.851416
	UBERON:zone of skin of foot            0.775132        0.594466
	UBERON:zone of skin of hand            0.749243        0.575334
	UBERON:zone of skin of head            0.475962        0.714064
	UBERON:zone of skin of knee            0.975991        0.769935
	UBERON:zone of skin of outer ear       0.303128        0.562905
									 UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                           0.332000              0.532000
	UBERON:feces                               0.630000              0.947000
	UBERON:glans penis                         0.032000              0.079000
	UBERON:hair                                0.529000              0.913000
	UBERON:labia minora                        0.734000              0.953000
	UBERON:mouth                               0.349000              0.658000
	UBERON:nose                                0.113000              0.289000
	UBERON:nostril                             0.958000              0.577000
	UBERON:nostrils                            0.830000              0.810000
	UBERON:skin of arm                                               0.615000
	UBERON:skin of finger                      0.592914                      
	UBERON:skin of forearm                     0.371933              0.779539
	UBERON:tongue                              0.098531              0.248242
	UBERON:urine                               0.965699              0.541627
	UBERON:zone of skin of abdomen             0.994167              0.679596
	UBERON:zone of skin of foot                0.728657              0.358493
	UBERON:zone of skin of hand                0.703681              0.341737
	UBERON:zone of skin of head                0.508274              0.855446
	UBERON:zone of skin of knee                0.929763              0.547416
	UBERON:zone of skin of outer ear           0.336826              0.699842
									 UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                               0.593000      0.908000
	UBERON:feces                                   0.769000      0.210000
	UBERON:glans penis                             0.026000      0.395000
	UBERON:hair                                    0.905000      0.386000
	UBERON:labia minora                            0.915000      0.572000
	UBERON:mouth                                   0.736000      0.642000
	UBERON:nose                                    0.243000      0.956000
	UBERON:nostril                                 0.352000      0.095000
	UBERON:nostrils                                0.596000      0.245000
	UBERON:skin of arm                             0.368000      0.095000
	UBERON:skin of finger                          0.810000      0.263000
	UBERON:skin of forearm                                       0.310000
	UBERON:tongue                                  0.320265              
	UBERON:urine                                   0.357202      0.063404
	UBERON:zone of skin of abdomen                 0.447449      0.208862
	UBERON:zone of skin of foot                    0.229412      0.028155
	UBERON:zone of skin of hand                    0.217340      0.026109
	UBERON:zone of skin of head                    0.962142      0.369137
	UBERON:zone of skin of knee                    0.381047      0.082314
	UBERON:zone of skin of outer ear               0.908972      0.347478
									 UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                     0.270000                       0.447000
	UBERON:feces                         0.564000                       0.748000
	UBERON:glans penis                   0.042000                       0.048000
	UBERON:hair                          0.551000                       0.580000
	UBERON:labia minora                  0.716000                       0.762000
	UBERON:mouth                         0.333000                       0.427000
	UBERON:nose                          0.127000                       0.170000
	UBERON:nostril                       0.980000                       0.948000
	UBERON:nostrils                      0.766000                       0.858000
	UBERON:skin of arm                   0.967000                       0.991000
	UBERON:skin of finger                0.536000                       0.694000
	UBERON:skin of forearm               0.327000                       0.438000
	UBERON:tongue                        0.059000                       0.211000
	UBERON:urine                                                        0.974000
	UBERON:zone of skin of abdomen       0.968769                               
	UBERON:zone of skin of foot          0.719219                       0.796552
	UBERON:zone of skin of hand          0.689617                       0.778030
	UBERON:zone of skin of head          0.433225                       0.624664
	UBERON:zone of skin of knee          0.954878                       0.943106
	UBERON:zone of skin of outer ear     0.246144                       0.484679
									 UBERON:zone of skin of foot
	UBERON:ear canal                                    0.178000
	UBERON:feces                                        0.330000
	UBERON:glans penis                                  0.036000
	UBERON:hair                                         0.419000
	UBERON:labia minora                                 0.624000
	UBERON:mouth                                        0.242000
	UBERON:nose                                         0.102000
	UBERON:nostril                                      0.796000
	UBERON:nostrils                                     0.621000
	UBERON:skin of arm                                  0.735000
	UBERON:skin of finger                               0.366000
	UBERON:skin of forearm                              0.222000
	UBERON:tongue                                       0.034000
	UBERON:urine                                        0.697000
	UBERON:zone of skin of abdomen                      0.794000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                         0.965504
	UBERON:zone of skin of head                         0.267260
	UBERON:zone of skin of knee                         0.785651
	UBERON:zone of skin of outer ear                    0.110386
									 UBERON:zone of skin of hand
	UBERON:ear canal                                    0.177000
	UBERON:feces                                        0.314000
	UBERON:glans penis                                  0.041000
	UBERON:hair                                         0.405000
	UBERON:labia minora                                 0.613000
	UBERON:mouth                                        0.235000
	UBERON:nose                                         0.099000
	UBERON:nostril                                      0.744000
	UBERON:nostrils                                     0.595000
	UBERON:skin of arm                                  0.706000
	UBERON:skin of finger                               0.351000
	UBERON:skin of forearm                              0.195000
	UBERON:tongue                                       0.025000
	UBERON:urine                                        0.671000
	UBERON:zone of skin of abdomen                      0.780000
	UBERON:zone of skin of foot                         0.970000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                         0.253792
	UBERON:zone of skin of knee                         0.757453
	UBERON:zone of skin of outer ear                    0.101739
									 UBERON:zone of skin of head
	UBERON:ear canal                                    0.629000
	UBERON:feces                                        0.809000
	UBERON:glans penis                                  0.185000
	UBERON:hair                                         0.957000
	UBERON:labia minora                                 0.959000
	UBERON:mouth                                        0.800000
	UBERON:nose                                         0.466000
	UBERON:nostril                                      0.496000
	UBERON:nostrils                                     0.706000
	UBERON:skin of arm                                  0.529000
	UBERON:skin of finger                               0.868000
	UBERON:skin of forearm                              0.972000
	UBERON:tongue                                       0.383000
	UBERON:urine                                        0.428000
	UBERON:zone of skin of abdomen                      0.620000
	UBERON:zone of skin of foot                         0.290000
	UBERON:zone of skin of hand                         0.235000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                         0.445794
	UBERON:zone of skin of outer ear                    0.873684
									 UBERON:zone of skin of knee
	UBERON:ear canal                                    0.293000
	UBERON:feces                                        0.561000
	UBERON:glans penis                                  0.069000
	UBERON:hair                                         0.565000
	UBERON:labia minora                                 0.702000
	UBERON:mouth                                        0.367000
	UBERON:nose                                         0.161000
	UBERON:nostril                                      0.979000
	UBERON:nostrils                                     0.778000
	UBERON:skin of arm                                  0.929000
	UBERON:skin of finger                               0.589000
	UBERON:skin of forearm                              0.378000
	UBERON:tongue                                       0.079000
	UBERON:urine                                        0.953000
	UBERON:zone of skin of abdomen                      0.936000
	UBERON:zone of skin of foot                         0.797000
	UBERON:zone of skin of hand                         0.761000
	UBERON:zone of skin of head                         0.465000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                    0.261342
									 UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.639
	UBERON:feces                                                0.628
	UBERON:glans penis                                          0.138
	UBERON:hair                                                 0.859
	UBERON:labia minora                                         0.866
	UBERON:mouth                                                0.865
	UBERON:nose                                                 0.430
	UBERON:nostril                                              0.318
	UBERON:nostrils                                             0.578
	UBERON:skin of arm                                          0.354
	UBERON:skin of finger                                       0.715
	UBERON:skin of forearm                                      0.915
	UBERON:tongue                                               0.363
	UBERON:urine                                                0.245
	UBERON:zone of skin of abdomen                              0.489
	UBERON:zone of skin of foot                                 0.138
	UBERON:zone of skin of hand                                 0.111
	UBERON:zone of skin of head                                 0.884
	UBERON:zone of skin of knee                                 0.287
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are insignificant.
	
Test 3
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	R --slave --args -d Whole\ Body/unweighted_unifrac_dm_shuffled_2.txt -m Whole\ Body/map.txt -c BODY_SITE -o betadisper_negative_2 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value Pr(>F)
	Groups     19 0.007599 0.00039994  0.7714 0.7423
	Residuals 565 0.292927 0.00051845               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
	Groups     19 0.007599 0.00039994 0.7714    999  0.741
	Residuals 565 0.292927 0.00051845                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
									 UBERON:ear canal UBERON:feces
	UBERON:ear canal                                      0.802000
	UBERON:feces                             0.809460             
	UBERON:glans penis                       0.630392     0.439281
	UBERON:hair                              0.432323     0.232982
	UBERON:labia minora                      0.839904     0.683421
	UBERON:mouth                             0.831786     0.606633
	UBERON:nose                              0.890293     0.915637
	UBERON:nostril                           0.359574     0.142120
	UBERON:nostrils                          0.935456     0.691673
	UBERON:skin of arm                       0.811866     0.959382
	UBERON:skin of finger                    0.775908     0.473111
	UBERON:skin of forearm                   0.688989     0.397107
	UBERON:tongue                            0.396479     0.139000
	UBERON:urine                             0.844537     0.905967
	UBERON:zone of skin of abdomen           0.773034     0.942461
	UBERON:zone of skin of foot              0.432269     0.109658
	UBERON:zone of skin of hand              0.514257     0.162353
	UBERON:zone of skin of head              0.946979     0.804400
	UBERON:zone of skin of knee              0.772119     0.994388
	UBERON:zone of skin of outer ear         0.735179     0.901716
									 UBERON:glans penis UBERON:hair
	UBERON:ear canal                           0.623000    0.427000
	UBERON:feces                               0.468000    0.234000
	UBERON:glans penis                                     0.951000
	UBERON:hair                                0.950590            
	UBERON:labia minora                        0.779727    0.607250
	UBERON:mouth                               0.625495    0.397652
	UBERON:nose                                0.462580    0.231632
	UBERON:nostril                             0.959928    0.842376
	UBERON:nostrils                            0.547000    0.318343
	UBERON:skin of arm                         0.354905    0.143593
	UBERON:skin of finger                      0.693789    0.507098
	UBERON:skin of forearm                     0.739032    0.558205
	UBERON:tongue                              0.984249    0.894538
	UBERON:urine                               0.399236    0.189895
	UBERON:zone of skin of abdomen             0.381229    0.154025
	UBERON:zone of skin of foot                0.946796    0.847630
	UBERON:zone of skin of hand                0.815699    0.659079
	UBERON:zone of skin of head                0.510606    0.299189
	UBERON:zone of skin of knee                0.331318    0.131227
	UBERON:zone of skin of outer ear           0.378568    0.178251
									 UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                            0.849000     0.831000    0.891000
	UBERON:feces                                0.713000     0.630000    0.928000
	UBERON:glans penis                          0.788000     0.643000    0.486000
	UBERON:hair                                 0.632000     0.394000    0.219000
	UBERON:labia minora                                      0.938000    0.693000
	UBERON:mouth                                0.922316                 0.638000
	UBERON:nose                                 0.668968     0.625099            
	UBERON:nostril                              0.604914     0.373313    0.164165
	UBERON:nostrils                             0.823398     0.846372    0.765859
	UBERON:skin of arm                          0.597616     0.525545    0.931848
	UBERON:skin of finger                       0.987205     0.933117    0.613199
	UBERON:skin of forearm                      0.942072     0.822956    0.505294
	UBERON:tongue                               0.669089     0.451496    0.231599
	UBERON:urine                                0.672636     0.606852    0.975367
	UBERON:zone of skin of abdomen              0.523712     0.451234    0.834166
	UBERON:zone of skin of foot                 0.740107     0.541309    0.295130
	UBERON:zone of skin of hand                 0.831319     0.657473    0.352920
	UBERON:zone of skin of head                 0.770579     0.739934    0.926256
	UBERON:zone of skin of knee                 0.593879     0.507705    0.897088
	UBERON:zone of skin of outer ear            0.626102     0.528088    0.840951
									 UBERON:nostril UBERON:nostrils
	UBERON:ear canal                       0.351000        0.941000
	UBERON:feces                           0.126000        0.681000
	UBERON:glans penis                     0.973000        0.584000
	UBERON:hair                            0.841000        0.323000
	UBERON:labia minora                    0.643000        0.846000
	UBERON:mouth                           0.389000        0.836000
	UBERON:nose                            0.151000        0.762000
	UBERON:nostril                                         0.262000
	UBERON:nostrils                        0.263851                
	UBERON:skin of arm                     0.082189        0.659537
	UBERON:skin of finger                  0.476824        0.791254
	UBERON:skin of forearm                 0.558492        0.675255
	UBERON:tongue                          0.963563        0.323426
	UBERON:urine                           0.112238        0.718188
	UBERON:zone of skin of abdomen         0.095500        0.600444
	UBERON:zone of skin of foot            0.958232        0.390917
	UBERON:zone of skin of hand            0.707515        0.490054
	UBERON:zone of skin of head            0.215798        0.853394
	UBERON:zone of skin of knee            0.067361        0.621083
	UBERON:zone of skin of outer ear       0.094194        0.604039
									 UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                           0.820000              0.767000
	UBERON:feces                               0.964000              0.477000
	UBERON:glans penis                         0.390000              0.701000
	UBERON:hair                                0.144000              0.516000
	UBERON:labia minora                        0.631000              0.992000
	UBERON:mouth                               0.550000              0.949000
	UBERON:nose                                0.928000              0.605000
	UBERON:nostril                             0.078000              0.464000
	UBERON:nostrils                            0.671000              0.793000
	UBERON:skin of arm                                               0.471000
	UBERON:skin of finger                      0.470137                      
	UBERON:skin of forearm                     0.363554              0.891104
	UBERON:tongue                              0.114823              0.483947
	UBERON:urine                               0.951308              0.480896
	UBERON:zone of skin of abdomen             0.882940              0.498276
	UBERON:zone of skin of foot                0.143939              0.516905
	UBERON:zone of skin of hand                0.191491              0.676127
	UBERON:zone of skin of head                0.837001              0.637791
	UBERON:zone of skin of knee                0.955871              0.400612
	UBERON:zone of skin of outer ear           0.864564              0.376541
									 UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                               0.663000      0.397000
	UBERON:feces                                   0.380000      0.126000
	UBERON:glans penis                             0.747000      0.985000
	UBERON:hair                                    0.533000      0.884000
	UBERON:labia minora                            0.941000      0.678000
	UBERON:mouth                                   0.827000      0.463000
	UBERON:nose                                    0.516000      0.219000
	UBERON:nostril                                 0.536000      0.965000
	UBERON:nostrils                                0.674000      0.334000
	UBERON:skin of arm                             0.356000      0.126000
	UBERON:skin of finger                          0.887000      0.501000
	UBERON:skin of forearm                                       0.547000
	UBERON:tongue                                  0.570831              
	UBERON:urine                                   0.388207      0.116342
	UBERON:zone of skin of abdomen                 0.391266      0.168348
	UBERON:zone of skin of foot                    0.627143      0.924260
	UBERON:zone of skin of hand                    0.804528      0.678936
	UBERON:zone of skin of head                    0.541173      0.224593
	UBERON:zone of skin of knee                    0.308430      0.078869
	UBERON:zone of skin of outer ear               0.310128      0.087866
									 UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                     0.865000                       0.778000
	UBERON:feces                         0.913000                       0.946000
	UBERON:glans penis                   0.419000                       0.382000
	UBERON:hair                          0.187000                       0.149000
	UBERON:labia minora                  0.701000                       0.560000
	UBERON:mouth                         0.623000                       0.453000
	UBERON:nose                          0.974000                       0.838000
	UBERON:nostril                       0.118000                       0.108000
	UBERON:nostrils                      0.719000                       0.611000
	UBERON:skin of arm                   0.953000                       0.879000
	UBERON:skin of finger                0.490000                       0.497000
	UBERON:skin of forearm               0.395000                       0.371000
	UBERON:tongue                        0.111000                       0.153000
	UBERON:urine                                                        0.870000
	UBERON:zone of skin of abdomen       0.860518                               
	UBERON:zone of skin of foot          0.103726                       0.236692
	UBERON:zone of skin of hand          0.159310                       0.272320
	UBERON:zone of skin of head          0.868728                       0.786743
	UBERON:zone of skin of knee          0.898674                       0.920465
	UBERON:zone of skin of outer ear     0.791662                       0.993316
									 UBERON:zone of skin of foot
	UBERON:ear canal                                    0.420000
	UBERON:feces                                        0.105000
	UBERON:glans penis                                  0.951000
	UBERON:hair                                         0.832000
	UBERON:labia minora                                 0.750000
	UBERON:mouth                                        0.545000
	UBERON:nose                                         0.294000
	UBERON:nostril                                      0.959000
	UBERON:nostrils                                     0.366000
	UBERON:skin of arm                                  0.154000
	UBERON:skin of finger                               0.532000
	UBERON:skin of forearm                              0.625000
	UBERON:tongue                                       0.921000
	UBERON:urine                                        0.104000
	UBERON:zone of skin of abdomen                      0.243000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                         0.724254
	UBERON:zone of skin of head                         0.219958
	UBERON:zone of skin of knee                         0.079400
	UBERON:zone of skin of outer ear                    0.057227
									 UBERON:zone of skin of hand
	UBERON:ear canal                                    0.516000
	UBERON:feces                                        0.158000
	UBERON:glans penis                                  0.809000
	UBERON:hair                                         0.638000
	UBERON:labia minora                                 0.862000
	UBERON:mouth                                        0.657000
	UBERON:nose                                         0.341000
	UBERON:nostril                                      0.685000
	UBERON:nostrils                                     0.468000
	UBERON:skin of arm                                  0.199000
	UBERON:skin of finger                               0.657000
	UBERON:skin of forearm                              0.809000
	UBERON:tongue                                       0.669000
	UBERON:urine                                        0.176000
	UBERON:zone of skin of abdomen                      0.278000
	UBERON:zone of skin of foot                         0.731000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                         0.302969
	UBERON:zone of skin of knee                         0.119274
	UBERON:zone of skin of outer ear                    0.093463
									 UBERON:zone of skin of head
	UBERON:ear canal                                    0.952000
	UBERON:feces                                        0.799000
	UBERON:glans penis                                  0.531000
	UBERON:hair                                         0.297000
	UBERON:labia minora                                 0.794000
	UBERON:mouth                                        0.753000
	UBERON:nose                                         0.910000
	UBERON:nostril                                      0.209000
	UBERON:nostrils                                     0.850000
	UBERON:skin of arm                                  0.838000
	UBERON:skin of finger                               0.619000
	UBERON:skin of forearm                              0.532000
	UBERON:tongue                                       0.231000
	UBERON:urine                                        0.866000
	UBERON:zone of skin of abdomen                      0.775000
	UBERON:zone of skin of foot                         0.232000
	UBERON:zone of skin of hand                         0.285000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                         0.779349
	UBERON:zone of skin of outer ear                    0.700846
									 UBERON:zone of skin of knee
	UBERON:ear canal                                    0.755000
	UBERON:feces                                        0.993000
	UBERON:glans penis                                  0.332000
	UBERON:hair                                         0.129000
	UBERON:labia minora                                 0.632000
	UBERON:mouth                                        0.532000
	UBERON:nose                                         0.898000
	UBERON:nostril                                      0.086000
	UBERON:nostrils                                     0.617000
	UBERON:skin of arm                                  0.956000
	UBERON:skin of finger                               0.378000
	UBERON:skin of forearm                              0.313000
	UBERON:tongue                                       0.082000
	UBERON:urine                                        0.902000
	UBERON:zone of skin of abdomen                      0.915000
	UBERON:zone of skin of foot                         0.100000
	UBERON:zone of skin of hand                         0.112000
	UBERON:zone of skin of head                         0.793000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                    0.886053
									 UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.715
	UBERON:feces                                                0.906
	UBERON:glans penis                                          0.390
	UBERON:hair                                                 0.172
	UBERON:labia minora                                         0.671
	UBERON:mouth                                                0.541
	UBERON:nose                                                 0.817
	UBERON:nostril                                              0.093
	UBERON:nostrils                                             0.601
	UBERON:skin of arm                                          0.862
	UBERON:skin of finger                                       0.371
	UBERON:skin of forearm                                      0.306
	UBERON:tongue                                               0.073
	UBERON:urine                                                0.780
	UBERON:zone of skin of abdomen                              0.992
	UBERON:zone of skin of foot                                 0.057
	UBERON:zone of skin of hand                                 0.083
	UBERON:zone of skin of head                                 0.657
	UBERON:zone of skin of knee                                 0.897
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are insignificant.
	
Test 4
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	R --slave --args -d Whole\ Body/unweighted_unifrac_dm_shuffled_3.txt -m Whole\ Body/map.txt -c BODY_SITE -o betadisper_negative_3 < betadisper.r

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
			   Df   Sum Sq    Mean Sq F value Pr(>F)
	Groups     19 0.008864 0.00046655  0.8801 0.6084
	Residuals 565 0.299502 0.00053009               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
			   Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
	Groups     19 0.008864 0.00046655 0.8801    999  0.618
	Residuals 565 0.299502 0.00053009                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
									 UBERON:ear canal UBERON:feces
	UBERON:ear canal                                     0.6690000
	UBERON:feces                            0.6485757             
	UBERON:glans penis                      0.3513263    0.1549040
	UBERON:hair                             0.7101255    0.9495634
	UBERON:labia minora                     0.8684986    0.5943062
	UBERON:mouth                            0.4049724    0.5495694
	UBERON:nose                             0.2736159    0.0662598
	UBERON:nostril                          0.9403975    0.5525088
	UBERON:nostrils                         0.6447169    0.9233151
	UBERON:skin of arm                      0.5663604    0.7924728
	UBERON:skin of finger                   0.8889145    0.6659402
	UBERON:skin of forearm                  0.5725727    0.7805681
	UBERON:tongue                           0.2453426    0.4395189
	UBERON:urine                            0.5582436    0.8742373
	UBERON:zone of skin of abdomen          0.8175508    0.4572090
	UBERON:zone of skin of foot             0.3333797    0.4738408
	UBERON:zone of skin of hand             0.7928409    0.6740538
	UBERON:zone of skin of head             0.2696213    0.3498215
	UBERON:zone of skin of knee             0.2875977    0.4046836
	UBERON:zone of skin of outer ear        0.4078144    0.5771957
									 UBERON:glans penis UBERON:hair
	UBERON:ear canal                          0.3610000   0.7200000
	UBERON:feces                              0.1510000   0.9540000
	UBERON:glans penis                                    0.2670000
	UBERON:hair                               0.2540041            
	UBERON:labia minora                       0.5314495   0.6698494
	UBERON:mouth                              0.0746252   0.7011375
	UBERON:nose                               0.9749716   0.1590738
	UBERON:nostril                            0.1084239   0.6111217
	UBERON:nostrils                           0.1567516   0.9872483
	UBERON:skin of arm                        0.1591828   0.9017453
	UBERON:skin of finger                     0.2127800   0.7183918
	UBERON:skin of forearm                    0.1782171   0.8927804
	UBERON:tongue                             0.0095797   0.6246373
	UBERON:urine                              0.1125464   0.9694011
	UBERON:zone of skin of abdomen            0.4118793   0.5532235
	UBERON:zone of skin of foot               0.0572719   0.6947518
	UBERON:zone of skin of hand               0.1448509   0.7252143
	UBERON:zone of skin of head               0.0452218   0.5645821
	UBERON:zone of skin of knee               0.0417863   0.6202207
	UBERON:zone of skin of outer ear          0.0878269   0.7697631
									 UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                           0.8480000    0.3980000   0.2550000
	UBERON:feces                               0.5750000    0.5770000   0.0620000
	UBERON:glans penis                         0.5380000    0.0850000   0.9770000
	UBERON:hair                                0.6500000    0.6870000   0.1760000
	UBERON:labia minora                                     0.4030000   0.4750000
	UBERON:mouth                               0.4127632                0.0440000
	UBERON:nose                                0.4983522    0.0373902            
	UBERON:nostril                             0.7439355    0.2192758   0.0701766
	UBERON:nostrils                            0.5954483    0.6587076   0.0875841
	UBERON:skin of arm                         0.5512837    0.7455973   0.0735994
	UBERON:skin of finger                      0.7535290    0.3538976   0.1255757
	UBERON:skin of forearm                     0.5630486    0.7706877   0.0850692
	UBERON:tongue                              0.2407563    0.9139570   0.0027092
	UBERON:urine                               0.5234979    0.6048419   0.0410132
	UBERON:zone of skin of abdomen             0.9977165    0.2639228   0.3645394
	UBERON:zone of skin of foot                0.3614583    0.8646884   0.0127110
	UBERON:zone of skin of hand                0.6661219    0.3156342   0.0594514
	UBERON:zone of skin of head                0.3091188    0.9182877   0.0117821
	UBERON:zone of skin of knee                0.3178699    0.9725923   0.0098219
	UBERON:zone of skin of outer ear           0.4261082    0.8133024   0.0248113
									 UBERON:nostril UBERON:nostrils
	UBERON:ear canal                      0.9410000       0.6210000
	UBERON:feces                          0.5580000       0.9240000
	UBERON:glans penis                    0.1150000       0.1850000
	UBERON:hair                           0.6330000       0.9930000
	UBERON:labia minora                   0.7470000       0.5990000
	UBERON:mouth                          0.2330000       0.6510000
	UBERON:nose                           0.0620000       0.0890000
	UBERON:nostril                                        0.5460000
	UBERON:nostrils                       0.5260959                
	UBERON:skin of arm                    0.4278797       0.9002309
	UBERON:skin of finger                 0.9048184       0.6546939
	UBERON:skin of forearm                0.4343749       0.8892125
	UBERON:tongue                         0.0929435       0.5783329
	UBERON:urine                          0.4324461       0.9833172
	UBERON:zone of skin of abdomen        0.6559662       0.4646113
	UBERON:zone of skin of foot           0.1785103       0.6676765
	UBERON:zone of skin of hand           0.7766800       0.6634468
	UBERON:zone of skin of head           0.1129581       0.5133439
	UBERON:zone of skin of knee           0.1313337       0.5786226
	UBERON:zone of skin of outer ear      0.2527929       0.7525041
									 UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                          0.5570000             0.8900000
	UBERON:feces                              0.7870000             0.6660000
	UBERON:glans penis                        0.1730000             0.2270000
	UBERON:hair                               0.9080000             0.7240000
	UBERON:labia minora                       0.5440000             0.7470000
	UBERON:mouth                              0.7420000             0.3650000
	UBERON:nose                               0.0890000             0.1270000
	UBERON:nostril                            0.4350000             0.8980000
	UBERON:nostrils                           0.9060000             0.6580000
	UBERON:skin of arm                                              0.5370000
	UBERON:skin of finger                     0.5401149                      
	UBERON:skin of forearm                    0.9834464             0.5415090
	UBERON:tongue                             0.7048538             0.2052761
	UBERON:urine                              0.8857238             0.5497821
	UBERON:zone of skin of abdomen            0.4075976             0.6628269
	UBERON:zone of skin of foot               0.7686770             0.2607727
	UBERON:zone of skin of hand               0.5072174             0.8960123
	UBERON:zone of skin of head               0.5902142             0.1948869
	UBERON:zone of skin of knee               0.6708797             0.2176967
	UBERON:zone of skin of outer ear          0.8562840             0.3427562
									 UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                              0.5700000     0.2580000
	UBERON:feces                                  0.7950000     0.4330000
	UBERON:glans penis                            0.1850000     0.0150000
	UBERON:hair                                   0.8960000     0.6510000
	UBERON:labia minora                           0.5360000     0.2260000
	UBERON:mouth                                  0.7800000     0.9160000
	UBERON:nose                                   0.0810000     0.0020000
	UBERON:nostril                                0.4460000     0.0860000
	UBERON:nostrils                               0.8780000     0.5640000
	UBERON:skin of arm                            0.9830000     0.6840000
	UBERON:skin of finger                         0.5250000     0.2100000
	UBERON:skin of forearm                                      0.7340000
	UBERON:tongue                                 0.7380587              
	UBERON:urine                                  0.8693053     0.5196464
	UBERON:zone of skin of abdomen                0.4206314     0.1134016
	UBERON:zone of skin of foot                   0.7976490     0.9088106
	UBERON:zone of skin of hand                   0.5029245     0.1752866
	UBERON:zone of skin of head                   0.6212949     0.7770731
	UBERON:zone of skin of knee                   0.7017056     0.9303723
	UBERON:zone of skin of outer ear              0.8815777     0.8230433
									 UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                    0.5640000                      0.8160000
	UBERON:feces                        0.8670000                      0.4650000
	UBERON:glans penis                  0.1200000                      0.4000000
	UBERON:hair                         0.9730000                      0.5590000
	UBERON:labia minora                 0.5200000                      0.9980000
	UBERON:mouth                        0.6020000                      0.2570000
	UBERON:nose                         0.0410000                      0.3600000
	UBERON:nostril                      0.4670000                      0.6530000
	UBERON:nostrils                     0.9860000                      0.4680000
	UBERON:skin of arm                  0.8740000                      0.3950000
	UBERON:skin of finger               0.5240000                      0.6500000
	UBERON:skin of forearm              0.8460000                      0.4380000
	UBERON:tongue                       0.5010000                      0.1100000
	UBERON:urine                                                       0.3900000
	UBERON:zone of skin of abdomen      0.3742275                               
	UBERON:zone of skin of foot         0.5699519                      0.2030125
	UBERON:zone of skin of hand         0.5296818                      0.5456379
	UBERON:zone of skin of head         0.4059798                      0.1609112
	UBERON:zone of skin of knee         0.4786967                      0.1675283
	UBERON:zone of skin of outer ear    0.6793253                      0.2661091
									 UBERON:zone of skin of foot
	UBERON:ear canal                                   0.3350000
	UBERON:feces                                       0.4430000
	UBERON:glans penis                                 0.0580000
	UBERON:hair                                        0.7000000
	UBERON:labia minora                                0.3520000
	UBERON:mouth                                       0.8700000
	UBERON:nose                                        0.0110000
	UBERON:nostril                                     0.1990000
	UBERON:nostrils                                    0.6890000
	UBERON:skin of arm                                 0.7680000
	UBERON:skin of finger                              0.2550000
	UBERON:skin of forearm                             0.7850000
	UBERON:tongue                                      0.9030000
	UBERON:urine                                       0.5660000
	UBERON:zone of skin of abdomen                     0.1890000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                        0.1816522
	UBERON:zone of skin of head                        0.6989050
	UBERON:zone of skin of knee                        0.8383950
	UBERON:zone of skin of outer ear                   0.8953861
									 UBERON:zone of skin of hand
	UBERON:ear canal                                   0.8010000
	UBERON:feces                                       0.7040000
	UBERON:glans penis                                 0.1610000
	UBERON:hair                                        0.7560000
	UBERON:labia minora                                0.6440000
	UBERON:mouth                                       0.3140000
	UBERON:nose                                        0.0700000
	UBERON:nostril                                     0.7600000
	UBERON:nostrils                                    0.6760000
	UBERON:skin of arm                                 0.5100000
	UBERON:skin of finger                              0.8970000
	UBERON:skin of forearm                             0.5190000
	UBERON:tongue                                      0.1670000
	UBERON:urine                                       0.5300000
	UBERON:zone of skin of abdomen                     0.5520000
	UBERON:zone of skin of foot                        0.1840000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                        0.1321384
	UBERON:zone of skin of knee                        0.1537181
	UBERON:zone of skin of outer ear                   0.2638695
									 UBERON:zone of skin of head
	UBERON:ear canal                                   0.2510000
	UBERON:feces                                       0.3540000
	UBERON:glans penis                                 0.0460000
	UBERON:hair                                        0.5820000
	UBERON:labia minora                                0.3070000
	UBERON:mouth                                       0.9260000
	UBERON:nose                                        0.0130000
	UBERON:nostril                                     0.1160000
	UBERON:nostrils                                    0.5050000
	UBERON:skin of arm                                 0.6040000
	UBERON:skin of finger                              0.2000000
	UBERON:skin of forearm                             0.6380000
	UBERON:tongue                                      0.7720000
	UBERON:urine                                       0.3980000
	UBERON:zone of skin of abdomen                     0.1400000
	UBERON:zone of skin of foot                        0.6830000
	UBERON:zone of skin of hand                        0.1430000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                        0.8496472
	UBERON:zone of skin of outer ear                   0.6400064
									 UBERON:zone of skin of knee
	UBERON:ear canal                                   0.2780000
	UBERON:feces                                       0.4210000
	UBERON:glans penis                                 0.0410000
	UBERON:hair                                        0.6360000
	UBERON:labia minora                                0.3030000
	UBERON:mouth                                       0.9730000
	UBERON:nose                                        0.0100000
	UBERON:nostril                                     0.1330000
	UBERON:nostrils                                    0.5950000
	UBERON:skin of arm                                 0.6660000
	UBERON:skin of finger                              0.2080000
	UBERON:skin of forearm                             0.6960000
	UBERON:tongue                                      0.9210000
	UBERON:urine                                       0.4620000
	UBERON:zone of skin of abdomen                     0.1590000
	UBERON:zone of skin of foot                        0.8410000
	UBERON:zone of skin of hand                        0.1520000
	UBERON:zone of skin of head                        0.8280000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                   0.7577955
									 UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.397
	UBERON:feces                                                0.603
	UBERON:glans penis                                          0.092
	UBERON:hair                                                 0.778
	UBERON:labia minora                                         0.414
	UBERON:mouth                                                0.798
	UBERON:nose                                                 0.025
	UBERON:nostril                                              0.265
	UBERON:nostrils                                             0.769
	UBERON:skin of arm                                          0.863
	UBERON:skin of finger                                       0.338
	UBERON:skin of forearm                                      0.876
	UBERON:tongue                                               0.808
	UBERON:urine                                                0.704
	UBERON:zone of skin of abdomen                              0.250
	UBERON:zone of skin of foot                                 0.897
	UBERON:zone of skin of hand                                 0.251
	UBERON:zone of skin of head                                 0.635
	UBERON:zone of skin of knee                                 0.767
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are insignificant.

References
----------

[1]
http://www.stat.auckland.ac.nz/~mja/Programs.htm

[2]
http://www.stat.auckland.ac.nz/~mja/prog/PERMDISP_UserNotes.pdf

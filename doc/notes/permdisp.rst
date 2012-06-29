=============================================================
Permutational Analysis of Multivariate Dispersions (PERMDISP)
=============================================================

Synopsis
--------

PERMDISP is a procedure for the analysis of multivariate homogeneity of group dispersions (variances). Permutations can be utilized to measure the dissimilatities between groups.

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

	compare_categories.py --method permdisp -i datasets/glen_canyon/unweighted_unifrac_dm.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o glen_canyon_positive_permdisp_results

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		  Df   Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.003154 0.0031536   1.762 0.1877
	Residuals 92 0.164659 0.0017898               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		  Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
	Groups     1 0.003154 0.0031536 1.762    999  0.186
	Residuals 92 0.164659 0.0017898                    

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		 No   Yes
	No          0.195
	Yes 0.18766          
	
The p-value indicates that the results are more significant than the negative tests.
	
Test 2
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o glen_canyon_negative_permdisp_results_1

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		  Df  Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.00162 0.0016154  0.1863  0.667
	Residuals 92 0.79769 0.0086705               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		  Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.00162 0.0016154 0.1863    999  0.666
	Residuals 92 0.79769 0.0086705                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		 No   Yes
	No          0.666
	Yes 0.66701      

The p-value indicates that the results are insignificant.
	
Test 3
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o glen_canyon_negative_permdisp_results_2

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		  Df  Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.00004 0.0000404  0.0048 0.9452
	Residuals 92 0.78170 0.0084967               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		  Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.00004 0.0000404 0.0048    999  0.942
	Residuals 92 0.78170 0.0084967                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		 No  Yes
	No          0.94
	Yes 0.94518        

The p-value indicates that the results insignificant.
	
Test 4
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the CurrentlyWet category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet -o glen_canyon_negative_permdisp_results_3
	
**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		  Df  Sum Sq   Mean Sq F value Pr(>F)
	Groups     1 0.00516 0.0051590  0.5948 0.4425
	Residuals 92 0.79792 0.0086731               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		  Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     1 0.00516 0.0051590 0.5948    999  0.443
	Residuals 92 0.79792 0.0086731                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		 No   Yes
	No          0.445
	Yes 0.44254             
	
The p-value indicates that the results are insignificant.
	
Keyboard
^^^^^^^^

Test 1
~~~~~~

**Description:**

This test uses the original unweighted unifrac distance matrix and the HOST_SUBJECT_ID category as a positive control.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o keyboard_positive_permdisp_results

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table
	
	Response: Distances
		   Df  Sum Sq  Mean Sq F value    Pr(>F)    
	Groups     10 0.81031 0.081031  36.751 < 2.2e-16 ***
	Residuals 104 0.22931 0.002205                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
	Groups     10 0.81031 0.081031 36.751    999  0.001 ***
	Residuals 104 0.22931 0.002205                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		   F1 L1 L3         M1         M2         M3         M9 R1 U1 U2 U3
	F1                  3.9800e-01 1.9100e-01 2.6700e-01 9.1000e-01            
	L1                                                                         
	L3                                                                         
	M1 3.7767e-01                  1.6100e-01 1.8000e-02 1.2000e-01            
	M2 2.1166e-01       2.1787e-01            1.0000e-03 4.0000e-03            
	M3 3.1179e-01       1.4740e-03 1.0126e-08            8.2000e-02            
	M9 9.3542e-01       1.1021e-01 3.5387e-03 6.7247e-02                       
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

	compare_categories.py --method permdisp -i datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o keyboard_negative_permdisp_results_1

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df  Sum Sq Mean Sq F value    Pr(>F)    
	Groups     10 1.07966 0.10797  41.200 < 2.2e-16 ***
	Residuals 104 0.27253 0.00262                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
	Groups     10 1.07966 0.10797 41.200    999  0.001 ***
	Residuals 104 0.27253 0.00262                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		   F1 L1 L3         M1         M2         M3         M9 R1 U1 U2 U3
	F1                  3.2200e-01 1.0000e-02 2.6000e-02 2.0000e-02            
	L1                                                                         
	L3                                                                         
	M1 3.2691e-01                  5.0000e-03 9.0000e-03 8.0000e-03            
	M2 2.6932e-03       3.8059e-05            5.5900e-01 4.0900e-01            
	M3 1.2583e-02       4.0305e-04 5.5574e-01            8.1100e-01            
	M9 9.7747e-03       1.8902e-04 3.7458e-01 7.9932e-01                       
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

	compare_categories.py --method permdisp -i datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o keyboard_negative_permdisp_results_2

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df  Sum Sq  Mean Sq F value    Pr(>F)    
	Groups     10 1.08893 0.108893  42.434 < 2.2e-16 ***
	Residuals 104 0.26688 0.002566                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
	Groups     10 1.08893 0.108893 42.434    999  0.001 ***
	Residuals 104 0.26688 0.002566                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		   F1 L1 L3         M1         M2         M3         M9 R1 U1 U2 U3
	F1                  1.9300e-01 1.2600e-01 5.7000e-02 1.4400e-01            
	L1                                                                         
	L3                                                                         
	M1 2.0791e-01                  5.0000e-03 3.0000e-03 7.0000e-03            
	M2 1.4750e-01       1.9806e-06            2.0200e-01 6.0000e-01            
	M3 5.9335e-02       3.8031e-07 1.9174e-01            5.7100e-01            
	M9 1.7403e-01       5.8180e-05 5.9153e-01 5.6234e-01                       
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

	compare_categories.py --method permdisp -i datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o keyboard_negative_permdisp_results_3

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df  Sum Sq  Mean Sq F value    Pr(>F)    
	Groups     10 1.06059 0.106059  46.827 < 2.2e-16 ***
	Residuals 104 0.23555 0.002265                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
	Groups     10 1.06059 0.106059 46.827    999  0.001 ***
	Residuals 104 0.23555 0.002265                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		  F1 L1 L3        M1        M2        M3        M9 R1 U1 U2 U3
	F1                 0.8600000 0.0160000 0.0190000 0.0230000            
	L1                                                                    
	L3                                                                    
	M1 0.8341713                 0.0360000 0.0450000 0.0350000            
	M2 0.0029092       0.0085542           0.6310000 0.6030000            
	M3 0.0033429       0.0089829 0.6263920           0.3360000            
	M9 0.0034822       0.0099070 0.5808124 0.3261826                      
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

	compare_categories.py --method permdisp -i datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c BODY_SITE -o  whole_body_positive_permdisp_results

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df Sum Sq  Mean Sq F value    Pr(>F)    
	Groups     19 2.5948 0.136569  42.752 < 2.2e-16 ***
	Residuals 565 1.8048 0.003194                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
	Groups     19 2.5948 0.136569 42.752    999  0.001 ***
	Residuals 565 1.8048 0.003194                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		                         UBERON:ear canal UBERON:feces
	UBERON:ear canal                                    4.4200e-01
	UBERON:feces                           4.6703e-01             
	UBERON:glans penis                     9.2269e-01   5.8869e-01
	UBERON:hair                            2.9577e-01   1.0502e-01
	UBERON:labia minora                    4.6547e-01   1.2309e-02
	UBERON:mouth                           3.0444e-05   8.1281e-22
	UBERON:nose                            9.6453e-01   3.6454e-01
	UBERON:nostril                         9.3498e-01   3.1264e-01
	UBERON:nostrils                        1.9969e-01   6.3745e-06
	UBERON:skin of arm                     5.8680e-01   8.6330e-01
	UBERON:skin of finger                  7.3413e-01   2.0283e-02
	UBERON:skin of forearm                 8.7177e-01   2.5248e-01
	UBERON:tongue                          4.0432e-10   1.2830e-31
	UBERON:urine                           1.6228e-02   1.2239e-04
	UBERON:zone of skin of abdomen         7.6712e-01   7.4643e-01
	UBERON:zone of skin of foot            1.0107e-04   4.8215e-11
	UBERON:zone of skin of hand            3.3643e-01   7.8337e-01
	UBERON:zone of skin of head            3.3284e-01   3.3459e-01
	UBERON:zone of skin of knee            3.3201e-01   5.9994e-01
	UBERON:zone of skin of outer ear       9.4726e-03   3.7176e-05
		                         UBERON:glans penis UBERON:hair
	UBERON:ear canal                         9.2500e-01  2.8600e-01
	UBERON:feces                             5.5400e-01  9.6000e-02
	UBERON:glans penis                                   2.7700e-01
	UBERON:hair                              2.9565e-01            
	UBERON:labia minora                      3.7729e-01  4.1527e-02
	UBERON:mouth                             9.8813e-07  1.6117e-11
	UBERON:nose                              9.2368e-01  1.4849e-01
	UBERON:nostril                           9.3114e-01  1.0019e-01
	UBERON:nostrils                          7.0899e-02  1.0141e-04
	UBERON:skin of arm                       7.1483e-01  4.2161e-01
	UBERON:skin of finger                    5.4488e-01  6.5167e-03
	UBERON:skin of forearm                   9.6913e-01  3.2223e-02
	UBERON:tongue                            4.1893e-10  3.7905e-17
	UBERON:urine                             5.7014e-02  1.5190e-01
	UBERON:zone of skin of abdomen           8.3580e-01  2.0885e-01
	UBERON:zone of skin of foot              6.1790e-04  4.5561e-03
	UBERON:zone of skin of hand              4.7629e-01  1.0440e-01
	UBERON:zone of skin of head              4.4564e-01  6.0786e-01
	UBERON:zone of skin of knee              3.6348e-01  1.0568e-01
	UBERON:zone of skin of outer ear         2.0701e-02  2.9813e-01
		                         UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                          4.9400e-01   1.0000e-03  9.6200e-01
	UBERON:feces                              8.0000e-03   1.0000e-03  3.5900e-01
	UBERON:glans penis                        3.8600e-01   1.0000e-03  9.1800e-01
	UBERON:hair                               4.3000e-02   1.0000e-03  1.4200e-01
	UBERON:labia minora                                    2.0000e-03  2.7600e-01
	UBERON:mouth                              6.4891e-04               1.0000e-03
	UBERON:nose                               2.8123e-01   1.1077e-08            
	UBERON:nostril                            1.9200e-01   5.4186e-11  9.7457e-01
	UBERON:nostrils                           9.0577e-01   3.6008e-11  4.2741e-02
	UBERON:skin of arm                        1.3607e-01   1.0669e-09  5.3749e-01
	UBERON:skin of finger                     2.0815e-01   1.2757e-14  5.6159e-01
	UBERON:skin of forearm                    7.1926e-02   5.7802e-17  8.8394e-01
	UBERON:tongue                             1.6838e-06   1.3896e-02  5.5519e-14
	UBERON:urine                              4.8876e-03   7.1717e-14  6.8408e-03
	UBERON:zone of skin of abdomen            1.4692e-01   2.7055e-11  6.9658e-01
	UBERON:zone of skin of foot               3.0948e-06   1.3613e-27  3.4290e-06
	UBERON:zone of skin of hand               3.8491e-03   4.2418e-27  2.4330e-01
	UBERON:zone of skin of head               4.0950e-02   3.2674e-13  2.4563e-01
	UBERON:zone of skin of knee               2.0380e-03   8.2393e-28  1.7760e-01
	UBERON:zone of skin of outer ear          1.1439e-04   8.2544e-26  1.3170e-03
		                         UBERON:nostril UBERON:nostrils
	UBERON:ear canal                     9.3900e-01      2.0200e-01
	UBERON:feces                         3.1700e-01      1.0000e-03
	UBERON:glans penis                   9.1500e-01      5.4000e-02
	UBERON:hair                          1.0000e-01      1.0000e-03
	UBERON:labia minora                  1.9400e-01      9.1700e-01
	UBERON:mouth                         1.0000e-03      1.0000e-03
	UBERON:nose                          9.7300e-01      4.6000e-02
	UBERON:nostril                                       2.8000e-02
	UBERON:nostrils                      2.2153e-02                
	UBERON:skin of arm                   4.5306e-01      9.3573e-03
	UBERON:skin of finger                4.7183e-01      2.0926e-02
	UBERON:skin of forearm               9.0753e-01      4.7778e-04
	UBERON:tongue                        3.8781e-19      3.0294e-15
	UBERON:urine                         4.2379e-04      2.6698e-06
	UBERON:zone of skin of abdomen       6.9015e-01      2.7637e-03
	UBERON:zone of skin of foot          1.0358e-08      3.0026e-14
	UBERON:zone of skin of hand          1.7919e-01      4.0282e-07
	UBERON:zone of skin of head          1.5769e-01      3.6122e-04
	UBERON:zone of skin of knee          1.5286e-01      5.8629e-09
	UBERON:zone of skin of outer ear     9.9823e-05      7.8190e-11
		                         UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                         5.7900e-01            7.1700e-01
	UBERON:feces                             8.6000e-01            1.3000e-02
	UBERON:glans penis                       6.9400e-01            5.1500e-01
	UBERON:hair                              4.3200e-01            4.0000e-03
	UBERON:labia minora                      1.4200e-01            2.1300e-01
	UBERON:mouth                             1.0000e-03            1.0000e-03
	UBERON:nose                              5.3300e-01            5.5500e-01
	UBERON:nostril                           4.7700e-01            4.7500e-01
	UBERON:nostrils                          1.2000e-02            2.4000e-02
	UBERON:skin of arm                                             1.3000e-01
	UBERON:skin of finger                    1.4376e-01                      
	UBERON:skin of forearm                   4.4993e-01            2.5820e-01
	UBERON:tongue                            1.1189e-17            9.3587e-22
	UBERON:urine                             9.0740e-03            1.8513e-05
	UBERON:zone of skin of abdomen           7.8918e-01            2.1091e-01
	UBERON:zone of skin of foot              1.4563e-05            1.9922e-12
	UBERON:zone of skin of hand              9.7285e-01            5.3697e-03
	UBERON:zone of skin of head              6.3052e-01            1.9369e-02
	UBERON:zone of skin of knee              9.2382e-01            1.6985e-03
	UBERON:zone of skin of outer ear         1.3609e-02            9.5388e-08
		                         UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                             8.7200e-01    1.0000e-03
	UBERON:feces                                 2.6100e-01    1.0000e-03
	UBERON:glans penis                           9.7400e-01    1.0000e-03
	UBERON:hair                                  2.8000e-02    1.0000e-03
	UBERON:labia minora                          6.3000e-02    1.0000e-03
	UBERON:mouth                                 1.0000e-03    1.5000e-02
	UBERON:nose                                  8.9900e-01    1.0000e-03
	UBERON:nostril                               9.1000e-01    1.0000e-03
	UBERON:nostrils                              2.0000e-03    1.0000e-03
	UBERON:skin of arm                           4.6900e-01    1.0000e-03
	UBERON:skin of finger                        2.3500e-01    1.0000e-03
	UBERON:skin of forearm                                     1.0000e-03
	UBERON:tongue                                8.3044e-23              
	UBERON:urine                                 3.3783e-04    1.3054e-25
	UBERON:zone of skin of abdomen               6.2517e-01    3.2151e-15
	UBERON:zone of skin of foot                  5.4854e-10    2.7671e-43
	UBERON:zone of skin of hand                  1.4111e-01    1.4827e-40
	UBERON:zone of skin of head                  1.2625e-01    1.4433e-22
	UBERON:zone of skin of knee                  5.5579e-02    1.7496e-34
	UBERON:zone of skin of outer ear             1.1628e-05    1.1267e-39
		                         UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                   1.6000e-02                     7.6500e-01
	UBERON:feces                       1.0000e-03                     7.5300e-01
	UBERON:glans penis                 6.1000e-02                     8.4900e-01
	UBERON:hair                        1.5500e-01                     2.2800e-01
	UBERON:labia minora                7.0000e-03                     1.4300e-01
	UBERON:mouth                       1.0000e-03                     1.0000e-03
	UBERON:nose                        5.0000e-03                     7.1400e-01
	UBERON:nostril                     1.0000e-03                     6.9900e-01
	UBERON:nostrils                    1.0000e-03                     2.0000e-03
	UBERON:skin of arm                 6.0000e-03                     7.9100e-01
	UBERON:skin of finger              1.0000e-03                     2.1400e-01
	UBERON:skin of forearm             1.0000e-03                     6.3100e-01
	UBERON:tongue                      1.0000e-03                     1.0000e-03
	UBERON:urine                                                      2.7000e-02
	UBERON:zone of skin of abdomen     2.0556e-02                               
	UBERON:zone of skin of foot        5.0033e-01                     2.3277e-05
	UBERON:zone of skin of hand        1.1372e-05                     6.0591e-01
	UBERON:zone of skin of head        1.2425e-02                     4.4204e-01
	UBERON:zone of skin of knee        2.3048e-04                     4.2903e-01
	UBERON:zone of skin of outer ear   1.1364e-01                     5.5141e-03
		                         UBERON:zone of skin of foot
	UBERON:ear canal                                  1.0000e-03
	UBERON:feces                                      1.0000e-03
	UBERON:glans penis                                2.0000e-03
	UBERON:hair                                       6.0000e-03
	UBERON:labia minora                               1.0000e-03
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       1.0000e-03
	UBERON:nostril                                    1.0000e-03
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                1.0000e-03
	UBERON:skin of finger                             1.0000e-03
	UBERON:skin of forearm                            1.0000e-03
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      4.9900e-01
	UBERON:zone of skin of abdomen                    1.0000e-03
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                       1.8751e-13
	UBERON:zone of skin of head                       1.4751e-05
	UBERON:zone of skin of knee                       4.5478e-11
	UBERON:zone of skin of outer ear                  8.4740e-04
		                         UBERON:zone of skin of hand
	UBERON:ear canal                                  2.9500e-01
	UBERON:feces                                      7.9900e-01
	UBERON:glans penis                                4.5500e-01
	UBERON:hair                                       1.1200e-01
	UBERON:labia minora                               5.0000e-03
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       2.4400e-01
	UBERON:nostril                                    1.7600e-01
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                9.6800e-01
	UBERON:skin of finger                             4.0000e-03
	UBERON:skin of forearm                            1.2600e-01
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      1.0000e-03
	UBERON:zone of skin of abdomen                    6.4100e-01
	UBERON:zone of skin of foot                       1.0000e-03
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                       3.5927e-01
	UBERON:zone of skin of knee                       7.9743e-01
	UBERON:zone of skin of outer ear                  7.9520e-06
		                         UBERON:zone of skin of head
	UBERON:ear canal                                  3.2900e-01
	UBERON:feces                                      3.6100e-01
	UBERON:glans penis                                4.3600e-01
	UBERON:hair                                       6.3800e-01
	UBERON:labia minora                               3.9000e-02
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       2.5700e-01
	UBERON:nostril                                    1.7400e-01
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                6.5700e-01
	UBERON:skin of finger                             1.8000e-02
	UBERON:skin of forearm                            1.2800e-01
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      1.0000e-02
	UBERON:zone of skin of abdomen                    4.6200e-01
	UBERON:zone of skin of foot                       1.0000e-03
	UBERON:zone of skin of hand                       3.6500e-01
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                       4.7827e-01
	UBERON:zone of skin of outer ear                  3.3027e-02
		                         UBERON:zone of skin of knee
	UBERON:ear canal                                  3.0600e-01
	UBERON:feces                                      5.9400e-01
	UBERON:glans penis                                3.6600e-01
	UBERON:hair                                       1.0900e-01
	UBERON:labia minora                               5.0000e-03
	UBERON:mouth                                      1.0000e-03
	UBERON:nose                                       1.8700e-01
	UBERON:nostril                                    1.5300e-01
	UBERON:nostrils                                   1.0000e-03
	UBERON:skin of arm                                9.2000e-01
	UBERON:skin of finger                             2.0000e-03
	UBERON:skin of forearm                            4.7000e-02
	UBERON:tongue                                     1.0000e-03
	UBERON:urine                                      1.0000e-03
	UBERON:zone of skin of abdomen                    4.4100e-01
	UBERON:zone of skin of foot                       1.0000e-03
	UBERON:zone of skin of hand                       7.8900e-01
	UBERON:zone of skin of head                       4.8500e-01
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                  5.2485e-05
		                         UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.009
	UBERON:feces                                                0.001
	UBERON:glans penis                                          0.024
	UBERON:hair                                                 0.316
	UBERON:labia minora                                         0.001
	UBERON:mouth                                                0.001
	UBERON:nose                                                 0.001
	UBERON:nostril                                              0.001
	UBERON:nostrils                                             0.001
	UBERON:skin of arm                                          0.013
	UBERON:skin of finger                                       0.001
	UBERON:skin of forearm                                      0.001
	UBERON:tongue                                               0.001
	UBERON:urine                                                0.114
	UBERON:zone of skin of abdomen                              0.009
	UBERON:zone of skin of foot                                 0.001
	UBERON:zone of skin of hand                                 0.001
	UBERON:zone of skin of head                                 0.033
	UBERON:zone of skin of knee                                 0.001
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are significant.

Test 2
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m datasets/whole_body/map.txt -c BODY_SITE -o  whole_body_negative_permdisp_results_1

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df Sum Sq   Mean Sq F value    Pr(>F)    
	Groups     19 0.2600 0.0136828  2.3748 0.0009275 ***
	Residuals 565 3.2554 0.0057617                      
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df Sum Sq   Mean Sq      F N.Perm Pr(>F)    
	Groups     19 0.2600 0.0136828 2.3748    999  0.001 ***
	Residuals 565 3.2554 0.0057617                         
	---
	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		                         UBERON:ear canal UBERON:feces
	UBERON:ear canal                                     0.0580000
	UBERON:feces                            0.0618004             
	UBERON:glans penis                      0.6888362    0.0205773
	UBERON:hair                             0.3423893    0.5054131
	UBERON:labia minora                     0.9581384    0.1463686
	UBERON:mouth                            0.8883601    0.0174021
	UBERON:nose                             0.8259411    0.0404160
	UBERON:nostril                          0.1103653    0.9278009
	UBERON:nostrils                         0.1019801    0.7814930
	UBERON:skin of arm                      0.1191182    0.8698623
	UBERON:skin of finger                   0.2467010    0.3757375
	UBERON:skin of forearm                  0.7187030    0.0274418
	UBERON:tongue                           0.9091165    0.0122678
	UBERON:urine                            0.1112210    0.7961795
	UBERON:zone of skin of abdomen          0.1561569    0.7513293
	UBERON:zone of skin of foot             0.0355112    0.6186005
	UBERON:zone of skin of hand             0.0514663    0.8192744
	UBERON:zone of skin of head             0.7732816    0.0279813
	UBERON:zone of skin of knee             0.0196908    0.3222223
	UBERON:zone of skin of outer ear        0.2592380    0.1402208
		                         UBERON:glans penis UBERON:hair
	UBERON:ear canal                          0.7010000   0.3390000
	UBERON:feces                              0.0200000   0.4960000
	UBERON:glans penis                                    0.1150000
	UBERON:hair                               0.1063323            
	UBERON:labia minora                       0.5996063   0.3957491
	UBERON:mouth                              0.6893747   0.1663916
	UBERON:nose                               0.2620375   0.2276860
	UBERON:nostril                            0.0345755   0.5873161
	UBERON:nostrils                           0.0212736   0.4274816
	UBERON:skin of arm                        0.0319633   0.6175236
	UBERON:skin of finger                     0.0831564   0.9748348
	UBERON:skin of forearm                    0.2960381   0.2858409
	UBERON:tongue                             0.5134089   0.2310853
	UBERON:urine                              0.0533612   0.6681297
	UBERON:zone of skin of abdomen            0.0408069   0.4570704
	UBERON:zone of skin of foot               0.0211119   0.3495896
	UBERON:zone of skin of hand               0.0296219   0.4472614
	UBERON:zone of skin of head               0.4246734   0.3283786
	UBERON:zone of skin of knee               0.0077288   0.1910723
	UBERON:zone of skin of outer ear          0.1049619   0.7551278
		                         UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                           0.9630000    0.8950000   0.8430000
	UBERON:feces                               0.1380000    0.0160000   0.0350000
	UBERON:glans penis                         0.6230000    0.6740000   0.2560000
	UBERON:hair                                0.3990000    0.1670000   0.2390000
	UBERON:labia minora                                     0.8360000   0.8870000
	UBERON:mouth                               0.8381557                0.6350000
	UBERON:nose                                0.8763900    0.6024995            
	UBERON:nostril                             0.1933516    0.0360349   0.0641959
	UBERON:nostrils                            0.1432283    0.0269572   0.0312012
	UBERON:skin of arm                         0.1929068    0.0374467   0.0632978
	UBERON:skin of finger                      0.3537723    0.1106086   0.2059961
	UBERON:skin of forearm                     0.8082748    0.5045231   0.8683008
	UBERON:tongue                              0.9841120    0.7413796   0.8784521
	UBERON:urine                               0.2385674    0.0451830   0.1009613
	UBERON:zone of skin of abdomen             0.1933384    0.0516743   0.0503213
	UBERON:zone of skin of foot                0.1269355    0.0113132   0.0308040
	UBERON:zone of skin of hand                0.1597249    0.0181752   0.0480237
	UBERON:zone of skin of head                0.8714306    0.5976492   0.9414590
	UBERON:zone of skin of knee                0.0679258    0.0040474   0.0086909
	UBERON:zone of skin of outer ear           0.4220288    0.1299103   0.3017238
		                         UBERON:nostril UBERON:nostrils
	UBERON:ear canal                      0.0980000       0.1000000
	UBERON:feces                          0.9180000       0.7780000
	UBERON:glans penis                    0.0360000       0.0270000
	UBERON:hair                           0.5870000       0.4240000
	UBERON:labia minora                   0.1800000       0.1410000
	UBERON:mouth                          0.0350000       0.0320000
	UBERON:nose                           0.0620000       0.0270000
	UBERON:nostril                                        0.7430000
	UBERON:nostrils                       0.7456414                
	UBERON:skin of arm                    0.9475908       0.6966236
	UBERON:skin of finger                 0.4795059       0.3329946
	UBERON:skin of forearm                0.0563056       0.0367184
	UBERON:tongue                         0.0345144       0.0314949
	UBERON:urine                          0.8898789       0.6572355
	UBERON:zone of skin of abdomen        0.7261029       0.9412218
	UBERON:zone of skin of foot           0.6094781       0.9196342
	UBERON:zone of skin of hand           0.7782185       0.9258778
	UBERON:zone of skin of head           0.0654480       0.0565963
	UBERON:zone of skin of knee           0.3430842       0.6198555
	UBERON:zone of skin of outer ear      0.2383897       0.1659350
		                         UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                          0.1030000             0.2570000
	UBERON:feces                              0.8820000             0.3830000
	UBERON:glans penis                        0.0320000             0.0890000
	UBERON:hair                               0.6280000             0.9770000
	UBERON:labia minora                       0.1790000             0.3520000
	UBERON:mouth                              0.0320000             0.1150000
	UBERON:nose                               0.0540000             0.2130000
	UBERON:nostril                            0.9530000             0.4600000
	UBERON:nostrils                           0.6750000             0.3200000
	UBERON:skin of arm                                              0.5210000
	UBERON:skin of finger                     0.5190176                      
	UBERON:skin of forearm                    0.0619593             0.2166864
	UBERON:tongue                             0.0407661             0.1378129
	UBERON:urine                              0.9446466             0.5510415
	UBERON:zone of skin of abdomen            0.6831720             0.3713897
	UBERON:zone of skin of foot               0.5667198             0.2012059
	UBERON:zone of skin of hand               0.7284500             0.2950442
	UBERON:zone of skin of head               0.0760855             0.2265783
	UBERON:zone of skin of knee               0.3111467             0.0886925
	UBERON:zone of skin of outer ear          0.2748402             0.7196108
		                         UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                              0.6960000     0.9020000
	UBERON:feces                                  0.0240000     0.0140000
	UBERON:glans penis                            0.2850000     0.5230000
	UBERON:hair                                   0.2940000     0.2470000
	UBERON:labia minora                           0.8130000     0.9850000
	UBERON:mouth                                  0.5090000     0.7460000
	UBERON:nose                                   0.8670000     0.8810000
	UBERON:nostril                                0.0500000     0.0360000
	UBERON:nostrils                               0.0390000     0.0380000
	UBERON:skin of arm                            0.0580000     0.0370000
	UBERON:skin of finger                         0.2430000     0.1590000
	UBERON:skin of forearm                                      0.7550000
	UBERON:tongue                                 0.7359425              
	UBERON:urine                                  0.0705184     0.0311412
	UBERON:zone of skin of abdomen                0.0635432     0.0611634
	UBERON:zone of skin of foot                   0.0126775     0.0035142
	UBERON:zone of skin of hand                   0.0231759     0.0070529
	UBERON:zone of skin of head                   0.9439340     0.7981608
	UBERON:zone of skin of knee                   0.0036323     0.0014552
	UBERON:zone of skin of outer ear              0.2939924     0.1473604
		                         UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                    0.1080000                      0.1800000
	UBERON:feces                        0.8050000                      0.7500000
	UBERON:glans penis                  0.0590000                      0.0400000
	UBERON:hair                         0.6800000                      0.4640000
	UBERON:labia minora                 0.2280000                      0.1840000
	UBERON:mouth                        0.0380000                      0.0580000
	UBERON:nose                         0.1020000                      0.0510000
	UBERON:nostril                      0.8790000                      0.7210000
	UBERON:nostrils                     0.6770000                      0.9330000
	UBERON:skin of arm                  0.9540000                      0.6710000
	UBERON:skin of finger               0.5620000                      0.3770000
	UBERON:skin of forearm              0.0650000                      0.0620000
	UBERON:tongue                       0.0240000                      0.0740000
	UBERON:urine                                                       0.6460000
	UBERON:zone of skin of abdomen      0.6560159                               
	UBERON:zone of skin of foot         0.4548529                      0.9960404
	UBERON:zone of skin of hand         0.6269023                      0.8766239
	UBERON:zone of skin of head         0.0619372                      0.0960377
	UBERON:zone of skin of knee         0.2377532                      0.7385285
	UBERON:zone of skin of outer ear    0.2573702                      0.2148275
		                         UBERON:zone of skin of foot
	UBERON:ear canal                                   0.0280000
	UBERON:feces                                       0.6140000
	UBERON:glans penis                                 0.0230000
	UBERON:hair                                        0.3780000
	UBERON:labia minora                                0.1180000
	UBERON:mouth                                       0.0130000
	UBERON:nose                                        0.0310000
	UBERON:nostril                                     0.5870000
	UBERON:nostrils                                    0.9180000
	UBERON:skin of arm                                 0.5680000
	UBERON:skin of finger                              0.1940000
	UBERON:skin of forearm                             0.0070000
	UBERON:tongue                                      0.0030000
	UBERON:urine                                       0.4570000
	UBERON:zone of skin of abdomen                     0.9940000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                        0.7765505
	UBERON:zone of skin of head                        0.0088502
	UBERON:zone of skin of knee                        0.6055076
	UBERON:zone of skin of outer ear                   0.0406187
		                         UBERON:zone of skin of hand
	UBERON:ear canal                                   0.0510000
	UBERON:feces                                       0.8290000
	UBERON:glans penis                                 0.0350000
	UBERON:hair                                        0.4630000
	UBERON:labia minora                                0.1440000
	UBERON:mouth                                       0.0220000
	UBERON:nose                                        0.0410000
	UBERON:nostril                                     0.7700000
	UBERON:nostrils                                    0.9010000
	UBERON:skin of arm                                 0.7190000
	UBERON:skin of finger                              0.3110000
	UBERON:skin of forearm                             0.0300000
	UBERON:tongue                                      0.0080000
	UBERON:urine                                       0.6220000
	UBERON:zone of skin of abdomen                     0.8790000
	UBERON:zone of skin of foot                        0.7790000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                        0.0167137
	UBERON:zone of skin of knee                        0.4397192
	UBERON:zone of skin of outer ear                   0.0792463
		                         UBERON:zone of skin of head
	UBERON:ear canal                                   0.7640000
	UBERON:feces                                       0.0270000
	UBERON:glans penis                                 0.4320000
	UBERON:hair                                        0.3540000
	UBERON:labia minora                                0.8540000
	UBERON:mouth                                       0.6080000
	UBERON:nose                                        0.9420000
	UBERON:nostril                                     0.0500000
	UBERON:nostrils                                    0.0580000
	UBERON:skin of arm                                 0.0670000
	UBERON:skin of finger                              0.2260000
	UBERON:skin of forearm                             0.9620000
	UBERON:tongue                                      0.8150000
	UBERON:urine                                       0.0530000
	UBERON:zone of skin of abdomen                     0.1020000
	UBERON:zone of skin of foot                        0.0030000
	UBERON:zone of skin of hand                        0.0230000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                        0.0037957
	UBERON:zone of skin of outer ear                   0.2618914
		                         UBERON:zone of skin of knee
	UBERON:ear canal                                   0.0200000
	UBERON:feces                                       0.3350000
	UBERON:glans penis                                 0.0110000
	UBERON:hair                                        0.2150000
	UBERON:labia minora                                0.0520000
	UBERON:mouth                                       0.0050000
	UBERON:nose                                        0.0100000
	UBERON:nostril                                     0.3370000
	UBERON:nostrils                                    0.6020000
	UBERON:skin of arm                                 0.3170000
	UBERON:skin of finger                              0.0860000
	UBERON:skin of forearm                             0.0030000
	UBERON:tongue                                      0.0020000
	UBERON:urine                                       0.2360000
	UBERON:zone of skin of abdomen                     0.7290000
	UBERON:zone of skin of foot                        0.6150000
	UBERON:zone of skin of hand                        0.4220000
	UBERON:zone of skin of head                        0.0040000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                   0.0143988
		                         UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.240
	UBERON:feces                                                0.126
	UBERON:glans penis                                          0.119
	UBERON:hair                                                 0.764
	UBERON:labia minora                                         0.405
	UBERON:mouth                                                0.118
	UBERON:nose                                                 0.301
	UBERON:nostril                                              0.253
	UBERON:nostrils                                             0.164
	UBERON:skin of arm                                          0.293
	UBERON:skin of finger                                       0.728
	UBERON:skin of forearm                                      0.294
	UBERON:tongue                                               0.154
	UBERON:urine                                                0.268
	UBERON:zone of skin of abdomen                              0.210
	UBERON:zone of skin of foot                                 0.044
	UBERON:zone of skin of hand                                 0.086
	UBERON:zone of skin of head                                 0.266
	UBERON:zone of skin of knee                                 0.010
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are significant.
	
Test 3
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m datasets/whole_body/map.txt -c BODY_SITE -o  whole_body_negative_permdisp_results_2

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df Sum Sq   Mean Sq F value Pr(>F)
	Groups     19 0.1035 0.0054482  0.9055 0.5763
	Residuals 565 3.3996 0.0060169               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     19 0.1035 0.0054482 0.9055    999  0.561
	Residuals 565 3.3996 0.0060169                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		                         UBERON:ear canal UBERON:feces
	UBERON:ear canal                                     0.9470000
	UBERON:feces                            0.9491718             
	UBERON:glans penis                      0.6632282    0.5907393
	UBERON:hair                             0.9588931    0.8902878
	UBERON:labia minora                     0.9527453    0.9760672
	UBERON:mouth                            0.4837764    0.2875361
	UBERON:nose                             0.9395697    0.8640260
	UBERON:nostril                          0.9994857    0.9323520
	UBERON:nostrils                         0.5769986    0.4797631
	UBERON:skin of arm                      0.7699735    0.6117978
	UBERON:skin of finger                   0.9981619    0.9337492
	UBERON:skin of forearm                  0.8136912    0.6666330
	UBERON:tongue                           0.3351867    0.1937985
	UBERON:urine                            0.8600845    0.8567278
	UBERON:zone of skin of abdomen          0.1783487    0.0687944
	UBERON:zone of skin of foot             0.5050372    0.3257261
	UBERON:zone of skin of hand             0.3702147    0.1874797
	UBERON:zone of skin of head             0.8467355    0.8410195
	UBERON:zone of skin of knee             0.5479070    0.3159727
	UBERON:zone of skin of outer ear        0.9734156    0.8800934
		                         UBERON:glans penis UBERON:hair
	UBERON:ear canal                          0.6630000   0.9620000
	UBERON:feces                              0.5950000   0.8790000
	UBERON:glans penis                                    0.6030000
	UBERON:hair                               0.6097327            
	UBERON:labia minora                       0.7598062   0.9158982
	UBERON:mouth                              0.3349791   0.4860120
	UBERON:nose                               0.5968217   0.9791971
	UBERON:nostril                            0.5835292   0.9472192
	UBERON:nostrils                           0.9530852   0.5129808
	UBERON:skin of arm                        0.4406562   0.8058809
	UBERON:skin of finger                     0.6284870   0.9555471
	UBERON:skin of forearm                    0.5088490   0.8480225
	UBERON:tongue                             0.8764168   0.2744082
	UBERON:urine                              0.6748019   0.8006441
	UBERON:zone of skin of abdomen            0.1137767   0.1624004
	UBERON:zone of skin of foot               0.9654872   0.4444102
	UBERON:zone of skin of hand               0.9357535   0.3111455
	UBERON:zone of skin of head               0.6912625   0.7872091
	UBERON:zone of skin of knee               0.2915441   0.5733335
	UBERON:zone of skin of outer ear          0.5663692   0.9745532
		                         UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                           0.9410000    0.4930000   0.9410000
	UBERON:feces                               0.9800000    0.2990000   0.8550000
	UBERON:glans penis                         0.7500000    0.3110000   0.5950000
	UBERON:hair                                0.9200000    0.4830000   0.9810000
	UBERON:labia minora                                     0.5540000   0.9120000
	UBERON:mouth                               0.5498655                0.5010000
	UBERON:nose                                0.9006067    0.5012596            
	UBERON:nostril                             0.9420900    0.3618682   0.9226619
	UBERON:nostrils                            0.7098940    0.1794403   0.4957653
	UBERON:skin of arm                         0.7668182    0.5409303   0.8298140
	UBERON:skin of finger                      0.9480286    0.4129452   0.9337049
	UBERON:skin of forearm                     0.8106983    0.5641927   0.8695295
	UBERON:tongue                              0.5089190    0.0547708   0.2604980
	UBERON:urine                               0.9552521    0.2571717   0.7763594
	UBERON:zone of skin of abdomen             0.2255571    0.5859020   0.1710448
	UBERON:zone of skin of foot                0.6852249    0.1024413   0.4261175
	UBERON:zone of skin of hand                0.5692975    0.0523750   0.2953040
	UBERON:zone of skin of head                0.9421462    0.2671659   0.7636669
	UBERON:zone of skin of knee                0.5993041    0.6817933   0.5967322
	UBERON:zone of skin of outer ear           0.9242949    0.3524776   0.9493880
		                         UBERON:nostril UBERON:nostrils
	UBERON:ear canal                      0.9990000       0.5770000
	UBERON:feces                          0.9250000       0.4500000
	UBERON:glans penis                    0.5920000       0.9480000
	UBERON:hair                           0.9420000       0.5090000
	UBERON:labia minora                   0.9560000       0.7170000
	UBERON:mouth                          0.3730000       0.1960000
	UBERON:nose                           0.9230000       0.4840000
	UBERON:nostril                                        0.4800000
	UBERON:nostrils                       0.4742374                
	UBERON:skin of arm                    0.7010992       0.2979962
	UBERON:skin of finger                 0.9969102       0.5193891
	UBERON:skin of forearm                0.7524657       0.3654254
	UBERON:tongue                         0.2086617       0.7397134
	UBERON:urine                          0.8107780       0.5891034
	UBERON:zone of skin of abdomen        0.0973733       0.0352175
	UBERON:zone of skin of foot           0.3548891       0.9746778
	UBERON:zone of skin of hand           0.2185954       0.8172248
	UBERON:zone of skin of head           0.7975624       0.6186635
	UBERON:zone of skin of knee           0.4207663       0.1406792
	UBERON:zone of skin of outer ear      0.9624945       0.4332690
		                         UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                          0.7610000             0.9980000
	UBERON:feces                              0.5960000             0.9420000
	UBERON:glans penis                        0.4330000             0.6480000
	UBERON:hair                               0.8150000             0.9640000
	UBERON:labia minora                       0.7830000             0.9550000
	UBERON:mouth                              0.5510000             0.4360000
	UBERON:nose                               0.8290000             0.9410000
	UBERON:nostril                            0.7100000             0.9980000
	UBERON:nostrils                           0.2990000             0.5400000
	UBERON:skin of arm                                              0.7170000
	UBERON:skin of finger                     0.7301521                      
	UBERON:skin of forearm                    0.9731041             0.7745033
	UBERON:tongue                             0.1000571             0.2511161
	UBERON:urine                              0.5211068             0.8192846
	UBERON:zone of skin of abdomen            0.1785814             0.1477113
	UBERON:zone of skin of foot               0.1897576             0.3761432
	UBERON:zone of skin of hand               0.1006000             0.2433273
	UBERON:zone of skin of head               0.5221982             0.8105041
	UBERON:zone of skin of knee               0.7161439             0.4593725
	UBERON:zone of skin of outer ear          0.7136199             0.9679383
		                         UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                              0.8010000     0.3150000
	UBERON:feces                                  0.6510000     0.1780000
	UBERON:glans penis                            0.4960000     0.8820000
	UBERON:hair                                   0.8460000     0.2940000
	UBERON:labia minora                           0.8180000     0.5120000
	UBERON:mouth                                  0.5790000     0.0590000
	UBERON:nose                                   0.8610000     0.2570000
	UBERON:nostril                                0.7400000     0.2090000
	UBERON:nostrils                               0.3840000     0.7490000
	UBERON:skin of arm                            0.9700000     0.1060000
	UBERON:skin of finger                         0.7790000     0.2640000
	UBERON:skin of forearm                                      0.1470000
	UBERON:tongue                                 0.1441337              
	UBERON:urine                                  0.5750651     0.2792883
	UBERON:zone of skin of abdomen                0.2272929     0.0055520
	UBERON:zone of skin of foot                   0.2300373     0.7154321
	UBERON:zone of skin of hand                   0.1334654     0.8852895
	UBERON:zone of skin of head                   0.5806936     0.3173780
	UBERON:zone of skin of knee                   0.7086137     0.0249034
	UBERON:zone of skin of outer ear              0.7587793     0.1608156
		                         UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                    0.8680000                      0.1800000
	UBERON:feces                        0.8420000                      0.0720000
	UBERON:glans penis                  0.6780000                      0.0900000
	UBERON:hair                         0.8240000                      0.1600000
	UBERON:labia minora                 0.9710000                      0.2240000
	UBERON:mouth                        0.2520000                      0.5980000
	UBERON:nose                         0.7840000                      0.1660000
	UBERON:nostril                      0.8090000                      0.0860000
	UBERON:nostrils                     0.5720000                      0.0300000
	UBERON:skin of arm                  0.4880000                      0.1790000
	UBERON:skin of finger               0.8490000                      0.1430000
	UBERON:skin of forearm              0.5570000                      0.2220000
	UBERON:tongue                       0.2770000                      0.0050000
	UBERON:urine                                                       0.0610000
	UBERON:zone of skin of abdomen      0.0675393                               
	UBERON:zone of skin of foot         0.4293941                      0.0230360
	UBERON:zone of skin of hand         0.2694508                      0.0078468
	UBERON:zone of skin of head         0.9746332                      0.0662097
	UBERON:zone of skin of knee         0.2488212                      0.2408571
	UBERON:zone of skin of outer ear    0.7394463                      0.1096576
		                         UBERON:zone of skin of foot
	UBERON:ear canal                                   0.5160000
	UBERON:feces                                       0.3060000
	UBERON:glans penis                                 0.9640000
	UBERON:hair                                        0.4680000
	UBERON:labia minora                                0.7130000
	UBERON:mouth                                       0.1010000
	UBERON:nose                                        0.4200000
	UBERON:nostril                                     0.3590000
	UBERON:nostrils                                    0.9660000
	UBERON:skin of arm                                 0.1970000
	UBERON:skin of finger                              0.3930000
	UBERON:skin of forearm                             0.2230000
	UBERON:tongue                                      0.7230000
	UBERON:urine                                       0.4390000
	UBERON:zone of skin of abdomen                     0.0260000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                        0.7788859
	UBERON:zone of skin of head                        0.4980717
	UBERON:zone of skin of knee                        0.0489356
	UBERON:zone of skin of outer ear                   0.2362975
		                         UBERON:zone of skin of hand
	UBERON:ear canal                                   0.3810000
	UBERON:feces                                       0.1840000
	UBERON:glans penis                                 0.9340000
	UBERON:hair                                        0.3370000
	UBERON:labia minora                                0.5810000
	UBERON:mouth                                       0.0590000
	UBERON:nose                                        0.2730000
	UBERON:nostril                                     0.2320000
	UBERON:nostrils                                    0.8160000
	UBERON:skin of arm                                 0.0940000
	UBERON:skin of finger                              0.2500000
	UBERON:skin of forearm                             0.1280000
	UBERON:tongue                                      0.8880000
	UBERON:urine                                       0.2830000
	UBERON:zone of skin of abdomen                     0.0060000
	UBERON:zone of skin of foot                        0.8070000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                        0.3327359
	UBERON:zone of skin of knee                        0.0178476
	UBERON:zone of skin of outer ear                   0.1276731
		                         UBERON:zone of skin of head
	UBERON:ear canal                                   0.8530000
	UBERON:feces                                       0.8300000
	UBERON:glans penis                                 0.6910000
	UBERON:hair                                        0.7770000
	UBERON:labia minora                                0.9490000
	UBERON:mouth                                       0.2730000
	UBERON:nose                                        0.7590000
	UBERON:nostril                                     0.8020000
	UBERON:nostrils                                    0.6120000
	UBERON:skin of arm                                 0.5180000
	UBERON:skin of finger                              0.8340000
	UBERON:skin of forearm                             0.5790000
	UBERON:tongue                                      0.3060000
	UBERON:urine                                       0.9770000
	UBERON:zone of skin of abdomen                     0.0690000
	UBERON:zone of skin of foot                        0.4950000
	UBERON:zone of skin of hand                        0.3430000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                        0.2672596
	UBERON:zone of skin of outer ear                   0.7391490
		                         UBERON:zone of skin of knee
	UBERON:ear canal                                   0.5420000
	UBERON:feces                                       0.2860000
	UBERON:glans penis                                 0.2920000
	UBERON:hair                                        0.6050000
	UBERON:labia minora                                0.6130000
	UBERON:mouth                                       0.6710000
	UBERON:nose                                        0.6000000
	UBERON:nostril                                     0.4430000
	UBERON:nostrils                                    0.1400000
	UBERON:skin of arm                                 0.7130000
	UBERON:skin of finger                              0.4690000
	UBERON:skin of forearm                             0.7190000
	UBERON:tongue                                      0.0260000
	UBERON:urine                                       0.2410000
	UBERON:zone of skin of abdomen                     0.2500000
	UBERON:zone of skin of foot                        0.0440000
	UBERON:zone of skin of hand                        0.0140000
	UBERON:zone of skin of head                        0.2720000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                   0.3955359
		                         UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.973
	UBERON:feces                                                0.889
	UBERON:glans penis                                          0.543
	UBERON:hair                                                 0.971
	UBERON:labia minora                                         0.928
	UBERON:mouth                                                0.361
	UBERON:nose                                                 0.949
	UBERON:nostril                                              0.971
	UBERON:nostrils                                             0.456
	UBERON:skin of arm                                          0.699
	UBERON:skin of finger                                       0.968
	UBERON:skin of forearm                                      0.768
	UBERON:tongue                                               0.148
	UBERON:urine                                                0.742
	UBERON:zone of skin of abdomen                              0.129
	UBERON:zone of skin of foot                                 0.252
	UBERON:zone of skin of hand                                 0.136
	UBERON:zone of skin of head                                 0.742
	UBERON:zone of skin of knee                                 0.392
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are insignificant.
	
Test 4
~~~~~~

**Description:**

This test uses a shuffled unweighted unifrac distance matrix and the BODY_SITE category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m datasets/whole_body/map.txt -c BODY_SITE -o whole_body_negative_permdisp_results_3

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df Sum Sq   Mean Sq F value Pr(>F)
	Groups     19 0.1508 0.0079357  1.3289 0.1582
	Residuals 565 3.3740 0.0059717               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups     19 0.1508 0.0079357 1.3289    999  0.151
	Residuals 565 3.3740 0.0059717                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		                         UBERON:ear canal UBERON:feces
	UBERON:ear canal                                    0.65600000
	UBERON:feces                           0.64369434             
	UBERON:glans penis                     0.13390939   0.22202619
	UBERON:hair                            0.76930293   0.94856194
	UBERON:labia minora                    0.26573188   0.30803348
	UBERON:mouth                           0.38825286   0.13410962
	UBERON:nose                            0.04644476   0.09139076
	UBERON:nostril                         0.77144109   0.31164904
	UBERON:nostrils                        0.80914684   0.40998186
	UBERON:skin of arm                     0.81005310   0.81753053
	UBERON:skin of finger                  0.97422033   0.49272998
	UBERON:skin of forearm                 0.90995260   0.47479702
	UBERON:tongue                          0.40041373   0.08394161
	UBERON:urine                           0.86874942   0.33972854
	UBERON:zone of skin of abdomen         0.52830558   0.79547296
	UBERON:zone of skin of foot            0.75818075   0.22111901
	UBERON:zone of skin of hand            0.93036257   0.51231458
	UBERON:zone of skin of head            0.42023494   0.08561553
	UBERON:zone of skin of knee            0.47005280   0.08974326
	UBERON:zone of skin of outer ear       0.52914309   0.09432887
		                         UBERON:glans penis UBERON:hair
	UBERON:ear canal                         0.13500000  0.75900000
	UBERON:feces                             0.23400000  0.94800000
	UBERON:glans penis                                   0.31100000
	UBERON:hair                              0.31485022            
	UBERON:labia minora                      0.97642427  0.43157394
	UBERON:mouth                             0.02118084  0.28666784
	UBERON:nose                              0.99663888  0.16913109
	UBERON:nostril                           0.03273600  0.50578420
	UBERON:nostrils                          0.08688087  0.58948395
	UBERON:skin of arm                       0.21829905  0.91868905
	UBERON:skin of finger                    0.05666118  0.66846192
	UBERON:skin of forearm                   0.12901681  0.66289356
	UBERON:tongue                            0.00446918  0.23059933
	UBERON:urine                             0.06026260  0.56537298
	UBERON:zone of skin of abdomen           0.27585622  0.80429340
	UBERON:zone of skin of foot              0.05339055  0.46839008
	UBERON:zone of skin of hand              0.08584173  0.72201918
	UBERON:zone of skin of head              0.01838243  0.25437056
	UBERON:zone of skin of knee              0.03276043  0.28277540
	UBERON:zone of skin of outer ear         0.03066872  0.30356640
		                         UBERON:labia minora UBERON:mouth UBERON:nose
	UBERON:ear canal                          0.28700000   0.38800000  0.04700000
	UBERON:feces                              0.32200000   0.14600000  0.09800000
	UBERON:glans penis                        0.97500000   0.01100000  0.99700000
	UBERON:hair                               0.45100000   0.27300000  0.16800000
	UBERON:labia minora                                    0.08300000  0.97300000
	UBERON:mouth                              0.08234783               0.00200000
	UBERON:nose                               0.96595058   0.00258717            
	UBERON:nostril                            0.09288672   0.39987815  0.00450330
	UBERON:nostrils                           0.17972658   0.50771736  0.02158724
	UBERON:skin of arm                        0.31353155   0.25072265  0.09085265
	UBERON:skin of finger                     0.13448037   0.27172046  0.01087135
	UBERON:skin of forearm                    0.21252941   0.43154476  0.03805394
	UBERON:tongue                             0.02619457   0.68689433  0.00020157
	UBERON:urine                              0.11553034   0.34021622  0.01016977
	UBERON:zone of skin of abdomen            0.45492294   0.12130731  0.15474301
	UBERON:zone of skin of foot               0.09691715   0.42433681  0.00767909
	UBERON:zone of skin of hand               0.14408642   0.20852617  0.01769832
	UBERON:zone of skin of head               0.05282882   0.81181553  0.00154299
	UBERON:zone of skin of knee               0.07009782   0.78592682  0.00356419
	UBERON:zone of skin of outer ear          0.06343217   0.65257392  0.00298993
		                         UBERON:nostril UBERON:nostrils
	UBERON:ear canal                     0.76900000      0.81600000
	UBERON:feces                         0.31300000      0.41700000
	UBERON:glans penis                   0.02800000      0.07600000
	UBERON:hair                          0.51800000      0.62000000
	UBERON:labia minora                  0.08800000      0.18600000
	UBERON:mouth                         0.39100000      0.52700000
	UBERON:nose                          0.00700000      0.03100000
	UBERON:nostril                                       0.98500000
	UBERON:nostrils                      0.98945448                
	UBERON:skin of arm                   0.51493588      0.59483323
	UBERON:skin of finger                0.73013131      0.77215159
	UBERON:skin of forearm               0.87122857      0.88834381
	UBERON:tongue                        0.49012080      0.58969183
	UBERON:urine                         0.86525105      0.88109265
	UBERON:zone of skin of abdomen       0.26363915      0.36662477
	UBERON:zone of skin of foot          0.97499842      0.99082569
	UBERON:zone of skin of hand          0.58724157      0.65288167
	UBERON:zone of skin of head          0.46527782      0.56683955
	UBERON:zone of skin of knee          0.51833317      0.61318240
	UBERON:zone of skin of outer ear     0.62419076      0.70425342
		                         UBERON:skin of arm UBERON:skin of finger
	UBERON:ear canal                         0.81500000            0.97300000
	UBERON:feces                             0.79000000            0.50100000
	UBERON:glans penis                       0.21100000            0.05000000
	UBERON:hair                              0.91700000            0.67100000
	UBERON:labia minora                      0.30800000            0.12800000
	UBERON:mouth                             0.23500000            0.27800000
	UBERON:nose                              0.08800000            0.01400000
	UBERON:nostril                           0.50200000            0.73300000
	UBERON:nostrils                          0.58900000            0.76200000
	UBERON:skin of arm                                             0.71200000
	UBERON:skin of finger                    0.71343830                      
	UBERON:skin of forearm                   0.67885778            0.90435132
	UBERON:tongue                            0.20571054            0.28715526
	UBERON:urine                             0.57677297            0.85368451
	UBERON:zone of skin of abdomen           0.69147681            0.39604507
	UBERON:zone of skin of foot              0.44806320            0.69776866
	UBERON:zone of skin of hand              0.78039956            0.86537734
	UBERON:zone of skin of head              0.21094626            0.29103360
	UBERON:zone of skin of knee              0.22911370            0.33493676
	UBERON:zone of skin of outer ear         0.25239010            0.40169838
		                         UBERON:skin of forearm UBERON:tongue
	UBERON:ear canal                             0.91000000    0.37900000
	UBERON:feces                                 0.46500000    0.09000000
	UBERON:glans penis                           0.12900000    0.00500000
	UBERON:hair                                  0.67400000    0.23200000
	UBERON:labia minora                          0.20600000    0.03100000
	UBERON:mouth                                 0.43300000    0.71100000
	UBERON:nose                                  0.02900000    0.00100000
	UBERON:nostril                               0.87900000    0.49500000
	UBERON:nostrils                              0.87500000    0.58000000
	UBERON:skin of arm                           0.68200000    0.21700000
	UBERON:skin of finger                        0.90500000    0.27500000
	UBERON:skin of forearm                                     0.46700000
	UBERON:tongue                                0.46292243              
	UBERON:urine                                 0.97698710    0.37927308
	UBERON:zone of skin of abdomen               0.45176747    0.07543195
	UBERON:zone of skin of foot                  0.83538476    0.50745925
	UBERON:zone of skin of hand                  0.77366898    0.18917276
	UBERON:zone of skin of head                  0.44195087    0.88440190
	UBERON:zone of skin of knee                  0.47432668    0.94269494
	UBERON:zone of skin of outer ear             0.54562590    0.87752123
		                         UBERON:urine UBERON:zone of skin of abdomen
	UBERON:ear canal                   0.88400000                     0.53900000
	UBERON:feces                       0.33500000                     0.79500000
	UBERON:glans penis                 0.04700000                     0.26200000
	UBERON:hair                        0.55300000                     0.80600000
	UBERON:labia minora                0.10200000                     0.47500000
	UBERON:mouth                       0.33900000                     0.12600000
	UBERON:nose                        0.00900000                     0.17000000
	UBERON:nostril                     0.88000000                     0.27200000
	UBERON:nostrils                    0.88000000                     0.35900000
	UBERON:skin of arm                 0.57700000                     0.69600000
	UBERON:skin of finger              0.85100000                     0.40400000
	UBERON:skin of forearm             0.97200000                     0.45300000
	UBERON:tongue                      0.38200000                     0.07400000
	UBERON:urine                                                      0.33400000
	UBERON:zone of skin of abdomen     0.33966368                               
	UBERON:zone of skin of foot        0.81735802                     0.27999512
	UBERON:zone of skin of hand        0.68035246                     0.46561810
	UBERON:zone of skin of head        0.34651958                     0.11675393
	UBERON:zone of skin of knee        0.37557566                     0.15319375
	UBERON:zone of skin of outer ear   0.45661293                     0.16561517
		                         UBERON:zone of skin of foot
	UBERON:ear canal                                  0.74900000
	UBERON:feces                                      0.22600000
	UBERON:glans penis                                0.05100000
	UBERON:hair                                       0.47800000
	UBERON:labia minora                               0.09900000
	UBERON:mouth                                      0.44500000
	UBERON:nose                                       0.00900000
	UBERON:nostril                                    0.97100000
	UBERON:nostrils                                   0.98900000
	UBERON:skin of arm                                0.42400000
	UBERON:skin of finger                             0.71400000
	UBERON:skin of forearm                            0.83500000
	UBERON:tongue                                     0.51700000
	UBERON:urine                                      0.80800000
	UBERON:zone of skin of abdomen                    0.29200000
	UBERON:zone of skin of foot                                 
	UBERON:zone of skin of hand                       0.48840738
	UBERON:zone of skin of head                       0.44572926
	UBERON:zone of skin of knee                       0.47135924
	UBERON:zone of skin of outer ear                  0.57948514
		                         UBERON:zone of skin of hand
	UBERON:ear canal                                  0.92600000
	UBERON:feces                                      0.50800000
	UBERON:glans penis                                0.08300000
	UBERON:hair                                       0.71400000
	UBERON:labia minora                               0.14300000
	UBERON:mouth                                      0.22500000
	UBERON:nose                                       0.02100000
	UBERON:nostril                                    0.58700000
	UBERON:nostrils                                   0.65900000
	UBERON:skin of arm                                0.77600000
	UBERON:skin of finger                             0.86300000
	UBERON:skin of forearm                            0.76200000
	UBERON:tongue                                     0.18200000
	UBERON:urine                                      0.69500000
	UBERON:zone of skin of abdomen                    0.47600000
	UBERON:zone of skin of foot                       0.50600000
	UBERON:zone of skin of hand                                 
	UBERON:zone of skin of head                       0.17117513
	UBERON:zone of skin of knee                       0.18027667
	UBERON:zone of skin of outer ear                  0.21450973
		                         UBERON:zone of skin of head
	UBERON:ear canal                                  0.42200000
	UBERON:feces                                      0.09400000
	UBERON:glans penis                                0.01300000
	UBERON:hair                                       0.26100000
	UBERON:labia minora                               0.04700000
	UBERON:mouth                                      0.82000000
	UBERON:nose                                       0.00500000
	UBERON:nostril                                    0.46100000
	UBERON:nostrils                                   0.58400000
	UBERON:skin of arm                                0.18600000
	UBERON:skin of finger                             0.29300000
	UBERON:skin of forearm                            0.43800000
	UBERON:tongue                                     0.87400000
	UBERON:urine                                      0.31700000
	UBERON:zone of skin of abdomen                    0.11700000
	UBERON:zone of skin of foot                       0.41500000
	UBERON:zone of skin of hand                       0.15200000
	UBERON:zone of skin of head                                 
	UBERON:zone of skin of knee                       0.95046038
	UBERON:zone of skin of outer ear                  0.77640842
		                         UBERON:zone of skin of knee
	UBERON:ear canal                                  0.46700000
	UBERON:feces                                      0.10100000
	UBERON:glans penis                                0.01800000
	UBERON:hair                                       0.27900000
	UBERON:labia minora                               0.06200000
	UBERON:mouth                                      0.80000000
	UBERON:nose                                       0.01000000
	UBERON:nostril                                    0.53200000
	UBERON:nostrils                                   0.62100000
	UBERON:skin of arm                                0.24900000
	UBERON:skin of finger                             0.35600000
	UBERON:skin of forearm                            0.49500000
	UBERON:tongue                                     0.95000000
	UBERON:urine                                      0.38300000
	UBERON:zone of skin of abdomen                    0.15100000
	UBERON:zone of skin of foot                       0.49200000
	UBERON:zone of skin of hand                       0.19300000
	UBERON:zone of skin of head                       0.94500000
	UBERON:zone of skin of knee                                 
	UBERON:zone of skin of outer ear                  0.82348257
		                         UBERON:zone of skin of outer ear
	UBERON:ear canal                                            0.536
	UBERON:feces                                                0.092
	UBERON:glans penis                                          0.025
	UBERON:hair                                                 0.289
	UBERON:labia minora                                         0.059
	UBERON:mouth                                                0.652
	UBERON:nose                                                 0.005
	UBERON:nostril                                              0.643
	UBERON:nostrils                                             0.701
	UBERON:skin of arm                                          0.232
	UBERON:skin of finger                                       0.415
	UBERON:skin of forearm                                      0.534
	UBERON:tongue                                               0.894
	UBERON:urine                                                0.433
	UBERON:zone of skin of abdomen                              0.168
	UBERON:zone of skin of foot                                 0.553
	UBERON:zone of skin of hand                                 0.224
	UBERON:zone of skin of head                                 0.767
	UBERON:zone of skin of knee                                 0.834
	UBERON:zone of skin of outer ear                                 

The p-value indicates that the results are insignificant.

Test 5
~~~~~~

**Description:**

This test uses the unweighted unifrac distance matrix and the SEX category to perform a negative control test.

**Command:** ::

	compare_categories.py --method permdisp -i datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c SEX -o  whole_body_negative_permdisp_results

**Results:**

The following results were written to the output file: ::

	Analysis of Variance Table

	Response: Distances
		   Df Sum Sq   Mean Sq F value Pr(>F)
	Groups      1 0.0119 0.0118769  2.0989 0.1479
	Residuals 583 3.2990 0.0056587               

	Permutation test for homogeneity of multivariate dispersions

	No. of permutations: 999  

	**** STRATA ****
	Permutations are unstratified

	**** SAMPLES ****
	Permutation type: free 
	Mirrored permutations for Samples?: No 

	Response: Distances
		   Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
	Groups      1 0.0119 0.0118769 2.0989    999  0.145
	Residuals 583 3.2990 0.0056587                     

	Pairwise comparisons:
	(Observed p-value below diagonal, permuted p-value above diagonal)
		female  male
	female         0.153
	male   0.14794      

The p-value indicates that the results are insignificant.

References
----------

[1]
http://www.stat.auckland.ac.nz/~mja/Programs.htm

[2]
http://www.stat.auckland.ac.nz/~mja/prog/PERMDISP_UserNotes.pdf

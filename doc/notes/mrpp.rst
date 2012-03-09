==============================================
Multiple Response Permutation Procedure (MRPP)
==============================================

Synopsis
--------
MRPP is a method that tests whether two or more categories are significantly
different. You can specify a category in the metadata mapping file to separate
samples into groups and then test whether there are significant differences
between those groups. For example, you might test whether `Control` samples are
significantly different from `Fast` samples. Since MRPP is non-parametric,
significance is determined through permutations.

Introduction
------------
The Multiple Response Permutation Procedure (MRPP) tests whether there is a
significant difference between two or more groups of samples. It is similar to
ANOSIM except that it uses the original distances instead of ranks,
so it may be more sensitive to outliers. It is also conceptually similar to
ANOVA in that it compares distances within and between groups. MRPP is
non-parametric. The original MRPP paper can be referenced here
[:ref:`2 <mrppref2>`].

MRPP calculates a delta statistic, which is the overall weighted mean of
within-group means of the pairwise distances among samples. The samples and
their pairwise distances are then permuted and delta is calculated for each
permutation of the data. The significance test is the fraction of permuted
deltas that are less than the observed delta. If two groups of samples are
really different, the average of the within-group distances should be less than
the average of the distances between two random collection of samples drawn from
the entire population [:ref:`1 <mrppref1>`].

MRPP also calculates an A statistic, which is:

A = 1 - (observed delta / expected delta)

Thus, an observed A value that is close to zero seems to indicate that there
isn't evidence of clustering because the observed value of our delta statistic
was very close to the expected delta value (which is obtained from
permutations).

MRPP tests whether there are differences between two or more groups of samples,
but it may not always find differences in groups due to differences in means.
Instead, MRPP may find differences in groups based on spread (differences in
within-group distance). For example, it might find that two groups are
significantly different because one group has greater within-group distances
than the other. Thus, it has been recommended to use adonis instead of MRPP when
possible because adonis doesn't have this problem [:ref:`1 <mrppref1>`].

A helpful introduction to MRPP can be found here [:ref:`4 <mrppref4>`].

Existing Implementations
------------------------
There are several existing implementations of MRPP in the following statistical
packages:

* vegan package for R

* MRPP macro for SPSS package [:ref:`3 <mrppref3>`]

* FORTRAN implementation from original author

* possibly others...

As the vegan package has an MRPP implementation and it is free and open-source,
this implementation will be tested. A simple R script wrapping vegan's mrpp
function has been checked into the Qiimeutils repository under
:file:`microbiogeo/r/examples/`. The following sections of the document will
explain how to run the R script.

System Setup and Required Dependencies
--------------------------------------
:note: The following instructions have been tested on 64-bit Linux Mint (essentially Debian). However, they `should` work across different Linux distros and on Macs, though some commands may need to be tweaked, or different package names might have to be used. The instructions assume you use bash as your shell.

The first step is to install R. The following command downloaded and installed R
(for me, it was R version 2.13.1): ::

    sudo apt-get install r-base

Next, you must install the vegan and optparse packages in R. Run the following
commands: ::

    sudo R
    install.packages("vegan")
    install.packages("optparse")
    q()

The install process for the packages will prompt you to choose a mirror to
download them from. Other than that, it is completely automated. On my system, I
ended up with vegan version 2.0-2 and optparse version 0.9.4.

Next, your system must have a version of QIIME installed (I used the latest
version of QIIME in SVN). The MRPP script uses some R utility functions in QIIME
to load data.

Next, you must define an environment variable to tell the MRPP script where to
look for the R utility functions in QIIME. Run the following command, changing
the path to point to the location of your QIIME install: ::

    export QIIME_DIR=/home/jrideout/qiime/trunk

If you don't want to have to perform this step each time you open a new
terminal, run the following command to add it to your .bashrc: ::

    echo "export QIIME_DIR=/home/jrideout/qiime/trunk" >> ~/.bashrc
    source ~/.bashrc

Next, run the following command to test if you can run the MRPP script: ::

    R --slave --args -h < r/examples/mrpp.r

This should run the script in "help" mode. If instructions for how to run the
script are printed, you have successfully configured your system.

Input Files
-----------
The MRPP script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix from the QIIME overview tutorial. You can get the distance
matrix :download:`here <../downloads/overview_unweighted_unifrac_dm.txt>` and
the mapping file :download:`here <../downloads/Fasting_Map.txt>`.

Next, run the following command to execute the MRPP script: ::

    R --slave --args -d overview_unweighted_unifrac_dm.txt -m Fasting_Map.txt -c Treatment < r/examples/mrpp.r

The -c option specifies which column in the mapping file will be used to group
the samples. The `Treatment` column has two values: `Control` and `Fast`. Thus,
MRPP will be used to calculate the dissimilarity between the control and fast
groups.

Output Files
------------
The command in the previous section creates a single output file in the current
directory named :file:`mrpp_results.txt`. The resulting file should look like
this: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 
    
    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          Control Fast  
    delta 0.6237  0.6243
    n     5       4     

    Chance corrected within-group agreement A: 0.07164 
    Based on observed delta 0.624 and expected delta 0.6721 

    Significance of delta: 0.008 
    Based on  999  permutations

The second from the last line contains the p-value of the observed delta
statistic, which is 0.008. This indicates that the differences between `Control`
and `Fast` sample groups is significant, based on 999 permutations.

Testing Results
---------------
This section will describe different tests that were run on the MRPP script.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**

This test uses the `BODY_SITE` category as a positive control. We expect there
to be significant clustering due to previous analysis done on the Whole Body
dataset.

**Command:** ::

    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c BODY_SITE < r/examples/mrpp.r

**Results:**

The following output file is created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 
    
    Dissimilarity index: 
    Weights for groups:  n 
    
    Class means and counts:
    
          UBERON:ear canal UBERON:feces UBERON:glans penis UBERON:hair
    delta 0.6182           0.6209       0.6405             0.6716     
    n     13               43           7                  14         
          UBERON:labia minora UBERON:mouth UBERON:nose UBERON:nostril
    delta 0.5899              0.3782       0.6197      0.6081        
    n     6                   14           14          28            
          UBERON:nostrils UBERON:skin of arm UBERON:skin of finger
    delta 0.5549           0.63              0.5937               
    n     18              26                 28                   
          UBERON:skin of forearm UBERON:tongue UBERON:urine
    delta 0.613                  0.3132        0.7013      
    n     25                     32            46          
          UBERON:zone of skin of abdomen UBERON:zone of skin of foot
    delta 0.6365                         0.7141                     
    n     12                             64                         
          UBERON:zone of skin of hand UBERON:zone of skin of head
    delta 0.6237                      0.6426                     
    n     64                          32                         
          UBERON:zone of skin of knee UBERON:zone of skin of outer ear
    delta 0.6286                      0.6663                          
    n     41                          58                              
    
    Chance corrected within-group agreement A: 0.1524 
    Based on observed delta 0.6188 and expected delta 0.7301 
    
    Significance of delta: 0.001 
    Based on  999  permutations

The p-value of 0.001 indicates that body sites are significantly different (i.e.
there is clustering). This is a result that we would expect.

Test 2
~~~~~~
**Description:**

This test uses the `SEX` category as a negative control. We don't expect to see
significant clustering due to previous analysis done on the Whole Body dataset.

**Command:** ::

    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm.txt -m datasets/whole_body/map.txt -c SEX < r/examples/mrpp.r

**Results:**

The following output file is created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          female male  
    delta 0.7364 0.7221
    n     234    351   

    Chance corrected within-group agreement A: 0.003149 
    Based on observed delta 0.7278 and expected delta 0.7301 

    Significance of delta: 0.001 
    Based on  999  permutations

The p-value of 0.001 indicates that there is significant clustering based on sex
of the subjects. This result isn't something that we'd expect to see. The A
statistic (chance corrected within-group agreement) is pretty close to zero,
though, so this indicates that there might not be clustering (the p-value
doesn't back up this claim, though). This seems to be an issue with ANOSIM as
well, where the p-value claims significance but the test statistic says
otherwise.

Test 3
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `BODY_SITE` category to
perform three negative control tests. Since the labels of the distance matrices
are shuffled, we don't expect to see clustering any more on this category.

**Command:** ::

    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm_shuffled_1.txt -m datasets/whole_body/map.txt -c BODY_SITE < r/examples/mrpp.r
    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm_shuffled_2.txt -m datasets/whole_body/map.txt -c BODY_SITE < r/examples/mrpp.r
    R --slave --args -d datasets/whole_body/unweighted_unifrac_dm_shuffled_3.txt -m datasets/whole_body/map.txt -c BODY_SITE < r/examples/mrpp.r

**Results:**

The following output files are created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          UBERON:ear canal UBERON:feces UBERON:glans penis UBERON:hair
    delta 0.6971           0.7409       0.6915              0.74      
    n     13               43           7                  14         
          UBERON:labia minora UBERON:mouth UBERON:nose UBERON:nostril
    delta 0.7333              0.6822       0.6976      0.7427        
    n     6                   14           14          28            
          UBERON:nostrils UBERON:skin of arm UBERON:skin of finger
    delta 0.762           0.7419             0.7241               
    n     18              26                 28                   
          UBERON:skin of forearm UBERON:tongue UBERON:urine
    delta 0.6865                 0.6824        0.7347      
    n     25                     32            46          
          UBERON:zone of skin of abdomen UBERON:zone of skin of foot
    delta 0.7775                         0.7491                     
    n     12                             64                         
          UBERON:zone of skin of hand UBERON:zone of skin of head
    delta 0.7431                      0.6877                     
    n     64                          32                         
          UBERON:zone of skin of knee UBERON:zone of skin of outer ear
    delta 0.7629                      0.709                           
    n     41                          58                              

    Chance corrected within-group agreement A: 0.002724 
    Based on observed delta 0.7281 and expected delta 0.7301 

    Significance of delta: 0.015 
    Based on  999  permutations

::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          UBERON:ear canal UBERON:feces UBERON:glans penis UBERON:hair
    delta 0.7398           0.7219       0.8002             0.7346     
    n     13               43           7                  14         
          UBERON:labia minora UBERON:mouth UBERON:nose UBERON:nostril
    delta 0.7877              0.7013       0.7375      0.7229        
    n     6                   14           14          28            
          UBERON:nostrils UBERON:skin of arm UBERON:skin of finger
    delta 0.7574          0.7163             0.727                
    n     18              26                 28                   
          UBERON:skin of forearm UBERON:tongue UBERON:urine
    delta 0.7179                 0.754         0.7279      
    n     25                     32            46          
          UBERON:zone of skin of abdomen UBERON:zone of skin of foot
    delta 0.6832                         0.7418                     
    n     12                             64                         
          UBERON:zone of skin of hand UBERON:zone of skin of head
    delta 0.7462                      0.7315                     
    n     64                          32                         
          UBERON:zone of skin of knee UBERON:zone of skin of outer ear
    delta 0.6993                      0.7187                          
    n     41                          58                              

    Chance corrected within-group agreement A: 0.0002254 
    Based on observed delta 0.7299 and expected delta 0.7301 

    Significance of delta: 0.407 
    Based on  999  permutations

::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          UBERON:ear canal UBERON:feces UBERON:glans penis UBERON:hair
    delta 0.7443           0.7066       0.691              0.7308     
    n     13               43           7                  14         
          UBERON:labia minora UBERON:mouth UBERON:nose UBERON:nostril
    delta 0.7075              0.7792       0.6631      0.7378        
    n     6                   14           14          28            
          UBERON:nostrils UBERON:skin of arm UBERON:skin of finger
    delta 0.7474          0.7207             0.7274               
    n     18              26                 28                   
          UBERON:skin of forearm UBERON:tongue UBERON:urine
    delta 0.7342                 0.7486        0.7283      
    n     25                     32            46          
          UBERON:zone of skin of abdomen UBERON:zone of skin of foot
    delta 0.7187                         0.7307                     
    n     12                             64                         
          UBERON:zone of skin of hand UBERON:zone of skin of head
    delta 0.7178                      0.752                      
    n     64                          32                         
          UBERON:zone of skin of knee UBERON:zone of skin of outer ear
    delta 0.7505                      0.7419                          
    n     41                          58                              

    Chance corrected within-group agreement A: -0.00158 
    Based on observed delta 0.7313 and expected delta 0.7301 

    Significance of delta: 0.915 
    Based on  999  permutations

The p-values from the last two tests are very large, indicating that there isn't
significant clustering (this is what we would expect for our shuffled data). The
first test has a smallish p-value of 0.015, but this may be able to be thrown
out due to a bad shuffling of the data (this is why we are doing three shuffled
tests).

Keyboard
^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `HOST_SUBJECT_ID` category as a positive control. We expect
there to be significant clustering on host subjects due to previous analysis
done on the keyboard study dataset.

**Command:** ::

    R --slave --args -d datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID < r/examples/mrpp.r

**Results:**

The following output file is created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          F1     L1    L3    M1     M2     M3     M9    R1    U1    U2    U3   
    delta 0.6344   NaN   NaN 0.5936 0.4754 0.5614 0.529   NaN   NaN   NaN   NaN
    n     3      1     1     2      40     33     31    1     1     1     1    

    Chance corrected within-group agreement A: 0.1407 
    Based on observed delta 0.5232 and expected delta 0.6089 

    Significance of delta: 0.001 
    Based on  999  permutations

The p-value of 0.001 indicates that samples taken from different hosts
are significantly different (i.e. there is clustering). The observed value of
the A statistic also confirms this because it is not sitting around zero. This
is a result that we would expect.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `HOST_SUBJECT_ID`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    R --slave --args -d datasets/keyboard/unweighted_unifrac_dm_shuffled_1.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID < r/examples/mrpp.r
    R --slave --args -d datasets/keyboard/unweighted_unifrac_dm_shuffled_2.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID < r/examples/mrpp.r
    R --slave --args -d datasets/keyboard/unweighted_unifrac_dm_shuffled_3.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID < r/examples/mrpp.r

**Results:**

The following output files are created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          F1     L1    L3    M1     M2    M3     M9     R1    U1    U2    U3   
    delta 0.5839   NaN   NaN 0.5409 0.615 0.6074 0.6031   NaN   NaN   NaN   NaN
    n     3      1     1     2      40    33     31     1     1     1     1    

    Chance corrected within-group agreement A: 0.002931 
    Based on observed delta 0.6071 and expected delta 0.6089 

    Significance of delta: 0.259 
    Based on  999  permutations

::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          F1     L1    L3    M1     M2     M3     M9     R1    U1    U2    U3   
    delta 0.6525   NaN   NaN 0.4728 0.5983 0.6186 0.6115   NaN   NaN   NaN   NaN
    n     3      1     1     2      40     33     31     1     1     1     1    

    Chance corrected within-group agreement A: 0.002417 
    Based on observed delta 0.6074 and expected delta 0.6089 

    Significance of delta: 0.308 
    Based on  999  permutations

::

  Call:
  mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

  Dissimilarity index: 
  Weights for groups:  n 

  Class means and counts:

        F1    L1    L3    M1     M2     M3     M9     R1    U1    U2    U3   
  delta 0.585   NaN   NaN 0.6644 0.6087 0.6028 0.6186   NaN   NaN   NaN   NaN
  n     3     1     1     2      40     33     31     1     1     1     1    

  Chance corrected within-group agreement A: -0.002009 
  Based on observed delta 0.6101 and expected delta 0.6089 

  Significance of delta: 0.631 
  Based on  999  permutations

The p-values from the three tests are all very large, indicating that there is
not significant clustering, which is what we would expect from using shuffled
distance matrices. The three A statistics are sitting around zero as well, which
also confirms the lack of clustering.

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~
**Description:**

This test uses the `CurrentlyWet` category as a positive control. We expect
there to be significant clustering on this category due to previous analysis
done on the Glen Canyon dataset.

**Command:** ::

    R --slave --args -d datasets/glen_canyon/unweighted_unifrac_dm.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet < r/examples/mrpp.r

**Results:**

The following output file is created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          No     Yes   
    delta 0.5012 0.4908
    n     79     15    

    Chance corrected within-group agreement A: 0.1083 
    Based on observed delta 0.4996 and expected delta 0.5603 

    Significance of delta: 0.001 
    Based on  999  permutations

The p-value of 0.001 indicates that samples taken from wet and dry environments
are significantly different (i.e. there is clustering), which is what we'd
expect. The A statistic also confirms this result because it isn't very close to
zero.

Test 2
~~~~~~
**Description:**

This test uses three shuffled distance matrices and the `CurrentlyWet`
category to perform three negative control tests. Since the labels of the
distance matrices are shuffled, we don't expect to see clustering any more on
this category.

**Command:** ::

    R --slave --args -d datasets/glen_canyon/unweighted_unifrac_dm_shuffled_1.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet < r/examples/mrpp.r
    R --slave --args -d datasets/glen_canyon/unweighted_unifrac_dm_shuffled_2.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet < r/examples/mrpp.r
    R --slave --args -d datasets/glen_canyon/unweighted_unifrac_dm_shuffled_3.txt -m datasets/glen_canyon/map_25Jan2012.txt -c CurrentlyWet < r/examples/mrpp.r

**Results:**

The following output files are created: ::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          No     Yes   
    delta 0.5553 0.5852
    n     79     15    

    Chance corrected within-group agreement A: 0.0002716 
    Based on observed delta 0.5601 and expected delta 0.5603 

    Significance of delta: 0.362 
    Based on  999  permutations

::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          No    Yes   
    delta 0.558 0.5744
    n     79    15    

    Chance corrected within-group agreement A: -0.0006878 
    Based on observed delta 0.5606 and expected delta 0.5603 

    Significance of delta: 0.554 
    Based on  999  permutations

::

    Call:
    mrpp(dat = as.dist(qiime.data$distmat), grouping = qiime.data$map[[opts$category]]) 

    Dissimilarity index: 
    Weights for groups:  n 

    Class means and counts:

          No     Yes   
    delta 0.5629 0.5522
    n     79     15    

    Chance corrected within-group agreement A: -0.001617 
    Based on observed delta 0.5612 and expected delta 0.5603 

    Significance of delta: 0.812 
    Based on  999  permutations

The three p-values are very large, indicating that samples taken from wet vs.
dry environments are not significantly different, which is what we would expect.
The A statistics are all near zero, which also indicates that we didn't see
clustering.

References
----------
.. _mrppref1:

[1] R help page for vegan function mrpp

.. _mrppref2:

[2] http://www.jstor.org/stable/1940409

.. _mrppref3:

[3] http://lcai.bol.ucla.edu/programs.html

.. _mrppref4:

[4] http://people.oregonstate.edu/~mccuneb/Chapter24.ppt

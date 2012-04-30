.. _distance_matrix_comparison:

===========================
Comparing Distance Matrices 
===========================

Introduction
------------
This tutorial explains how to use several different distance matrix comparison techniques that are available in `compare_distance_matrices.py <../scripts/compare_distance_matrices.html>`_. All of the currently available comparison techniques are based on the Mantel test, which is a non-parametric statistical method that computes the correlation between two distance matrices. In addition to this statistical method, QIIME also provides the partial Mantel test and Mantel correlogram. Each of these methods will be described in greater detail below.

One common application of distance matrix comparison techniques is to determine if correlation exists between a community distance matrix (e.g. UniFrac distance matrix) and a second matrix derived from a gradient (e.g. difference in pH, temperature, or geographical location). For example, one might be interested in seeing if communities that are at dissimilar pH levels are more different from one another than communities that are at very similar pH levels. If so, this indicates positive gradient correlation. To create a gradient distance matrix from a continuous category in your mapping file, please refer to the `distance_matrix_from_mapping.py <../scripts/distance_matrix_from_mapping.html>`_ script documentation.

Please note that this tutorial does not attempt to cover every possible option that can be used in the distance matrix comparison script. Instead, it attempts to provide useful examples to give you an idea of how to use these statistical methods in your own analysis, as well as customize some of the output to your liking. For a complete listing of the available options, please refer to the `compare_distance_matrices.py <../scripts/compare_distance_matrices.html>`_ script documentation.

.. _inputfiles:

Input Files
^^^^^^^^^^^
You can obtain the files used in this tutorial `here <https://s3.amazonaws.com/s3-qiime_tutorial_files/88_soils/88_soils.zip>`_. The files are taken from the study performed by (Lauber et al., 2009) where 88 soil samples were collected at various regions around the world. pH was recorded for each of the soil samples. Using `distance_matrix_from_mapping.py <../scripts/distance_matrix_from_mapping.html>`_, we created a distance matrix containing differences in pH between each pair of samples.

.. _manteltest:

Mantel Test
-----------
The Mantel test tests the correlation between two distance matrices. It is non-parametric and computes the significance of the correlation through permutations of the rows and columns of one of the input distance matrices. The test statistic is the Pearson product-moment correlation coefficient `r`. `r` falls in the range of -1 to +1, where being close to -1 indicates strong negative correlation and +1 indicates strong positive correlation. An `r` value of 0 indicates no correlation.

To illustrate the use of the Mantel test, we will determine if there is significant correlation between the unweighted UniFrac distance matrix and the pH distance matrix. Run the following command:

[insert command]

This command will create a new output directory named :file:`tutorial_output`, which will contain a single text file called :file:`mantel_results.txt`. Open up :file:`mantel_results.txt` to see the results of the test:

[insert results here]

The Mantel r statistic of [foo] indicates that there is strong positive correlation between the distance and pH matrices. The p-value of [bar] indicates that our results are statistically significant at an alpha of 0.05. We determined the p-value by specifying 999 permutations with the -n option.

Partial Mantel Test
-------------------
The partial Mantel test is used to estimate the correlation between two matrices, A and B, while controlling for the effect of a control matrix C. The partial Mantel test is a first-order correlation analysis that utilizes three distance (dissimilarity) matrices. This builds on the simple Mantel test by adding a third "control" matrix. The goal is to test the correlation between matrices A and B while controlling the effect of a third matrix C, in order to remove spurious correlations. Only the first distance matrix is permuted so that the correlation structure between the second and third matrices is kept constant.

To illustrate the use of the partial Mantel test, we will determine if there is significant correlation between the unweighted and weighted UniFrac distance matrices, using the pH distance matrix as the control matrix. Run the following command:

[insert command]

This command will create a new output directory named :file:`tutorial_output`, which will contain a single text file called :file:`partial_mantel_results.txt`. Open up :file:`partial_mantel_results.txt` to see the results of the test:

[insert results here]

The Mantel r statistic of [foo] indicates that there is strong positive correlation between the unweighted and weighted distance matrices while controlling for differences in pH. The p-value of [bar] indicates that our results are statistically significant at an alpha of 0.05.

Mantel Correlogram
------------------
Mantel correlogram is a method that tests whether there is correlation between two distance matrices by examining the correlation between matrices for each distance class. Mantel correlogram performs a Mantel test on each distance class and generates a correlogram with distance classes on the x-axis and their corresponding Mantel test statistic on the y-axis. Thus, Mantel correlogram allows you to see where the correlation exists between the two matrices by providing a higher-resolution view than a traditional Mantel test. For example, you might want to see if there is correlation between a Unifrac distance matrix and a spatial distance matrix. Mantel Correlogram will let you see what the correlation is at different ranges of spatial distances (distance classes).

The Mantel correlogram method computes a Mantel statistic for each geographic distance class that can be derived from the input. It tests for significance of genetic/community distance versus geographic distance (or some other type of distance). The null hypothesis that is tested is that there is no association of geographic distance to community distance for each distance class.

Sturge’s rule is used to determine how many distance classes to use based on the number of pairwise comparisons you have. These distance classes can be thought of as bins (as used in histograms). For each distance class, a Mantel test is performed and a Mantel statisic is computed. A corrected p-value (i.e. Bonferroni, FDR, Holm, etc.) is also computed for each test. The results of this method are usually visualized in a correlogram, which is a graph with the geographic distance classes on the x-axis and the Mantel statistics on the y-axis.

This method is very similar to the Mantel method, so the resulting Mantel statistics can be interpreted in the same way as you would for a traditional Mantel test (i.e. a positive value indicates positive spatial correlation). p-values are obtained in the same way as well (i.e. through permutations).

To illustrate the use of the Mantel correlogram test, we will determine if there is significant correlation between the unweighted UniFrac distance matrix and the pH distance matrix. Run the following command:

[insert command]

This command will create a new output directory named :file:`tutorial_output`, which will contain two files called :file:`mantel_correlogram_results.txt` and :file:`mantel_correlogram.pdf`. Open up :file:`mantel_correlogram_results.txt` to see the results of the test:

[insert results here]

[write something about the results text file]

Open up :file:`mantel_correlogram.pdf` to view the Mantel correlogram:

[insert correlgoram image here]

The correlogram is a visual representation of the results in the output text file. Points that are filled in (black) are statistically significant at an alpha of 0.05 (this is configurable as a parameter to the script). Points that are not filled in (white) are not statistically significant at the specified alpha level. By examining the correlogram, we see that positive correlation exists at closer distance classes, while the strength of the correlation decreases as the distance classes increase.

References
----------
Pyrosequencing-based assessment of soil pH as a predictor of soil bacterial community structure at the continental scale.  Lauber CL, Hamady M, Knight R, Fierer N.  Appl Environ Microbiol. 2009 Aug;75(15):5111-20.

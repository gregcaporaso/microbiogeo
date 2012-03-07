===========================================
Microbiogeo Method Analysis Progress Report
===========================================

Introduction
------------
The purpose of this document is to summarize the progress made on our
statistical methods research and testing. We have successfully processed a
majority of the statistical methods assigned, though there are still a few
methods that have not been completely tested. Levi is currently in contact 
with the author of MultiCoLA to obtain all of the necessary scripts and expects them in a week.
Damien is finishing up PERMDISP after finding out the original
implementation that he picked would not work for our tests. He is currently
testing an implementation of PERMDISP in R. A few of the lower-priority methods
were recently assigned and are currently being evaluated.

Testing Strategy
----------------
The statistical methods can be categorized into two groups: clustering and
correlation. We tested the clustering methods with the keyboard, Glen Canyon,
and whole body datasets. Each of the clustering methods were tested with three
Unifrac distance matrices with shuffled labels, which served as negative control
tests. For the keyboard study, we tested clustering on individuals as a
positive control. For the whole body study, we used body site as a
positive control and sex as a negative control (in addition to the shuffled
distance matrices negative control test described above). For the Glen Canyon
study, we used the “currently wet” category as a positive control.

We tested the correlation methods with the 88 soils, keyboard, and Glen Canyon
datasets. The correlation methods accept either two or three distance matrices
as input and compute the correlation between the first two matrices. Negative
control tests were performed using three Unifrac distance matrices with shuffled
labels. For the 88 soils study, we created a second distance matrix containing
differences in pH and used this as a positive control. We used a distance matrix
containing differences in latitude for an additional negative control test. For
the keyboard study, we filtered the original Unifrac distance matrix to include
only samples taken from keys, and then we created a second distance matrix of
the euclidean distances between keys. An additional third distance matrix of
median distances between individuals was also created for testing Partial
Mantel. For the Glen Canyon study, we created a second distance matrix of
differences in “time since last submerged” to use as a positive control.

Summary of Results
------------------
The following table summarizes the testing results that were obtained for each
of the higher-priority statistical methods, as well as a couple of the
lower-priority ones. The methods have been ordered by their apparent
usefulness, which is determined by how well they responded to positive and
negative controls. We plan on implementing the methods in the order presented in
this table and will continue working on testing the remaining lower-priority
methods and implementing them in an appropriate order as they are finished and
conclusions are made regarding their usefulness.

:note: Click on a method name to see more detailed information about the method itself and the specific testing results.

+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| Method Name                                                           | Method Type               | Meaningful Results | Notes                                                                             |
+=======================================================================+===========================+====================+===================================================================================+
| `LSA <lsa.html>`_                                                     | Correlation               | N/A                | Failed on data provided by developer                                              |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Mantel <mantel.html>`_                                               | Correlation               | Yes                | Responds well to positive and negative controls                                   |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Mantel Correlogram <mantel_correlogram.html>`_                       | Correlation               | Yes                | Responds well to positive and negative controls                                   |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Moran's I <morans_i.html>`_                                          | Correlation               | Yes                |                                                                                   |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `RDA <rda.html>`_                                                     | Clustering/Ordination     | Yes                | Responds well to positive and negative controls                                   |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `ANOSIM <anosim.html>`_                                               | Clustering                | Unsure             | Low specificity**                                                                 |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Adonis <adonis.html>`_                                               | Clustering                | Unsure             | Low specificity**                                                                 |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `BEST <best.html>`_                                                   | Clustering                | Unsure             | Shuffled matrices neg. control didn't work                                        |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `MRPP <mrpp.html>`_                                                   | Clustering                | Unsure             | Low specificity**                                                                 |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `MultiCoLA <MultiCoLA.html>`_                                         | Clustering                | N/A                | Missing some scripts, in contact with author                                      |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Partial Mantel <partial_mantel.html>`_                               | Correlation               | Unsure             | Needs better positive control test*                                               |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `PERMANOVA <permanova.html>`_                                         | Clustering                | Unsure             | Low specificity**                                                                 |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `PERMDISP <permdisp.html>`_                                           | Clustering                | N/A                | Testing in progress                                                               |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+
| `Repeated Measures PERMANOVA <repeated_measures_permanova.html>`_     | Clustering                | Unsure             | Time fields in Glen Canyon (or any used as proxy) do not yield meaningful results |
+-----------------------------------------------------------------------+---------------------------+--------------------+-----------------------------------------------------------------------------------+

\* It is hard to judge whether Partial Mantel gives biologically meaningful
results because we don't have a good positive control test. We tested it with
the keyboard dataset to see if there was spatial correlation between the Unifrac
distance matrix and the euclidean distance matrix (with effects partialled out
by a third distance matrix of median Unifrac distances between individuals). The
results of the test didn't indicate significant spatial correlation, but we
don't really know what to expect for this test because it has never been done
before. We'd like to see significant spatial correlation, but we don't actually
have a way to verify whether Partial Mantel is functioning correctly or giving
meaningful results. It does, however, give the same results that Mantel and
Mantel correlogram did for this test (the third distance matrix was not used for
those tests, though).

\** Several of the clustering methods (i.e. ANOSIM, Adonis, PERMANOVA, and MRPP)
seem to be very sensitive to clustering. These methods were tested on the whole
body dataset using the sex category to determine the grouping of samples. We
thought this would be a good negative control to test as we did not expect to
see significant clustering on sex. Much to our surprise, the methods returned
very small p-values (around 0.001 with 999 permutations). We also tested the
methods on other categories that shouldn't indicate clustering (e.g. DOB for the
overview tutorial, Hour or Day for Glen Canyon) and still received significant
p-values. Greg suggested that we shuffle the individual's sexes in the whole
body mapping file and try the tests again. Adonis, PERMANOVA, and MRPP still
gave p-values of 0.001, but ANOSIM gave a large p-value (0.201). Thus, ANOSIM
might be detecting very weak clustering on sex, but it is hard to say what the
others are detecting. When we tested these methods using the three shuffled
distance matrices as negative controls, they computed large p-values, which is
what we expected to see.

Some sources say that ANOSIM, MRPP, and \*ANOVA sometimes detect differences in
spread (i.e. variability) within groups and report that the groups are
significantly different. We also came across warnings saying that MRPP can be
sensitive to outliers. Thus, our strange results may be related to these issues.
The next step is to test these methods on simulated data in order to help us
decide whether these methods are worth implementing or not.


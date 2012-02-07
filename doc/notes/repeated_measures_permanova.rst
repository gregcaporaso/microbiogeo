========================================================
Repeated Measures PERMANOVA Statistical Method Reference
========================================================

Introduction
------------
Repeated measures PERMANOVA is an application of traditional PERMANOVA to
time series data. It is useful if communities have been sampled multiple times
in a person or a plot over time.

This statistical method is very elusive and hard to track down. There doesn't
really seem to be a consensus regarding how time should be treated in these
types of data. Many people simply use traditional PERMANOVA, but this ignores
the correlation between space/time, thus inflating p-values
[:ref:`2 <rmpermref2>`]. Others use time as a fixed factor
[:ref:`1 <rmpermref1>`], while others still use time as a random factor
[:ref:`4 <rmpermref4>`]. [:ref:`5 <rmpermref5>`] is a good reference for an
explanation of fixed versus random factors.

Existing Implementations
------------------------
The only existing implementation of repeated measures PERMANOVA I could find was
in PRIMER v6 with PERMANOVA+ extension [:ref:`3 <rmpermref3>`].

PRIMER must be purchased, and it looks pretty expensive. The adonis function in
R's vegan package, which is used to do normal PERMANOVAs, might be used to
create our own repeated measures PERMANOVA implementation, though we must be
careful in how we set up the data and inputs to the adonis function in order to
account for time correctly. [:ref:`2 <rmpermref2>`] presents an R script that
supposedly runs repeated measures PERMANOVA over data taken at three timepoints.
We may be able to use this implementation as a guide to our own implementation
of this statistical method.

The rest of the document will explain how to set up your system to run the
repeated measures PERMANOVA example that I found online. You may obtain the
script from the Qiimeutils repository under
:file:`microbiogeo/r/repeated_measures_permanova.r`.

System Setup and Required Dependencies
--------------------------------------
:note: The following instructions have been tested on 64-bit Linux Mint (essentially Debian). However, they `should` work across different Linux distros and on Macs, though some commands may need to be tweaked, or different package names might have to be used. The instructions assume you use bash as your shell.

The first step is to install R. The following command downloaded and installed R
(for me, it was R version 2.13.1): ::

    sudo apt-get install r-base

Next, you must install the vegan package in R. Run the following commands: ::

    sudo R
    install.packages("vegan")
    q()

The install process for vegan will prompt you to choose a mirror to download
vegan from. Other than that, it is completely automated. On my system, I ended
up with vegan version 2.0-2.

Input Files
-----------
The script requires a distance matrix file (i.e. the result of
beta_diversity.py) and a metadata mapping file. I used the unweighted Unifrac
distance matrix and a modified mapping file from the QIIME overview tutorial.
You can get the distance matrix
:download:`here <../downloads/overview_unweighted_unifrac_dm.txt>` and the
modified mapping file
:download:`here <../downloads/Fasting_Map_rep_meas_perm.txt>`. The mapping file
was modified to include a "Time" column. `The values in the "Time" column are
completely contrived to use as example input to this script.`

Execute the following command to run the script: ::

    R --slave --args -d overview_unweighted_unifrac_dm.txt -m Fasting_Map_rep_meas_perm.txt -c Time < r/repeated_measures_permanova.r


Output Files
------------
There are no output files from this script as it prints all of its information
to stdout. It prints the true R2 value (i.e. no permutations) and then prints
the p-value based on the permutations that it computes. When I ran this example
a few times, the p-value was sitting around 0.1, which may indicate that the
"Time" category is a good indicator of variability in the samples. More
extensive testing will have to be done on real time series data.

References
----------
.. _rmpermref1:

[1] http://aspenface.mtu.edu/pdfs/Andrew%20and%20Lilleskov.pdf

.. _rmpermref2:

[2] http://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html#more

.. _rmpermref3:

[3] http://www.cfc.umt.edu/biogeochemistry/Pdfs/Nemergut_SBB_2010.pdf

.. _rmpermref4:

[4] http://www.talkstats.com/showthread.php/16088-PERMANOVA-in-R-adonis-function

.. _rmpermref5:

[5] http://www.jerrydallal.com/LHSP/fixran.htm

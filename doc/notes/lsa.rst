.. _lsa:

======================================================
Local Similarity Analysis
======================================================

Introduction
-------------
LSA [:ref:`1 <lsaref1>`] a python package for local similarity alignment of 
sequence data and constructing time-shifted interaction map

LSA attempts to capture the time-dependent associations (possibly time-shifted) 
between microbes and between microbe and environmental factors

Unfortunately the LSA software seems to be, at best, unstable. In my tests with
the example data provided by the developer, not only were the results excruciatingly
slow in their generation, but a majority of the "p-values" were reported as NaN. 
The algorithm took roughly one hour to perform a 99 permuatation analysis, the default
and recommended number of permutations is 1000, taking roughly 14 hours to run. 

Although the developer is somewhat responsive to questions, the general documentation
on the method is exceedingly sparse, what is available seems to have been thrown together
in haste. I feel, that, although this method does have potential it is far too underdeveloped
both in software and documentation at this time to be feasible for inclusion in the project
package.

Existing Implementations
------------------------
The only implementation that I am aware of is that found on the LSA repo[:ref:`2 <lsaref2>`]

System Setup and Required Dependencies
--------------------------------------
Currently the package works only for Ubuntu Linux and Mac OS X 10.6+ platforms.
The following was all tested in Mac OS X 10.7.3 with a local install of QIIME 1.4.0-dev.

In order to use LSA it must be installed. The LSA repository website contains
directions for installation [:ref:`2 <lsaref2>`]

LSA dependencies(several, so be sure to read the install instructions):

* SWIG
* Nose
* GFortran
* scipy
* numpy (already included if you have QIIME)

Additionally, there are virtual disk images that can be used inside VirtualBox,
if you so choose. Directions for this are also on the LSA repo.

To run the analysis, I executed the following on the example dataset: ::

  lsa_compute newdata.txt test.out -p 99
  lsa_compute newdata.txt test.out -p 1000 

Input Files
-----------

the input data file is an m by (r * s), tab delimited text file(QIIME map file?);
top left cell start with '#' to mark this is the header line;m is number of 
variables, r is number of replicates and s is the number of time spots; 
first row: #header s1r1 s1r2 s2r1 s2r2; 
second row: x ?.?? ?.??  ?.?? ?.??; for a 1 by (2*2) data

Output Files
------------

An excerpt from 'test.out': ::

  X  	Y	   LS   	lowCI  	upCI	 Xs Ys Len	Delay	P	      PCC	   Ppcc	   SPCC	  Pspcc	  Q	      Qpcc
  402	399	0.2304	0.2304	0.2304	1	1	 2	   0	  1.0000	nan	   1.0000  	nan	  nan	    1.0000	1.0000
  405	399	0.3370	0.3370	0.3370	1	1	 4	   0	  1.0000	nan	   1.0000  	nan	  1.0000  1.0000	1.0000
  405	402	0.1771	0.1771	0.1771	4	1	 1	   3	  1.0000	0.2098 0.7902	  nan	  nan	    1.0000	1.0000
  408	399	0.7083	0.7083	0.7083	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000  1.0000	1.0000
  408	402	0.2304	0.2304	0.2304	1	1	 2	   0	  0.3737	nan	   1.0000	  nan	  nan	    1.0000	1.0000
  408	405	0.3370	0.3370	0.3370	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000  1.0000	1.0000
  411	399	0.7083	0.7083	0.7083	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000  1.0000	1.0000
  411	402	0.2304	0.2304	0.2304	1	1	 2	   0	  0.4545	nan	   1.0000	  nan	  nan	    1.0000	1.0000
  411	405	0.3370	0.3370	0.3370	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000	1.0000	1.0000
  411	408	0.7083	0.7083	0.7083	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000	1.0000	1.0000
  414	399	0.7083	0.7083	0.7083	1	1	 4	   0	  1.0000	nan	   1.0000	  nan	  1.0000	1.0000	1.0000
  


References
----------
.. _lsaref1:

[1] http://meta.usc.edu/softs/lsa

.. _lsaref2:

[1] https://bitbucket.org/charade/lsa

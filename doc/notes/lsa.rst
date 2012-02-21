.. _lsa:

======================================================
Local Similarity Analysis
======================================================

Introduction
------------

LSA [:ref:`1 <lsaref1>`] 


Existing Implementations
------------------------
The only implementation that I am aware of is that found on the LSA repo[:ref:`2 <lsaref2>`]

System Setup and Required Dependencies
--------------------------------------
Currently the package works only for Ubuntu Linux and Mac OS X 10.6+ platforms.
The folowwing was all tested in Mac OS X 10.7.3 with a local install of QIIME 1.4.0-dev.

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


Input Files
-----------

the input data file is an m by (r * s), tab delimited text file(QIIME map file?);
top left cell start with '#' to mark this is the header line;m is number of 
variables, r is number of replicates and s it number of time spots; 
first row: #header s1r1 s1r2 s2r1 s2r2; 
second row: x ?.?? ?.??  ?.?? ?.??; for a 1 by (2*2) data


Output Files
------------


References
----------
.. _lsaref1:

[1] http://meta.usc.edu/softs/lsa

.. _lsaref2:

[1] https://bitbucket.org/charade/lsa

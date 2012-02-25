==============================================================
MultiCoLA
==============================================================

Introduction
------------


Existing Implementations
------------------------


System Setup and Required Dependencies
--------------------------------------

:note: The following instructions have been tested on 64-bit Fedora (essentially Redhat) using Python 2.7.2. However, they `should` work across different linux distros and on Macs. The instructions assume you use bash as your shell.

First, your system must have R installed, you can get R by running the command: ::

	yum install R
	
Then launch R by running the command: ::

	R
	
You should now be in the R console

Now install vegan and mass with the commands, this will take a long time: ::

	install.packages("vegan", dependencies=TRUE)
	install.packages("MASS", dependencies=TRUE)
	
Change working directory to where the input files and scripts are: ::
	
	setwd("/home/aragorn/MultiCoLA.1.3")
	
Useful commands
---------------
	
Read and store input file in a variable: ::

	M<-read.table("input.txt",header=TRUE,row.names=1)
	
Load script(s): ::

	source("script.r")
	
Run script and store output in a variable: ::

	results<-script(M)
	
	
My progress understanding MultiCoLA
-----------------------------------
Ran the commands: :
	M<-read.table("input.txt",header=TRUE,row.names=1)
	source("taxa.pooler.1.3.r")
	all_taxa_pooled<-taxa.pooler(M)
	
Answered the questions:

	Number of samples? (e.g. 16)... 16
	Number of taxonomic levels? (e.g. phylum+class+order+family+genus=5)... 5
	Presence/absence tables as output? (y/n) y
	Output as text files? (y/n)... y
	
The following files were created:

:download:`Phylum.matrix.txt <../downloads/MultiCoLA/Phylum.matrix.txt>`
:download:`Class.matrix.txt <../downloads/MultiCoLA/Class.matrix.txt>`	
:download:`Order.matrix.txt <../downloads/MultiCoLA/Order.matrix.txt>`
:download:`Family.matrix.txt <../downloads/MultiCoLA/Family.matrix.txt>`
:download:`Genus.matrix.txt <../downloads/MultiCoLA/Genus.matrix.txt>`
:download:`MultiCoLA.RData <../downloads/MultiCoLA/MultiCoLA.RData>`
:download:`OTUs_completeDS.matrix.txt <../downloads/MultiCoLA/OTUs_completeDS.matrix.txt>`
:download:`OTUs_wholeDS.matrix.txt <../downloads/MultiCoLA/OTUs_wholeDS.matrix.txt>`


Ran the commands: ::

	source("COtables.1.3.r")
	truncated.DS.i<-COtables(all_taxa_pooled[[i]], Type="ADS",typem="dominant")
	
Got Error in rbind(ODS, apply(ODS, 2, sum)) : object 'i' not found.
Finding out what "i" should be.

	
	



Input Files
-----------


Output Files
------------


Testing Results
---------------
This section will describe different tests that were run on the ANOSIM script.
These tests will use empirical data from one of the several datasets that the
team has access to. These data files will not be included for download due to
their (usually) large size. Unless otherwise noted, the data files that were
used can be found under the datasets directory.

Whole Body
^^^^^^^^^^
Test 1
~~~~~~
**Description:**


**Command:** ::

**Results:**

The following output file is created: ::



Test 2
~~~~~~

Test 3
~~~~~~

Keyboard
^^^^^^^^

Test 1
~~~~~~

Test 2
~~~~~~

Test 3
~~~~~~

Glen Canyon
^^^^^^^^^^^

Test 1
~~~~~~

Test 2
~~~~~~

Test 3
~~~~~~

References
----------

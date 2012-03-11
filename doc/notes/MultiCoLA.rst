==============================================
Multivariate Cutoff Level Analysis (MultiCoLA)
==============================================

Introduction
------------
Many microbiological studies have a problem with `rare` sequences in the raw data. The problem is: how do you define rare? Are the rare sequences organic, or simply an artifact of the sequencer?  Do the rare sequences have an effect on your community? MultiCoLA is a proposed solution to this problem.  The idea is that you take your dataset with no rarity cutoff, then you apply increasingly rigorous rarity cutoffs.  You can then study the results to see where significat changes took place, and where the data becomes unpredictable. When preforming the rarity cutoffs the rare sequences are discarded. MultiCoLA allows the definition of OTU to be changed, you can say that if a sequence is off by one base it is in a new OTU, or they can be the same OTU if they are in the same phylum, genus, ect.

Existing Implementations
------------------------
The only implementation seems to be by the authors who suggested this method.  The implementation is in R and is a collection of scripts rather than a all encompassing script.


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
	
	
Our progress understanding MultiCoLA
------------------------------------
**Ran the commands:** ::
	M<-read.table("input.txt",header=TRUE,row.names=1)
	source("taxa.pooler.1.3.r")
	all_taxa_pooled<-taxa.pooler(M)
	
**Answered the questions:**

	Number of samples? (e.g. 16)... 16
	Number of taxonomic levels? (e.g. phylum+class+order+family+genus=5)... 5
	Presence/absence tables as output? (y/n) y
	Output as text files? (y/n)... y
	
**The following files were created:**
:download:`Phylum.matrix.txt <../downloads/MultiCoLA/Phylum.matrix.txt>`
:download:`Class.matrix.txt <../downloads/MultiCoLA/Class.matrix.txt>`	
:download:`Order.matrix.txt <../downloads/MultiCoLA/Order.matrix.txt>`
:download:`Family.matrix.txt <../downloads/MultiCoLA/Family.matrix.txt>`
:download:`Genus.matrix.txt <../downloads/MultiCoLA/Genus.matrix.txt>`
:download:`MultiCoLA.RData <../downloads/MultiCoLA/MultiCoLA.RData>`
:download:`OTUs_completeDS.matrix.txt <../downloads/MultiCoLA/OTUs_completeDS.matrix.txt>`
:download:`OTUs_wholeDS.matrix.txt <../downloads/MultiCoLA/OTUs_wholeDS.matrix.txt>`


**Ran the commands:** ::

	source("COtables.1.3.r")
        truncated.DS.phylum<-COtables(all_taxa_pooled[[1]], Type="ADS",typem="dominant")
        truncated.DS.class<-COtables(all_taxa_pooled[[2]], Type="ADS",typem="dominant")
        truncated.DS.order<-COtables(all_taxa_pooled[[3]], Type="ADS",typem="dominant")
        truncated.DS.family<-COtables(all_taxa_pooled[[4]], Type="ADS",typem="dominant")
        truncated.DS.genus<-COtables(all_taxa_pooled[[5]], Type="ADS",typem="dominant")
        truncated.DS.OTUwholeDS<-COtables(all_taxa_pooled[[6]], Type="ADS",typem="dominant")
        source("cutoff.impact.1.3.r")
        corr.all<-cutoff.impact(all_taxa_pooled,Type="ADS",corcoef="spearman",typem="dominant")

**Answered the question:**

        Details of the NMDS calculations? (y/n)...      n

**Was warned:** ::

        There were 50 or more warnings (use warnings() to see the first 50)
        
**Ran the commands:** ::

        source("cutoff.impact.fig.1.3.r")
        output.all<-cutoff.impact.fig(corr.all)

**Answered the questions:**

        Output as text files? (y/n)...  y

        Plot the results? (y/n)...      y

**The files were created:**
:download:`abundance.txt <../downloads/MultiCoLA/abundance.txt>`
:download:`non-par.correlation.txt <../downloads/MultiCoLA/non-par.correlation.txt>`
:download:`procrustes.txt <../downloads/MultiCoLA/procrustes.txt>`

**This graph was displayed:**

.. image:: ../images/MultiCoLA/graph.png
      :align: center

**Ran the commands:** ::

        ENV<-read.table("env.txt",header=TRUE,row.names=1)
        source("VP.COL.1.3.r")
        VP.1.taxa<-VP.COL(all_taxa_pooled,ENV,Type="ADS")

**Answered the question:**

        Output as text files? (y/n)...  y
        Plot the results? (y/n)...      y
 
**The files were created:**
:download:`Phylum.sum.adjRsq.txt <../downloads/MultiCoLA/Phylum.sum.adjRsq.txt>`
:download:`Phylum.VarPart.txt <../downloads/MultiCoLA/Phylum.VarPart.txt>`
:download:`Class.sum.adjRsq.txt <../downloads/MultiCoLA/Class.sum.adjRsq.txt>`
:download:`Class.VarPart.txt <../downloads/MultiCoLA/Class.VarPart.txt>`
:download:`Order.sum.adjRsq.txt <../downloads/MultiCoLA/Order.sum.adjRsq.txt>`
:download:`Order.VarPart.txt <../downloads/MultiCoLA/Order.VarPart.txt>`
:download:`Family.sum.adjRsq.txt <../downloads/MultiCoLA/Family.sum.adjRsq.txt>`
:download:`Family.VarPart.txt <../downloads/MultiCoLA/Family.VarPart.txt>`
:download:`Genus.sum.adjRsq.txt <../downloads/MultiCoLA/Genus.sum.adjRsq.txt>`
:download:`Genus.VarPart.txt <../downloads/MultiCoLA/Genus.VarPart.txt>`

**This graph was displayed:**

.. image:: ../images/MultiCoLA/graph2.png
      :align: center


**Ran the command:** ::

        VP.1.taxa<-VP.COL(all_taxa_pooled,ENV,Type="ADS")

**Answered the questions:**

        Output as text files? (y/n)...  y
        Plot the results? (y/n)...      y
 
**Ran the commands:** ::
        
        source("corrcoeff.ENV.1.3.r")
        corrcoeff.table.ADS<-matrix(NA,21,5)
        row.names(corrcoeff.table.ADS)<-c(paste("CO_",c(0.01,seq(0.05,0.95,by=0.05),0.99),sep=""))
        colnames(corrcoeff.table.ADS)<-c("Sum",paste("RDA1.",colnames(ENV),sep=""))
        OTU.ADS<-VP.1.taxa[[c(7,3)]]

**Found problem with:** ::

        > SPE<-OTU.ADS[[1]];corrcoeff.table.ADS[1,]<-corrcoeff(SPE,ENV);rm(SPE)
        Error in if (nrow(SPE) < nrow(SPE)) { : argument is of length zero
        ...Researching


Input Files
-----------


Output Files
------------


Testing Results
---------------
This section will describe different tests that were run on the MultiCoLA script.
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

        The command

**Results:**

The following output file is created



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

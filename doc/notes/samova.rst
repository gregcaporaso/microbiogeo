===============================================
SAMOVA - Spatial Analysis of MOlecular VAriance
===============================================

Synopsis
--------

SAMOVA is a statistical method that takes information from samples with similar geographical data and partitions it into distinct groups that are different. A by-product of this process is that it lends itself to identifying these partitions.

The Problem Being Answered
^^^^^^^^^^^^^^^^^^^^^^^^^^

Is there any relationship between these samples when geographic data accounted for?

The Premise
^^^^^^^^^^^

If you're looking to use SAMOVA you need to be interested in trying to identify some relationship between samples when you have the geographic data of the samples as well.

If the following conditions are satisfied then you would be interested in using SAMOVA.

- If you have a windows operating system installed on your computer
- If you have a set of sample data:
    - and this sample data is in an arlequin file file format, (either DNA or Microsat)
    - and you have the sample data's geographical coordinates where they were taken
        - and these coordinates are distinct from all the other coordinate
        - and it's formatted in a .geo file format


Plain English
^^^^^^^^^^^^^


Introduction
------------


Selected Implementation
-----------------------

Of the available implementations there appear to be only one of notable interest. It's called SAMOVA 1.0. When simply performing search queries this is one of the top ranked search results. Furthermore it's returned on several different sets of search queries as well. There is an etoology website that points to SAMOVA 1.0. The best alternative to this was an application called Bayesian Analysis of Population Structure (BAPS).

Because of this SAMOVA 1.0 was selected as the application to be used.

Alternative Implementations
---------------------------

There don't appear to be many alternatives to SAMOVA 1.0. The only mentioned alternative was a program called BAPS, which can be found here:
http://web.abo.fi/fak/mnf/mate/jc/software/baps.html

System Setup and Required Dependencies
--------------------------------------

In order to run SAMOVA you first MUST HAVE WINDOWS XP, WINDOWS VISTA, OR WINDOWS 7.

Once you have verified you have the correct operating system you can download the SAMOVA 1.0 zip file from this link:
http://cmpg.unibe.ch/software/samova/SAMOVA_12_02.zip

If that link doesn't work you can download it from the main website at:
http://cmpg.unibe.ch/software/samova/

Once downloaded you can extract the contents of the file. If you extract the files from the default folder, you should expect a folder name of SAMOVA_12_02.

After that, you will find several files listed in the folder. There are only a few files that you will truly be concerned about. 

Input
-----

  -  A .geo file: This is a tab seperated file that contains the geographic data associated with the samples from the .arp file. THEY MUST BE IN THE SAME ORDER THATTHE ARP FILE IS. 


  -  A .arp file: This is an arlequin file format. It can be generated in multiple forms via the Arlequin program. (Indications from SAMOVA state that it only takes the DNA format and the Microsatellite data.)

When you get the initial default files provided by the SAMOVA 1.0 download there will be two default input files.

The files:

inputfile.arp - This file is the default file provided with the downloaded folder. It provides example usage for the .arp file type

inputfile.geo - This file is the default file provided with the downloaded folder. It provides example usage for the .geo file type

Both of these files must have a corresponding file of the other one so it can be used. In other words, this means there must be BOTH a .arp file and a .geo file. Furthermore, these files must be named the same as well. 

So let's pretend you want to add some sample data from some gut community samples. A good name for these might be gutCommunity. This would mean that their input files would have to be gutCommunity.arp, and gutCommunity.geo.

Output
------

Samova resturns five files back as output.

  -  SAMOVA_results_arlequin.txt: the genetic structure defined by SAMOVA as well as the fixation indices corresponding to this group structure and their significance level evaluated by 1,000 permutations of populations among groups.
 
  -  SAMOVA.log: this file contains all the steps done by SAMOVA and, in case of problems, the location of the problems.

  -  SAMOVA_finalstructure.arp: an arlequin project file created by appending the input arlequin project file with the genetic structure defined by SAMOVA.

  -  SAMOVA_results.ps: this files (eps) can be read with GSview for Windows or Adobe Illustrator 7.0; it contains a map of the sampling points and the barriers between the groups of populations defined by SAMOVA.

  -  Arlequin.log: this file is generated during the computation of the fixation indices corresponding to the genetic structure defined by SAMOVA. It contains all the run-time WARNINGS and ERRORS encountered during this computations. 


Testing Results
---------------

N/A

Results Analysis
----------------

N/A

References
----------

[1] SAMOVA

http://cmpg.unibe.ch/software/samova/

[2] etoology SAMOVA Reference

http://www.etoology.net/index.php/software/genetics/92-samova-10.html


[3] BAPS application

http://web.abo.fi/fak/mnf/mate/jc/software/baps.html

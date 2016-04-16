
This project is not currently active. We recommend referring to [GustaME](http://mb3is.megx.net/gustame), which provides useful documentation of many of the methods we were exploring here. 

microbiogeo
===========

Repository for code, data, and analysis results of the biogeographical
statistical methods comparison project (aka ```microbiogeo```).

To run the code in this repository, you will need the following dependencies
installed (versions tested against are in parentheses):

- [QIIME base install](http://qiime.org/) (1.8.0)
- [IPython](http://ipython.org/) (1.2.1)
- [pyzmq](http://zeromq.github.io/pyzmq/) (14.0.1)
- [nose](https://nose.readthedocs.org/en/latest/) (1.3.0)
- [R](http://www.r-project.org/) (3.0.2)
- R [optparse](http://cran.fhcrc.org/web/packages/optparse/index.html) package (1.0.2)
- R [vegan](http://cran.r-project.org/web/packages/vegan/index.html) package (2.0-10)
- R [ape](http://cran.r-project.org/web/packages/ape/index.html) package (3.0-11)

The easiest way to install the Python dependencies is via pip, e.g.:

    pip install numpy
    pip install qiime
    pip install ipython
    pip install pyzmq
    pip install nose

You will also need ```microbiogeo/code``` added to your ```PYTHONPATH``` and
```microbiogeo/code/scripts``` added to your ```PATH```.

To run the unit tests, assuming you are in the ```microbiogeo``` directory:

    nosetests code

To run the actual workflows, you will need to ```cd``` into the
```microbiogeo/code``` directory and start an IPython cluster with the number
of cores/processors you'd like parallel jobs to be executed on. For example,
the following command will start 4 IPython Engines:

    ipcluster start --n=4

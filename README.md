microbiogeo
===========
Repository for code, data, and analysis results of the biogeographical
statistical methods comparison project (aka ```microbiogeo```).

To run the code in this repository, you will need the latest development
versions of QIIME and IPython installed, with ```microbiogeo/code``` added to
your ```PYTHONPATH``` and ```microbiogeo/code/scripts``` added to your ```PATH```.

You can run all of the unit tests with the following command (assuming you are
in the ```microbiogeo``` directory):

    python code/tests/all_tests.py

To run the actual workflows, you will need to ```cd``` into the
```microbiogeo/code``` directory and start an IPython cluster with the number
of cores/processors you'd like parallel jobs to be executed on. For example,
the following command will start 4 IPython Engines:

    ipcluster start --n=4

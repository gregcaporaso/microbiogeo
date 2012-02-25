=======================
Simulation Code Example
=======================

Below is an short example of how to utilize the simulation code.

Open a python terminal and complete the following steps by typing them into the terminal: ::

	execfile('/Users/Damien/Documents/School/CS486C/qiimeutils/otu_simulation/flow_ms.py')

Where the file path points to the local file flow_ms.py.

:note: This may return the message:  "UserWarning: Not using MPI as mpi4py not found from cogent.util import parallel, terminal". This will not effect the functionality of the program.

::

	result = simulate_flow_with_rand_e_mat(16,4)

The first parameter represents OTUs (must be a power of 2) and the second is the environments (max of 26). So, in our example above we have 16 OTUs and 4 environments.

The function can also has some additional parameters that are set by default. These addition parameters are distribution_f, distro_mean, abundance_mode, ab_low, ab_high, and oneEnvAb. For this example we will not be overriding any of the defaulted parameters.

The result of this method call will be a tuple with two elements. The first element is a cogent tree object, and the second element is an abundance dictionary which gives the name of the node in the tree and the abundance of bacteria found at that node for each environment: ::

	tree = result[0]

Assign the first element in the list, the cogent tree object, to the variable tree: ::

	tree.root().Sequence

This list corresponds to how many microbes are in each environment: ::

	print tree.asciiArt()

From the results we can attempt to infer bacterial flow based off where the tips are found in the tree: ::

	tree.root().Q

The rates listed here correlate to the chance the microbe moves to a different environment.

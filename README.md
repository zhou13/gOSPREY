Introduction to gOSPREY
=======================
Welcome to the website for the gOSPREY!

gOSPREY is the abbreviation for <strong>O</strong>pen <strong>S</strong>ource
<strong>P</strong>rotein <strong>RE</strong>design for <strong>Y</strong>ou on a
<strong>G</strong>PU.  It enables the user to utilize a GPU to accelerate the
procedure of protein design in a massive parallel fashion.

gOSPREY is based on the software suite OSPREY, developed by the [Donald
Lab](http://www.cs.duke.edu/donaldlab/osprey.php) from Duke University.  Here is
a short introduction for OSPREY:

> OSPREY incorporates several different algorithmic modules for structure-based
> protein design, including a number of powerful Dead-End Elimination algorithms
> and the ensemble-based K\* algorithm for protein-ligand binding prediction.
> OSPREY allows the incorporation of continuous protein side-chain and
> continuous or discrete backbone flexibility, while maintaining provable
> guarantees with respect to the input model (input structure, rotamer library,
> energy function, and any backbone perturbations) for a given protein design
> problem. See full details of the different algorithmic modules in OSPREY.

Currently, the version of gOSPREY that you can download is based on OSPREY v2.1
beta with some improvements, including:

* A robust CUDA implementation of GA\*, a variant of A\* algorithm runs on any
  CUDA-compatible device in a massive parallel fashion.

* Improved computation of heuristic function for structure-based computation
  protein design, which is able to accelerate the protein design procedure
  significantly even without parallelization.

* An implementation of GSMA\* on CUDA, which is able to continue the design
  process even if the system is out of memory.

* An alpha release of OpenCL implementation of GA\*.  This is not packaged by
  default.  Advanced users with coding ability may want to try it in
  `native/additional` directory.

* A nice build/package system based on CMake so you should not worry about the
  ugly compilation process of CUDA (`nvcc`) and Java Native Interface!

Requirement
-----------

### Hardware Requirement
In order to use the function of GPU acceleration, the user must own a
CUDA-compatible video card, with the
[Compute-Capability](http://docs.nvidia.com/cuda/cuda-c-programming-guide/) at
least 1.2.

### Software Dependency
In order to run gOSPREY, a user must have the following installed:

*  A **Linux** operating system with **NVIDIA proprietary driver**

*  [NVIDIA CUDA SDK](http://developer.nvidia.com/cuda)

*  A JDK 1.7 implementation

*  A recent gcc/g++ release

*  CMake 2.8+

Installation
------------
1.  Check your environment setting.  Make sure you have CUDA installed on your
    Linux box.  And make sure the `deviceQuery` utility from
    `1_Utilities/deviceQuery` shipped by CUDA SDK returns normally.

2.  Download the source code from github.  If you have `git` installed, you can
    `cd` into your working directory and perform the clone operation:

        $ cd ~/src
        $ git clone https://github.com/zhou13/gOSPREY.git

    Alternatively, you can to download the tarball from the website and decompress
    it by yourself.

3.  Create a build directory for gOSPREY:

        $ cd gOSPREY
        $ mkdir build

4.  Run `cmake` to generate the Makefile.  Because gOSPREY uses Java Native
    Interface, a dynamic library must be installed under the
    `java.library.path`.  On most system, a prefix on `/usr` should do this job:

        $ cd build
        $ cmake -DCMAKE_INSTALL_PREFIX=/usr ..

5.  Compile gOSPREY and install its dynamic library:

        $ make
        $ sudo make install

Usage
-----
Thanks for the CMake and jar package system, the use of gOSPREY is pretty easy.
After the `make install` in the installation procedure, a file called
`osprey.jar` will be generated.  This contains all the Java classes needed by
gOSPREY.  You can copy/move this file to any place that make you feel
comfortable.

Let's use `ppi_GPU` as an example. You can find it under `doc/example/ppi_GPU/`.
Suppose you are still under the `build` directory.  Execute:

    $ cd ../doc/example/ppi_GPU
    $ java -jar ../../../build/osprey.jar -t 8 doDEE System.cfg DEE.cfg

Here, `-t 8` will make gOSPREY compute the energy matrix in 8 threads on your CPU.

Finally, hope that everything goes smooth for you!

Configuration
-------------
The document of original OSPREY can be found at `doc/manual.pdf`.  Besides that,
gOSPREY provided some additional parameter that a user need to configure.  You
can find an example under `doc/example/ppi_GPU/KStar.cfg`:

    enableAStarJava true
    enableAStarNativeC true
    enableAStarCUDA true
    maxNativeCPUMemory 5032706048
    maxNativeGPUMemory 5032706048
    numGPUWorkGroup 4
    numGPUWorkItem 192
    numGPUWorkItem2 192
    shrinkRatio 1

`enableAStarJava` determines whether the A\* module implemented by original
OSPREY will be enabled.  `enableAStarNativeC` determines whether the A\*
module implemented using native machine code through JNI with heuristic function
optimization will be used.  `enableAStarCUDA` determines whether GA\* will be
enabled through CUDA.  If more than one modules are enabled, gOSPREY will
compare the results returned by different modules to verily the correctness.

If CUDA is enable, `numGPUWorkGroup` \* `numGPUWorkItem` is the number of
parallel queues used in GA\*.  `numGPUWorkItem2` is the number of work items for
an individual work group when calculating the heuristic function in parallel.
When `shrinkRatio` is not equal to one, GSMA\* will be enabled.  In that case,
when the system runs out of memory, a fraction of nodes specified by `shrinkRatio`
will be dropped.  You may want to set it to `0.5`.

Support or Contact
------------------
Having trouble at installation or function? Feel free to contact the authors:
broken.zhou@gmail.com.

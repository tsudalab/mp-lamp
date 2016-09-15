# What is MP-LAMP?

MP-LAMP stands for Massive Parallel LAMP,
which is a parallel version of
[LAMP](http://a-terada.github.io/lamp/).
LAMP stands for Limitless-Arity Multiple-testing Procedure.


# Installation

## Installation to Amazon Web Service (AWS)

MP-LAMP will be ready by following the steps.

1. Create an Amazon EC2 using the Amazon Linux Image.
2. Download mp-lamp
3. Uncompress it.
4. Move to the top of the uncompressed directory.
5. Run the following command.

```text
$ bash aws/_installer_single.sh
```

## Installation to a local environment

Currently, Intel CPU and linux is assumed.
If you encounter troubles during the installation process,
please send us the error message and the environment.

### Prerequisite

| tools      | recommended version |
|------------|--------------------|
|Compiler    | g++ (4.3 or later) |
|MPI Library | [OpenMPI](https://www.open-mpi.org/), MPICH, MVAPICH or Intel MPI |
|build tool  | [SCons](http://scons.org/), 2.0.0 or higher (python is needed for SCons) |
|boost library | [boost library](http://www.boost.org/) 1.55.0 or later |
|gflags      | [gflags](http://gflags.github.io/gflags/) 2.0 or later |

Notes:
* For gcc, 4.9.3 or later is preferable.
Older gcc will produce slightly slower binary.

* The latest version of gflags requires CMAKE for the build tool.
For users not familiar with CMAKE,
we advise to use gflags v2.0 which could be installed by configure, make.
[gflags v2.0](https://github.com/gflags/gflags/archive/v2.0.tar.gz)

### Compilation

* Please satisfy the prerequisite.

* Copy local.sample.cfg to local.cfg and edit appropriately.
  * [compilers]
	* single: compiler for non-parallel code (g++ or icpc)
	* parallel: compiler for MPI (typically mpicxx)
	* options: additional options for compier (added to CXXFLAGS)
	* libs: additional options for library (added to LDFLAGS)
  * [paths]
	* include and library path
	* Not needed if there is not library in non-default location.

* Sample local.cfg

```text
[compilers]
single=g++
parallel=mpicxx
# an example for linux.
option=-DGTEST_USE_OWN_TR1_TUPLE=1 -DHAVE_CLOCK_GETTIME
libs=-lrt

[paths]
# include=/path/to/your/include_directory
# library=/path/to/your/include_library
```

* If local.cfg is ready, go to top directory of lamp_search and type

```text
$ scons
```

or for parallel build (for 4 threads)

```text
$ scons -j 4
```

It will take a while because it makes several types of binaries (opt, dbg, log)
for both sequential and parallel versions.
(In future releases, it will only make opt binary.)

* Parallel binary *mp_lamp* will be ready at

```text
lamp_search/mp_build/opt/mp-main/mp-lamp
```
	
## Usage

* *mp-lamp* could be used from command line.
For 32 processes,
$ mpiexec -hostfile ${machinefile} -np 32 ./mp-lamp --item item_file.csv --pos positive_file.csv --a 0.05 --show_progress --log
	* --item: item data file
	* -pos: positive data file
	* --a: significance level (default 0.05)
	* --show_progress:
	    It is adivsed to turn on show_progress for long jobs.
	* --log: Shows the breakdown of execution time. It is not needed for most users.
		It might be useful to find out problems when mp-lamp is unexpectedly slow.

## Sample Toy Data

* Item data file format.
By default, *mp-lamp* reads the following csv format item data.
It assumes that the first line includes the name of the items
and the rest of the lines have the name of the transactions at the beginning.

```text
#gene,TF1,TF2,TF3,TF4
A,1,1,1,0
B,1,1,1,0
C,1,0,0,1
D,0,0,0,0
E,1,1,1,0
F,1,0,0,0
G,1,1,1,1
H,0,0,0,0
I,0,1,0,1
J,0,0,1,0
K,0,0,0,1
L,0,0,0,1
M,0,0,0,1
N,1,1,1,0
O,0,0,0,0
```

* Positive data file format.
An example of the positive data format corresponding to the item data file
is shown below.
The first line is required to start with a "\#".
**Current version crashes if the number of lines does not match with the item file.**

```text
#gene,expression
A,1
B,1
C,0
D,0
E,1
F,0
G,1
H,1
I,0
J,0
K,0
L,0
M,1
N,1
O,0
```

## Sample usage and output

* Sample command and output of the 2-process parallel version solving the toy
  data. Do not forget to invoke the command using "mpiexec" or "mpirun".

```text
$ mpiexec -np 2 ./mp-lamp --item ./samples/sample_data/sample_item.csv --pos ./samples/sample_data/sample_expression_over1.csv --a 0.05 --show_progress
# item file    : ./samples/sample_data/sample_item.csv
# positive file: ./samples/sample_data/sample_expression_over1.csv
# # of transactions=          15	# of items=           4	# of total positives=           7	max freq=           7	max positive=           5	max items in trans.=           4
# preprocess end
# lambda=6	cs_thr[lambda]=               7	pmin_thr[lambda-1]=  0.00699301	num_expand=           1	elapsed_time=0.000616
# 1st phase start
# lambda=6	closed_set_num[n>=lambda]=           4	cs_thr[lambda]=               7	pmin_thr[lambda-1]=  0.00699301	num_expand=           1	elapsed_time=0.000661
# 1st phase end
# lambda=6	num_expand=           2	elapsed_time=0.001023
# 2nd phase start
# lambda=5	int_sig_lev=0.0125	elapsed_time=0.001052
# 2nd phase end
# closed_set_num=           5	sig_lev=0.01	num_expand=           3	elapsed_time=0.001165
# 3rd phase end
# sig_lev=0.01	elapsed_time=0.001564
# time all=    0.006031	time search=    0.001858
# min. sup=5	correction factor=5
# number of significant patterns=1
# pval (raw)    pval (corr)         freq     pos        # items items
0.00699301      0.034965               5       5	3	TF1	TF2	TF3
```

## Notes

* Current version does not work with "mpiexec -np 1".
  Please use at least two processes for the parallel version.

* Current version is only targeted for data with small number of transactions.
  For data with more than 100,000 transactions, please wait for the future updates.

## Contact

Please contact the following for bug reports, comments, or requests.

* yoshizoe(AT)acm.org

## License

MP-LAMP is an open source code project licensed under the Revised BSD license.


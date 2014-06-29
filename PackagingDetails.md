# Overview of the Package mechanisms for the nimble package

Basically, the package is now setup to be used as a regular R package, i.e.
via 
```shell
R CMD INSTALL nimble
```
from a local tar.gz file, or
```R
install.packages("nimble")
```
from some repository.

How did we organize this?

1. R code in the R directory,
2. a list of the R source files in the _DESCRIPTION_ file's *Collate* field
3. one C++ file in src/
4. configure script (generated from configure.ac) which determines the compilation 
   flags
  1. for include files in the Eigen distribution, 
  2. whether and how to compile the linking library libnimble, and
  3. if we are not using libnimble, we copy the code from inst/CppCode/ to src/ and compile that for the package DSO,
     (we can refine this in the future to only use the cpp files that the package needs. We can divide these cpp files into two groups - those needed by the package's .Call() invocations and those needed by the run-time generated DSOs)
  2. for how to compile the DSOs users of the nimble package generate at run-time, and
  4. how to find that library at run-time.
5. cleanup script to tidy up the src/ directory, .o files


The package uses symbol registration via the 
```R
useDynLib(nimble, .registration = TRUE)
```
directive in the *NAMESPACE* file.
The file *src/nimble.cpp* has a list of the routines that are (potentially)
called from R (via _.Call_()). 
The R package loading mechanism takes these, registers them to verify the
number of arguments, etc. and the registration arranges to create
R symbols/variables in the package's namespace corresponding to the 
registered names of these routines. We can then use this as, e.g.,
```R
.Call(getAvailableNames, a)
```
The value of _getAvailableNames_ is the resolved symbol in the nimble package's DSO.

When we use _.Call_() outside of the installed package mechanism, but
via *loadAllCode.R*, we have to deal with the different in the
registration.  We can define new variables corresponding to the
symbols that the R registration mechanism would create. The values of
these symbols could be either the resolve symbols in the DSO, or
alternatively, just the names of the routines.  The latter is
potentially dangerous, but probably sufficient for many uses.

(I am not certain what devtools does with registration.)


When R installs the package, it copies files and directories under the
*inst/* directory to the installation directory in which the package
will be available.  Therefore, R copies the *CppCode/* and *include*
directories to make them available at run-time for the package. These
are used when we create DSOs for the models nimble compiles.

(Note that using some of the facilities in devtools will cause these
directories not to be installed. That, and use in loadAllCode.R, is
why we need to be able specify where these directories are with
NimbleCodeDir and NimbleIncludeDir.)

In addition to the *CppCode/* and *include/* directories, we have a
*make/* directory under *inst/*. This contains two primary files -
*Makevars* and *Makevars_lib*.  When the configure script is run,
these are created from the corresponding .in files. These then contain
all of the necessary details for use in compiling the run-time DSOs.
When we compile the generated C++ code, we use R's `R CMD SHLIB`
command and that uses a Makevars directory to control the compilation
of the .cpp files and the linking of those files.  If we use
libnimble.so/dylib/dll, we use the *Makevars_lib* file.  Otherwise, we
use the *Makevars* file.  The *Makevars_lib* file takes care of being
able to find libnimble.so at run-time without having to set
`LD_LIBRARY_PATH` (or `DYLD_LIBRARY_PATH`).  This works for gcc/g++
and clang on Linux and OSX.  We'll have to do experiments on Windows,
but there we will be using gcc/g++.





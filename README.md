LMCLUS
======

LMCLUS is a linear manifold clustering algorithm, developed by Rave Harpaz. 
Current repository contains LMCLUS C++ library and R extensions.

Installation
============

LMCLUS is built via CMake (www.cmake.org). To build from the command line:

    $ cd lmclus
    $ cmake .
    $ make
    $ make test
  
If you want to create R extensions then 
    
    $ R CMD build lmclus
    $ R CMD INSTALL --build lmclus_1.0.0.tar.gz

Dependencies
------------

* gcc (>= 4.7.2)
* R (>= 2.15.0)
* BLAS (or MKL optional)

Licence
=======
Licensed under the LGPL, Version 3.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.gnu.org/copyleft/lesser.html

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

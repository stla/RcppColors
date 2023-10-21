## R CMD check results

OK


## Errors on Mac spotted by the CRAN results

The CRAN results show there are some C++ errors on Mac with the previous 
version. I added `#include <array>` in the header file and this has solved 
the problem (I checked with 'macbuilder').
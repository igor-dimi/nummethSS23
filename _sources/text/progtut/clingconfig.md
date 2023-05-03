
# Cing Config

To include the `hdnum` library in the jupyter notebook cling environment, following must be executed at the beginning of a notebook:

```cpp
#pragma cling add_include_path("/home/igor/Documents/uni/ss23/nummethSS23/hdnum/")
#include <iostream>
#include <cmath>
#include <complex>


#include "src/densematrix.hh"
#include "src/exceptions.hh"
#include "src/lr.hh"
#include "src/newton.hh"
#include "src/ode.hh"
#include "src/opcounter.hh"
#include "src/pde.hh"
#include "src/precision.hh"
#include "src/qr.hh"
#include "src/rungekutta.hh"
#include "src/sgrid.hh"
#include "src/timer.hh"
#include "src/vector.hh"
```
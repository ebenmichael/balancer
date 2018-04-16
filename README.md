# balancer
[![Project Status: WIP  Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)


## Overview
This package implements various forms of balancing weights which solve the general dual problem:

```
min f^*(X %*% theta) - X_t + h^*(theta)
```
This optimization problem is solved with accelerated proximal gradient descent.

## Installation
To install this package, first ensure that `devtools` is installed with

```
install.packages("devtools")
```

then install the package from GitHub with

```
devtools::install_github("ebenmichael/balancer")
```

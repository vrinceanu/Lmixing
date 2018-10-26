# Lmixing

__A Python module to calculate the rate coefficient for angular momentum mixing of Rydberg atoms in collisions with protons in astrophysical plasmas__

Distant collisions of protons change the angular momentum of highly excited hydrogenic Rydberg atoms with no energy transfer in a reaction: _p<sup>+</sup> + H(n, &#x2113;) &#x2192; p<sup>+</sup> + H(n,&#x2113;')_. This process is very efficient and it is responsible for creating the statistical equilibrium distribution of angular momentum states within the hydrogenic energy shell.

This module calculates the rate coefficients by using exact quantum mechanical results for low and moderate principal quantum numbers *n*, and by using a semiclassical approximation when the direct calculation is not practical.

In order to avoid truncation error and loss of precision due to subtraction of near-equal terms, the code exploits the unlimited integer Python type to calculate the Wigner's 6-j symbols, the *mpmath* package to extend the floating point precision in calculating high order polynomials.


## Installation

This module requires Pyhton 3.x and depends on [scipy](https://www.scipy.org/) for dealing with physical constants and unit conversion, and on [mpmath](http://mpmath.org) for extending the precision of floating point numbers beyond hardware capabilities.

You can download, clone and install the source from GitHub

```
$ git clone https://github.com/vrinceanu/Lmixing.git
$ cd Lmixing
$ python setup.py install --user
```

If you wish to know more in detail what git is and what you can do with it, the [github help page](https://help.github.com/articles/set-up-git) has all the references needed.

## Usage

As a quick example, we calculate the rate coefficient for

```Pyhton
>>> from Lmixing.rates import rate
>>> n = 10; L = 1; T = 10; ne= 100
>>> (rate(n,L,T,ne,method='semiclassical'), rate(n,L,T,ne,method='quantum'))
```

## Documentation



**Citing Lmixing:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1471776.svg)](https://doi.org/10.5281/zenodo.1471776)

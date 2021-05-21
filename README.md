# EigenApprox

Octave code to approximate kernel eigenfunctions in a Mercer expansion.
![](Figures/Cover_disk_contour.png)


## Background
Every strictly positive definite and continuous kernel <img src="https://render.githubusercontent.com/render/math?math=%24K%3A%5COmega%5Ctimes%5COmega%5Cto%5Cmathbb%7BR%7D%24"> on a bounded set <img src="https://render.githubusercontent.com/render/math?math=%24%5COmega%5Csubset%5Cmathbb%7BR%7D%24"> can be written as

<img src="https://render.githubusercontent.com/render/math?math=%24%0AK(x%2C%20y)%20%3D%20%5Csum_%7Bj%3D1%7D%5E%7B%5Cinfty%7D%20%5Clambda_j%20%5Cvarphi_j(x)%20%5Cvarphi_j(y)%2C%20%0A%24">

with 

<img src="https://render.githubusercontent.com/render/math?math=%5Clambda_j%20%5Cvarphi_j(x)%20%3D%20%5Cint_%7B%5COmega%7D%20%5Cvarphi_j(x)%20K(x%2C%20y)%20dy.%0A">

This code approximates the unknown eigenfunctions <img src="https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_1%2C%20%5Cdots%2C%20%5Cvarphi_n%2C%20%5Cdots%24"> and eigengvalues <img src="https://render.githubusercontent.com/render/math?math=%24%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_n%2C%20%5Cdots%24"> using the knowledge of <img src="https://render.githubusercontent.com/render/math?math=%24K%24"> and <img src="https://render.githubusercontent.com/render/math?math=%24%5COmega%24">.

The code implements the algorithm of
> G. Santin and R. Schaback, [Approximation of Eigenfunctions in Kernel-based Spaces](https://link.springer.com/article/10.1007/s10444-015-9449-5), Adv. Comput. Math., Vol. 42 (4), 973–993 (2016)
(see also the [preprint](https://arxiv.org/abs/1411.7656)).

## Quick start

You can start with one of the demos:
* [ApproximateEigenbasis1D.m](ApproximateEigenbasis1D.m): Example with the Mat&egrave;rn kernel on an interval.
* [ApproximateEigenbasis.m](ApproximateEigenbasis.m): Example with the Gaussian kernel on the unit disk.

Notice that both demos actually compute the eigenbasis elements as <img src="https://render.githubusercontent.com/render/math?math=%24%5Csqrt%7B%5Clambda_j%7D%20%5Cvarphi_j%24">, i.e., with a normalization that makes them orthonormal in the RKHS of the kernel.


## How to cite:
If you use this code in your work, please consider citing the paper

```bibtex:
@Article{Santin2016,
  author    = {Santin, Gabriele and Schaback, Robert},
  title     = {Approximation of eigenfunctions in kernel-based spaces},
  journal   = {Adv. Comput. Math.},
  year      = {2016},
  volume    = {42},
  number    = {4},
  pages     = {973--993},
  issn      = {1572-9044},
  doi       = {10.1007/s10444-015-9449-5},
}
```

Latex formulas are rendered using [https://jsfiddle.net/8ndx694g/](https://jsfiddle.net/8ndx694g/).

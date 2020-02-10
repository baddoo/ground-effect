# Exact ground effect

This repository contains the codes needed to produce the (theoretical) results in the paper ["Exact solutions for ground effect"](https://arxiv.org/abs/1912.02713). The basic setup is illustrated below:

![Schematic diagram](images/schematic.png?raw=true)

We solve the problem by conformally mapping an annulus to the physical domain, and solving the fluid problem inside the annulus. There are usually at least two quantities to compute when evaluating the solutions: the conformal map and the complex potential. Several choices for the conformal map are availabe within this code, as illustrated below:

![Conformal maps](images/conformal-maps.png?raw=true)

<p align="center"> 
<img src="images/conformal-maps.png?raw=true">
</p>

The repository includes the conformal maps:
* ```circularWing.m``` -- circular wing map,
* ```flatWng.m``` -- flat plate map,
* ```circularArcWing.m``` -- circular arc wing map, and
* ```centeredCircularArcWing.m``` -- centered circular arc wing map.

The repository includes the complex potential functions:
* ```circulation.m``` -- circulatory flows,
* ```vortex.m``` -- point vortices,
* ```uniform.m``` -- uniform flows,
* ```strain.m``` -- straining flows, and
* ```movement.m``` -- motion of the wing.

# Exact solutions for ground effect

This repository contains the codes used to produce the (theoretical) results in the paper ["Exact solutions for ground effect"](https://arxiv.org/abs/1912.02713). The basic setup is illustrated below:

<br/>
<p align="center"> 
<img src="images/schematic.png?raw=true" width="700px">
</p>
<br/>

We solve the problem by conformally mapping an annulus to the physical domain, and solving the fluid problem inside the annulus. There are usually at least two quantities to compute when evaluating the solutions: the conformal map and the complex potential. Several choices for the conformal map are availabe within this code, as illustrated below:

<br/>
<p align="center"> 
<img src="images/conformal-maps.png?raw=true" width="700px">
</p>
<br/>

The files correspond to different conformal maps as such:

File | Mapping
------------ | -------------
```circularWing.m``` | circular wing map
```flatWing.m``` | flat plate map
```circularArcWing.m```| circular arc wing map
```centeredCircularArcWing.m``` | centered circular arc wing map

A range of potential flow solutions are available, corresponding to various flow phenomena:

File | Flow phenomenon
------------ | -------------
```circulation.m``` | circulatory flow
 ```vortices.m``` | point vortices
```uniform.m```| uniform flow
```strain.m``` | straining flow
```movement.m```| motion of the wing

# BiophysicsTools
Tools for analyzing molecular and biophysical data

## Thermodynamical analysis of protein unfolding by circular dichroism

This code fits a protein unfolding curve recorded by circular dichroism (CD) to a two-state [Gibbs-Helmholtz equation](https://en.wikipedia.org/wiki/Gibbs–Helmholtz_equation) to extract thermodynamical information from it.

The fitting is performed as described by [Greenfield, N.J. (Nat Methods 2006)](https://www.nature.com/articles/nprot.2006.204)

The model relies on five parameters:
- T<sub>m</sub> = Melting temperature
- &Delta H = Enthalpy of unfolding
- C<sub>p</sub> = Heat capacity change of unfolding 
- &theta <sub>U</sub> = Ellipticity of the unfolded state
- &theta <sub>F</sub> = Ellipticity of the folded state

&theta; Hey!

The molar ellipticity &theta can be calculated and fitted to the following equations:

![equation0](https://latex.codecogs.com/gif.latex?\Delta&space;G&space;=&space;\Delta&space;H&space;\left&space;(&space;\frac{1-T}{T_m}&space;\right&space;)&space;-&space;\Delta&space;C_p&space;\left&space;[&space;\left&space;(&space;T_m&space;-&space;T&space;\right&space;)&space;&plus;&space;T\log{&space;\frac{T}{T_m}}&space;\right&space;])

![equation1](https://latex.codecogs.com/gif.latex?K&space;=&space;\exp{\left&space;(&space;\frac{-\Delta&space;G}{RT}&space;\right&space;)})

![equation2](https://latex.codecogs.com/gif.latex?\alpha&space;=&space;\frac{K}{1&plus;K})

![equation3](https://latex.codecogs.com/gif.latex?\theta&space;=&space;\alpha&space;\left&space;(&space;\theta_F&space;-&space;\theta_U&space;\right&space;)&space;&plus;&space;\theta_U)

The parameters are initially approximated by creating an instance of `CD_Unfolding` with some experimental data:

```python
Protein_Unfold = CD_Unfolding(T,Ellipticity,T_Units="C",Energy_Units="kcal")
```

![Before fitting](https://github.com/maximosanz/BiophysicsTools/blob/master/CircularDichroism_Unfolding/Before_Fitting.jpg)

And fitted by non-linear least-squares with the `Fit()` method:

```python
Protein_Unfold.Fit()
```

![After fitting](https://github.com/maximosanz/BiophysicsTools/blob/master/CircularDichroism_Unfolding/After_Fitting.jpg)

Key fitting information can be displayed by the `Info()` method:


```python
Protein_Unfold.Info()
```
```
Melting temperature =   44.71 C
Enthalpy            =   -44.1569 kcal / mol
deltaG at   25.00 C =   -12.3313 kcal / mol
```

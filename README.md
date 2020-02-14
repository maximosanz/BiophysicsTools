# BiophysicsTools
Tools for analyzing molecular and biophysical data

- [Circular dichroism of protein unfolding](#thermodynamical-analysis-of-protein-unfolding-by-circular-dichroism)
- [Molecular contacts from MD trajectories](#evaluating-3d-contacts-and-ring-current-effects-from-a-molecular-dynamics-trajectory)

## Thermodynamical analysis of protein unfolding by circular dichroism

#### CircularDichroism_Unfolding

This code fits a protein unfolding curve recorded by circular dichroism (CD) to a two-state [Gibbs-Helmholtz equation](https://en.wikipedia.org/wiki/Gibbs–Helmholtz_equation) to extract thermodynamical information from it.

The fitting is performed as described by [Greenfield, N.J. (Nat Methods 2006)](https://www.nature.com/articles/nprot.2006.204)

The model relies on five parameters:
- T<sub>m</sub> = Melting temperature
- &Delta;H = Enthalpy of unfolding
- C<sub>p</sub> = Heat capacity change of unfolding 
- &theta;<sub>U</sub> = Ellipticity of the unfolded state
- &theta;<sub>F</sub> = Ellipticity of the folded state


The molar ellipticity &theta; can be calculated and fitted to the following equations:

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

## Evaluating 3D contacts and ring-current effects from a molecular dynamics trajectory

#### MD_3DContacts

This code calculates atom- and residue-wise molecular contacts between two groups in an MD trajectory, returning contact matrices, optionally as a function of time.

Atom selection groups are specified in a selection file, which contains at two lines, describing each group separately, e.g.

```
protein and not name H*
name OW
```
to compute contacts between protein heavy atoms and water oxygens.

The option to calculate ring-current contacts on a group can be specified, in which case an atom selection is specified as the first line and multiple rings can be selected as further lines:

```
resid 50:253 and not name H*
resid 15 and resname TYR and name CG CD1 CD2 CE1 CE2 CZ
resid 23 and resname PHE and name CG CD1 CD2 CE1 CE2 CZ
resid 36 and resname TYR and name CG CD1 CD2 CE1 CE2 CZ
```
to compute ring-current contacts between protein heavy atoms from residue 50 onwards and three aromatic residues.

In addition to the contact cut-off distance ```epsilon```, ring-current effects take into account an angular cutoff. A vector ```v``` can be defined between the ring centroid and a target atom, and the angle &theta between ```v``` and the normal of the ring plane must be below the cut-off ```angle_epsilon```:

![Ring-current](https://github.com/maximosanz/BiophysicsTools/blob/master/MD_3DContacts/Ring_Current.png)

The code can be run in parallel using MPI, specifying the option -multi, e.g.

```
mpirun -np M python MD_3DContacts -s structFile -f trajFile -sel SelFile -cmat CMAT.npy -multi M
```
to split the job across ```M``` processors.

This code was used to analyse different molecular interactions (e.g. Hydrogen bonds and ring-current effects) between the small molecule EGCG and the enzyme mAcP in the following publication:

[Molecular determinants of the interaction of EGCG with ordered and disordered proteins. Fusco G., Sanz‐Hernandez M., Ruggeri F. S., Vendruscolo M., Dobson C. M., De Simone A. Biopolymers 2018](https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.23117)

![EGCG_Interactions](https://github.com/maximosanz/BiophysicsTools/blob/master/MD_3DContacts/EGCG_Interactions.png)

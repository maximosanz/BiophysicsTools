# BiophysicsTools
Tools for analyzing molecular and biophysical data

- [Circular dichroism of protein unfolding](#thermodynamical-analysis-of-protein-unfolding-by-circular-dichroism)
- [Generate elongated amyloid fibre structures](#making-large-amyloid-fibre-structures-from-a-pdb-file)
- [Principal component analysis of fluorescence emission data](#principal-component-analysis-of-fluorescence-spectra-to-fingerprint-protein-aggregates)
- [Molecular contacts from MD trajectories](#evaluating-3d-contacts-and-ring-current-effects-from-a-molecular-dynamics-trajectory)
- [Topological angle analysis in MD simulations](#computing-tilt-and-rotational-angles-of-&alpha;-helices-from-a-molecular-dynamics-trajectory)
- [Calculate oriented SSNMR observables from MD trajectories](#calculating-oriented-solid-state-nmr-observables-from-a-molecular-dynamics-trajectory)

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

An example jupyter notebook is provided [here](https://github.com/maximosanz/BiophysicsTools/blob/master/CircularDichroism_Unfolding/Thermodynamical%20analysis%20of%20protein%20unfolding%20by%20Circular%20Dichroism.ipynb).

## Making large amyloid fibre structures from a PDB file.

#### Fibre_maker

This scripts allows you to generate large fibrillar structures of a periodic amyloid assembly, incorporating a helical twist if specified (even if no cell unit is defined) or a curvature. Deposited PDB structures of amyloid fibrils generally represent a small cross-section of a long filament, which can be generated with this code:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fibre_maker/Twisted_fibre.png" width="800" title="Twisted_fibre">

Example input and output files corresponding to a [human prion fibre](https://www.nature.com/articles/s41594-020-0403-y) are provided. They can be reproduced using the command

```
python Fibre_maker.py -pdb 6uur.pdb -o 6uur_100fibre.pdb -unit A,B -repeated G,J -n_units 100 -pitch 126
```

Curved fibres can also be generated using the script, for example to reproduce an amyloid monolayer adhered to a curved hydrophobic surface:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fibre_maker/Curved_fibre.png" width="600" title="Curved_fibre">

An example curved fibre was generated using the following command:

```
python Fibre_maker.py -pdb 6uur.pdb -unit A -repeated G,J -n_units 40 -curvature 10 -o 6uur_curved.pdb
```

Curvature and helical pitch can be combined to create funky amyloid structures!

e.g. a protein [Möbius strip](https://en.wikipedia.org/wiki/M%C3%B6bius_strip).

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fibre_maker/mobius.gif" width="500" title="mobius">

## Principal component analysis of fluorescence spectra to fingerprint protein aggregates.

#### Fluorescence_PCA

There are a wide range of fluorescent molecules that have been described to bind to protein aggregates. Upon binding, their structural properties are modified leading to a change in their emission spectrum. A prominent example of such a dye is [Thioflavin-T](https://en.wikipedia.org/wiki/Thioflavin).

Some dyes show a different response upon binding to different protein aggregates, and therefore they can be exploited to specifically detect and fingerprint each aggregate type based on their modified emission spectra:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fluorescence_PCA/Figures/Emission.jpg" width="500" title="Emission">

A more powerful approach has recently emerged where the information from more than one fluorescent compound can be combined, providing a more specific signature for each protein aggregate. The combined analysis of different fluorescent dyes can be performed with principal component analysis (PCA).

This code performs a PCA on fluorescence data, using a ```.csv``` file as input. An example ```TEMPLATE.csv``` file is provided, with the following structure:

- The first column corresponds to the emission wavelength in nm.
- Each subsequent column contains the data for one Sample, with the sample name in the header row
- A new compound is introduced by a row with the keyword ```COMPOUND```, followed by the compound name

A ```Fluorescence_Data``` instance can be created from the ```.csv``` file:

```python
Data = Fluorescence_Data("TEMPLATE.csv")
```

And the PCA performed with the function ```do_PCA()```, which take some optional arguments:

```python
Fluorescence_Data.do_PCA(use_compounds=None,
               use_samples=None,
               exclude_buffer=True,
               N_Components=2,
               normalised=True,
               remove_columns=[]
```

- ```use_compounds```: a list with a subset of the compound names to use in the PCA. If ```None```, all compounds are used.
- ```use_samples```: a list with a subset of the sample names to use in the PCA. If ```None```, all samples are used.
- ```exclude_buffer```: Ignore samples labeled as ```buffer``` (case-insensitive) if present. Overriden by ```use_samples```.
- ```N_Components```: Number of principal components to calculate.
- ```normalised```: Whether to normalise the emission spectra to the ```[0,1]``` interval.
- ```remove_columns```: Ignore specific column indices from the PCA calculation.

The PCA results can be shown by plotting the first two principal components:

```python
Data.plot_PCA()
```

![PCA](https://github.com/maximosanz/BiophysicsTools/blob/master/Fluorescence_PCA/Figures/PCA.jpg)

The clustering of different proteins shows how the PCA analysis can capture structural differences in the aggregates.

The contribution of each compound to the PCA can also be shown:

```python
Data.plot_Compound_contribution()

```

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fluorescence_PCA/Figures/Compound_contribution_0.jpg" width="500" title="Compound_contribution_0">

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/Fluorescence_PCA/Figures/Compound_contribution_1.jpg" width="500" title="Compound_contribution_1">

An example jupyter notebook is provided [here](https://github.com/maximosanz/BiophysicsTools/blob/master/Fluorescence_PCA/Evaluate_PCA.ipynb).

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

![Ring-current](https://github.com/maximosanz/BiophysicsTools/blob/master/MD_3DContacts/Ring_Current_diagram.png)

The code can be run in parallel using MPI, specifying the option -multi, e.g.

```
mpirun -np M python MD_3DContacts -s structFile -f trajFile -sel SelFile -cmat CMAT.npy -multi M
```
to split the job across ```M``` processors.

This code was used to analyse different molecular interactions (e.g. Hydrogen bonds and ring-current effects) between the small molecule EGCG and the enzyme mAcP in the following publication:

[Molecular determinants of the interaction of EGCG with ordered and disordered proteins. Fusco G., Sanz‐Hernandez M., Ruggeri F. S., Vendruscolo M., Dobson C. M., De Simone A. Biopolymers 2018](https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.23117)

![EGCG_Interactions](https://github.com/maximosanz/BiophysicsTools/blob/master/MD_3DContacts/EGCG_Interactions.png)

## Computing tilt and rotational angles of &alpha;-helices from a molecular dynamics trajectory

#### MD_Angle_Analysis

This scripts can calculate the distribution of tilt and rotational angles for &alpha;-helices during a molecular dynamics simulation. The helical axis of a protein segment is defined by performing a least-squares fit on the 3D coordinates of its backbone atoms. The angles are then defined with respect to an external fixed vector **B<sub>0</sub>**, analogous to the external magnetic field in a nuclear magnetic resonance experiment.
  
The following figure from [De Simone et al., Biophys J. 2014](https://www.cell.com/biophysj/fulltext/S0006-3495(14)00326-9) shows how these angles are defined in a molecular system:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/MD_Angle_Analysis/Angles_definition.png" width="550" title="Angle_definition">

In addition, the code can be employed in oligomeric proteins, where more than one copy of the same molecule is present. A global angle of the oligomeric assemble can be defined too, allowing the computation of individual angles relative to the global vector.

This code was employed in the following publication, where a more detailed description of the method can be found:

[Accurate Determination of Conformational Transitions in Oligomeric Membrane Proteins. Sanz‐Hernandez M. et al. Scientific Reports 2016](https://www.nature.com/articles/srep23063)

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/MD_Angle_Analysis/PLN_Angles.png" width="800" title="PLN_Angles">

## Calculating oriented solid-state NMR observables from a molecular dynamics trajectory

#### MD_SSNMR_Calculator

Oriented solid-state NMR (oSSNMR) is a method that provides topological and dynamical information about the orientations of proteins in the membrane. Two key oSSNMR observables are the [chemical shift anisotropy](https://en.wikipedia.org/wiki/Solid-state_nuclear_magnetic_resonance#Examples_of_anisotropic_nuclear_interactions) (CSA) and [dipolar coupling](https://en.wikipedia.org/wiki/Magnetic_dipole%E2%80%93dipole_interaction) (DC). 

In the protein backbone, the CSA and DC can be calculated given the 3D coordinates of atoms, according to the following structural models:

<img src="https://github.com/maximosanz/BiophysicsTools/blob/master/MD_SSNMR_Calculator/SSNMR_models.png" width="800" title="SSNMR_Models">

The CSA is modelled as a rank-2 tensor centered on the amide nitrogen atom, with the following equation:

![equation0](https://quicklatex.com/cache3/8c/ql_7486fa6c36a70944847264b353f9418c_l3.png)

where &delta;<sub>11</sub>, &delta;<sub>22</sub> and &delta;<sub>33</sub> are respectively set to 64.0, 76.0, 216.9 ppm for non-glycine residues and 46.5, 66.3, 211.6 ppm for glycine. &alpha; and &beta; are the Euler angles (in degrees) used to transform from the laboratory frame to the principal axis frame.

The DC is only dependent on the length of the covalent bond and its angle &theta; with respect to the external magnetic
field B<sub>0</sub>:

![equation1](https://quicklatex.com/cache3/80/ql_e0124faadb679341e4a8ed338b3b2680_l3.png)

This script calculates the agreement between these experimental oSSNMR observables and MD trajectories in the form of Q-factors, as described in the following article.

[Accurate Determination of Conformational Transitions in Oligomeric Membrane Proteins. Sanz‐Hernandez M. et al. Scientific Reports 2016](https://www.nature.com/articles/srep23063)

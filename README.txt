Quantum Hall Edge Program

This program is for calculating the edge spectra for quantum Hall systems at
Laughlin 1/m fractional filling for fermions.

The relevant modules are:
  - IQHEDiag: Used for edge calculations in the integer case.
  - FQHEDiag: Used for edge calculations in the fractional case.
  - NBodyBasisMatrixElementCalc: for calculating matrix elements between two slater determinent states.
  - CoulombMatrixFunctions: for calculating coulomb matrix elements bwetween two particle states.
  - haldanePotentials: for calculating matrix elements of two particle states when the potentail is expressed as Haldane pseudo-potentials.
  - potentials: stoes some potentials as Haldane pseudo-potential parameters.
  - DiagHamInterface: used to slater decompose jack polynomial states using the DiagHam program.
  - waveFunctionClasses: module that defines a wavefunction class used in the frational case.
  - usefulTools: a collection of helpful functions used by more than one module.

Requirements:
  - This needs to be ran on Mac or Linux.
  - mpmath
  - sympy
  - matplotlib
  - DiagHam installation in the folder above this program.

Spectra.p is a pickle file for storing spectra of coulomb edge in integer case.
These are indexed by [(N, LMax)].

coulombMatrixElements.p stores matrix elements of the coulomb interaction between two particles as a pickle file.

All other undocumented modules are for quick testing in this project and so are not
documented here.

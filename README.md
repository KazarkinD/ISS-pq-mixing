# ISS-pq-mixing
This package implements the Inverse Seesaw model ISS(p,q) with parameters p (number of right-handed neutrinos) and q (number of singlets) in Wolfram Mathematica. It provides tools for numerical calculations of mass spectra, mixing matrices, and Yukawa parameter generation using the Casas-Ibarra parametrization.

## Features

- Construction of general ISS(p,q) mass matrices
- Numerical Takagi decomposition for symmetric matrices
- Casas-Ibarra parametrization enforcing observed neutrino masses and PMNS mixing
- Random Yukawa matrix generation (uniform and log-uniform distributions)
- Calculation of HNL mixing parameters $U_I^2 = \sum_\alpha |U_{\alpha I}|^2$
- Flavor-specific mixing parameters $U_\alpha^2 = \sum_I |U_{\alpha I}|^2$
- Mass spectrum visualization with logarithmic scales
- Customizable energy scales and PMNS parameters

## Requirements
- Wolfram Mathematica 12.0 or higher

## Installation
Clone the repository:

```bash
git clone https://github.com/KazarkinD/ISS-pq-mixing/
```

In your Mathematica notebook, load the package:

```mathematica
SetDirectory[NotebookDirectory[]];
Get["iss_functions.m"];
Get["iss_param.m"];
```

## Usage

Minimal example:

```mathematica
(* Generate random ISS(2,3) mass matrix with normal hierarchy *)
Omega = {{1, 0}, {0, 1}};
M = randomMassMatrixISSpq["NH", Omega, 0, 2, 3, issScales, "log-uniform"];

(* Compute mass spectrum and mixing matrix *)
{masses, Umix} = takagiDecomposition[M];

(* HNL mixing for lightest sterile state (dark matter candidate) *)
UUdm = getUUmixing[M, 50, 1];
(* Flavor-specific mixings *)
UUe = getUUalphaMixing[M, 50, 1];
UUmu = getUUalphaMixing[M, 50, 2];
UUtau = getUUalphaMixing[M, 50, 3];
```

## Functions

### Core Functions

| Function                                     | Description                                                     |
| -------------------------------------------- | --------------------------------------------------------------- |
| `massMatrixISSpq[p, q, scales, YD, YR, Ymu]` | Construct ISS(p,q) mass matrix                                  | 
| `takagiDecomposition[matrix, prec]`          | $\textrm{Takagi decomposition } M=U D U^T$                      |
| `getUUmixing[matrix, prec, Ihnl]`            | $\textrm{HNL mixing } U_I^2$​                                    |
| `getUUalphaMixing[matrix, prec, alpha]`      | $\textrm{Flavor mixing } U_\alpha^2$                            |
| `getMassMatrixBlocks[p, q, massMatrix]`      | Extract m_D,M_R,μ blocks                                        |
| `getThetaMatrix[p, q, massMatrix]`           | $\textrm{Compute } \tilde{\theta} \textrm{ mixing matrix}$      |
| `getMixMatrix[theta, Ua, Us]`                | Construct full mixing matrix                                    |

### Casas-Ibarra Parametrization

| Function                                                                      | Description                          |
| ----------------------------------------------------------------------------- | ------------------------------------ |
| `yukawaDiracCasasIbarra[hierarchy, mlight, p, q, scales, Omega, YR, Ymu]`     | Dirac Yukawa from CI parametrization |
| `massMatrixISSpqCasasIbarra[hierarchy, mlight, p, q, scales, Omega, YR, Ymu]` | Mass matrix with CI parametrization  |

### Random Generation

| Function                                                                            | Description                             |
| ----------------------------------------------------------------------------------- | --------------------------------------- |
| `randomYDirac[nrows, ncols, regime, zmin, zmax]`                                    | Random Dirac Yukawa matrix              |
| `randomYMajorana[nrows, regime, zmin, zmax]`                                        | Random symmetric Majorana Yukawa matrix |
| `randomMassMatrixISSpq[hierarchy, Omega, mlight, p, q, scales, regime, zmin, zmax]` | Random ISS mass matrix                  |

### Utilities

| Function                                                                             | Description                 |
| ------------------------------------------------------------------------------------ | --------------------------- |
| `massGenData[hierarchy, Omega, mlight, p, q, scales, regime, zmin, zmax, ngen]`      | Generate mass spectra array |
| `plotRandomMasses[hierarchy, Omega, mlight, p, q, scales, regime, zmin, zmax, ngen]` | Plot mass spectrum          |

## Parameters

Default parameters are defined in `iss_param.m`:

| Parameter       | Default                                                   | Description                                       |
| --------------- | ----------------------------------------------------------| ------------------------------------------------- |
| `issP`          | 2                                                         | $\textrm{Number of } \nu_R \textrm{ fields}$      |
| `issQ`          | 3                                                         | $\textrm{Number of} S_R \textrm{fields}$          |
| `mDscale`       | 172.2 GeV                                                 | Active neutrino Dirac mass scale                  |
| `mRscale`       | $10^6$ GeV                                                | Right-handed neutrino Dirac mass scale            |
| `muScale`       | $10^{−5}$ GeV                                             | Majorana mass scale                               |
| `zMIN`          | $10^{−6}+i\cdot10^{−6}$                                   | Minimum Yukawa value                              |
| `zMAX`          | $1+i$                                                     | Maximum Yukawa value                              |
| `massSquareDif` | $(7.5\times10^{−23},2.4\times10^{−21})\textrm{ GeV}^2$ | $\textrm{Solar and atmospheric } \Delta m^2$       |

PMNS parameters follow NuFIT best-fit values.

## Citation

If you use this package in your research, please cite:
```latex
M.~N.~Dubinin, E.~Y.~Fedotova and D.~M.~Kazarkin,
``Active-sterile neutrino mixing in the ISS(p,q) inverse seesaw models,''
[arXiv:2603.12987 [hep-ph]].
```

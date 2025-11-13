# Symbolic_Differential_Geometry_Library

Symbolic Python library for differential geometry and relativity: custom vectors/matrices, parametrized surfaces, general metrics, curvature, Gauss maps, fundamental forms, Christoffel, Ricci, Einstein and energy–momentum tensors, Kulkarni–Nomizu, Weyl/Schouten/Cotton tools, built on SymPy.

## Features

- **Vector & matrix types**
  - Custom `Vector` and `Matrix` classes.
  - Dot and cross products, norms, traces, determinants, inverses.
  - Basic linear algebra utilities (row reduction, augmented matrices, identity construction, etc.).

- **Parametrized surfaces in ℝⁿ**
  - First, second, and third fundamental forms: `firstFF`, `secondFF`, `thirdFF`.
  - Gauss map and its differential: `gauss_map`.
  - Gaussian curvature: `gaussianCurvature`.
  - Christoffel symbols on surfaces: `Christoffel_1`, `Christoffel_2`.
  - Riemann and Ricci curvature tensors from parametrizations: `Riemann_Curvature`, `Ricci_Curvature`, `Ricci_Scalar`.

- **Metric-based tensor calculus**
  - Christoffel symbols from a metric: `Christoffel_1_metric`, `Christoffel_2_metric`.
  - Riemann curvature tensor from a metric: `Riemann_Curvature_metric`.
  - Ricci tensor and scalar curvature from a metric: `Ricci_Curvature_metric`, `Ricci_Scalar_metric`.
  - Minkowski metric with signature (-, +, …, +): `Minkowski(n)`.

- **General relativity helpers**
  - Einstein tensor from a parametrized spacetime or metric: `Einstein_Tensor`, `Einstein_Tensor_metric`.
  - Energy–momentum tensor via Einstein’s equation: `EMT`, `EMT_metric`.
  - Kulkarni–Nomizu product for (0,2)-tensors: `Kulkarni_Nomizu`.
  - Weyl tensor, Schouten tensor, and Cotton tensor for both parametrized surfaces and metrics:
    - `Weyl_Tensor`, `Weyl_Tensor_metric`
    - `Schouten_Tensor`, `Schouten_Tensor_metric`
    - `Cotton_Tensor`, `Cotton_Tensor_metric`.

- **Constants**
  - `LAMBDA`: cosmological constant.
  - `KAPPA`: 8πG/c⁴ coupling constant.

All computations are symbolic and built on top of SymPy for clarity in learning and experimentation.

## Installation

This library is a single file, `Differential_Geometry.py`.

1. Download or copy `Differential_Geometry.py` into your project.
2. Make sure it is on your Python path.
3. Import what you need:

```python
import sympy as sym
from Differential_Geometry import Vector, Matrix
from Differential_Geometry import firstFF, gaussianCurvature
from Differential_Geometry import Minkowski, Einstein_Tensor_metric
```

## Quick start

### 1. Basic vectors and matrices

```python
import sympy as sym
from Differential_Geometry import Vector, Matrix, dot, cross, mag

v = Vector([1, 2, 3])
w = Vector([sym.Symbol('a'), 0, 1])

d = dot(v, w)      # symbolic dot product
c = cross(v, w)    # cross product in R^3
m = mag(w)         # symbolic magnitude
```

### 2. Parametrized surface and Gaussian curvature

```python
import sympy as sym
from Differential_Geometry import Vector, firstFF, secondFF, gaussianCurvature, gauss_map

u, v = sym.symbols('u v', real=True)

# Unit sphere parametrization in R^3
f = Vector([
    sym.cos(u) * sym.cos(v),
    sym.cos(u) * sym.sin(v),
    sym.sin(u)
])

g = firstFF(f, [u, v])             # first fundamental form (metric on the sphere)
h = secondFF(f, [u, v])            # second fundamental form
K = gaussianCurvature(f, [u, v])   # Gaussian curvature (should simplify to 1)
N = gauss_map(f, [u, v])           # Gauss map
```

### 3. Curvature from a metric

```python
import sympy as sym
from Differential_Geometry import Matrix
from Differential_Geometry import (
    Christoffel_2_metric,
    Riemann_Curvature_metric,
    Ricci_Curvature_metric,
    Ricci_Scalar_metric
)

x, y = sym.symbols('x y', real=True)

# Simple 2D metric: conformal factor e^{2x} on R^2
g = Matrix(2, 2, [
    sym.exp(2*x), 0,
    0, sym.exp(2*x)
])

Gamma = Christoffel_2_metric(g, [x, y])
Riemann = Riemann_Curvature_metric(g, [x, y])
Ric = Ricci_Curvature_metric(g, [x, y])
R = Ricci_Scalar_metric(g, [x, y])
```

### 4. Einstein tensor and energy–momentum tensor

```python
import sympy as sym
from Differential_Geometry import Matrix, Minkowski
from Differential_Geometry import Einstein_Tensor_metric, EMT_metric

# 3+1 dimensional Minkowski spacetime
eta = Minkowski(3)  # returns a 4x4 metric with signature (-, +, +, +)

t, x, y, z = sym.symbols('t x y z', real=True)
coords = [t, x, y, z]

Ein = Einstein_Tensor_metric(eta, coords)   # should be (symbolically) zero up to cosmological constant
T = EMT_metric(eta, coords)                 # corresponding energy–momentum tensor
```

## Notes and limitations

- Expressions can grow large; use SymPy’s simplification tools (`sym.simplify`, `sym.trigsimp`, etc.) as needed.
- Many routines assume that:
  - Vectors are instances of this module’s `Vector` class.
  - Metrics are instances of this module’s `Matrix` class.
- Sign conventions follow a mostly GR-style (-, +, …, +) signature but you should double-check for your specific application.

## Contributing

This library is primarily designed for learning and experimentation. If you extend or optimize it, consider documenting new functions in this README and adding small usage examples.

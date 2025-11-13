# Symbolic_Differential_Geometry_Library
A Python library for exploring differential geometry symbolically. Includes custom vectors and matrices, tools for parametrized surfaces, curvature, Gauss maps, fundamental forms, Christoffel symbols, and even Einstein tensors. Built on SymPy for clarity in learning and experimentation.

This project is designed for learning, experimentation, and research in differential geometry, general relativity, and geometric computation.

## Features

### Core Classes
- **Vector** class  
  Supports dot products, cross products, differentiation, magnitude, addition and subtraction, and compatibility checks.
- **Matrix** class  
  Supports determinant, trace, transpose, inverse, matrix multiplication, vector multiplication, augmented matrices, and row-reduction.

### Geometry Tools
- First, Second, and Third Fundamental Forms  
- Gauss Map  
- Gaussian Curvature  
- Differential Map  
- Christoffel Symbols (first and second kind)  
- Covariant Derivative  
- Riemann Curvature Tensor  
- Ricci Tensor and Ricci Scalar  

### Relativity Tools
- Minkowski metric generator  
- Einstein Tensor  
- Energy-Momentum Tensor  

## Example Usage

### Create a Vector
```python
from Differential_Geometry import Vector

v = Vector([1, 2, 3])
print(v)
```

### Create a 2×2 Matrix and compute determinant
```python
from Differential_Geometry import Matrix

A = Matrix(2, 2, [1, 2,
                  3, 4])
print(A.det())
```

### Gaussian Curvature of a Parametrized Surface
```python
import sympy as sym
from Differential_Geometry import Vector, gaussianCurvature

u, v = sym.symbols('u v')

# Parametric sphere (radius 1)
x = sym.cos(u) * sym.sin(v)
y = sym.sin(u) * sym.sin(v)
z = sym.cos(v)

f = Vector([x, y, z])

K = gaussianCurvature(f, [u, v])
print(K)
```

## Requirements

- Python 3.x  
- SymPy

Install SymPy:
```bash
pip install sympy
```

## Why This Library?

This project allows users to compute geometric quantities directly from their definitions without relying on black-box packages. It is useful for students, instructors, researchers, and anyone studying general relativity or surface theory.

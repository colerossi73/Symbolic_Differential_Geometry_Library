import sympy as sym
from Differential_Geometry import *

u = sym.symbols('u')
v = sym.symbols('v')
w = sym.symbols('w')
r = sym.symbols('r')
s = sym.symbols('s')
t = sym.symbols('t')


var_list = [u, v, w]

f1 = Vector([u*sym.cos(v), u*sym.sin(v), u*sym.cos(w), u*sym.sin(w)])
nu = gauss_map(f1, var_list)
print(nu)
print("")
g = firstFF(f1, var_list)
h = secondFF(f1, var_list)
printMatrix(g)
print("")
printMatrix(h)
print("")
Gamma = Christoffel_1(f1, var_list)

for k in range(0, 3):
    for i in range(0, 3):
        for j in range(0, 3):
            print("Gamma[{}][{}][{}] = {}".format(k, i, j, Gamma[k][i][j]))



#f2 = Vector([sym.cos(u)*sym.cos(v)*sym.cos(w), sym.sin(u)*sym.cos(v)*sym.cos(w), sym.sin(v)*sym.cos(w), sym.sin(w)])
#nu = gauss_map(f2, var_list)
#print(nu)
#g = firstFF(f2, var_list)
#h = secondFF(f2, var_list)
#printMatrix(g)
#printMatrix(h)

#K = gaussianCurvature(f, var_list)
#print(K)

#Df = differential(f, var_list)
#printMatrix(Df)

#g = firstFF(f, var_list)
#g_inv = inverse2(g)
#printMatrix(g)
#printMatrix(g_inv)
#Gamma = Christoffel_1(f, var_list)
#Gamma2 = Christoffel_2(f, var_list)
#print(Gamma2[0][0][0])
#print(Gamma2[0][0][1])
#print(Gamma2[0][1][0])
#print(Gamma2[0][1][1])
#print(Gamma2[1][0][0])
#print(Gamma2[1][0][1])
#print(Gamma2[1][1][0])
#print(Gamma2[1][1][1])



#M = Matrix(3, 3, [u, 0, 0, 0, v, 0, 0, 0, u*v])
#printMatrix(vm.RRE_AM(M))


#X = Vector([u, v])
#Y = Vector([2*u, 3*v])
#DelxY = Covariant_Deriv(f, var_list, X, Y)
#print(DelxY)

#var_list = [u, v]
#R = sym.symbols('R')
#f = Vector([R * sym.cos(u)*sym.sin(v), R * sym.sin(u)*sym.sin(v), R * sym.cos(v)])
#g = firstFF(f, var_list)
#printMatrix(g)
#Christoffel = Christoffel_2(f, var_list)
#Riemann = Riemann_Curvature(f, var_list)
#e = thirdFF(f, var_list)
#printMatrix(e)

#Ric = Ricci_Curvature(f, var_list)
#printMatrix(Ric)

#R_scalar = Ricci_Scalar(f, var_list)
#print(R_scalar)

##for i in range(0, 2):
##    for j in range(0, 2):
##        for k in range(0, 2):
##            print("Christoffel[{}][{}][{}] = {}".format(i, j, k, Christoffel[i][j][k]))

##for i in range(0, 2):
##    for j in range(0, 2):
##        for k in range(0, 2):
##            for l in range(0, 2):
##                print("R[{}][{}][{}][{}] = {}".format(i, j, k, l, Riemann[i][j][k][l]))


###
# Möbius Strip with n half-twists
###

#r = 1
#n = 1
##x = (r + v * sym.cos(n*u/2)) * sym.cos(u)
##y = (r + v * sym.cos(n*u/2)) * sym.sin(u)
##z = v * sym.sin(n*u/2)
##
##f = Vector([x, y, z])
##print("f")
##print(f)
##f1 = partial(f, u)
##print("f_u")
##print(f1)
##f2 = partial(f, v)
##print("f_v")
##print(f2)
##nu = gauss_map(f, u, v)
##print("nu")
##print(nu)
##g = firstFF(f, u, v)
##printMatrix(g)


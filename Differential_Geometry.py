###
# Program for Vectors, Matrices, and Parametrized Surfaces
# with symbolic manipulations. Note: Results will sometimes appear more complicated
# than they should be due to Python not knowing how to simplify square roots.
#
# @Author Cole Rossi
###

import sympy as sym
import re

###
# Defining Constants
###

LAMBDA = 1.1056e-52                 # m^-2          Cosmological Constant
C = 2.99792458e8                    # m/s           Speed of Light
G = 6.6743015e-11                   # N m^2 kg^-2   Gravitational Constant
KAPPA = (8 * sym.pi * G) / C**4     # Einstein's Gravitational Constant

###
# Vector and Matrix Classes
###

class Vector:

    def __init__(self, list_elem, row_vec = False):     #Default is a column vector
        self.components = tuple(list_elem)
        self.dim = len(list_elem)
        self.isrow = row_vec

    def __repr__(self):
        return str(self.components)


class Matrix:

    def __init__(self, n, m, list_matrix_inputs):
        if len(list_matrix_inputs) != n*m:
            print("Error: Must have {} inputs to initialize a {} x {} matrix.".format(n*m, n, m))

        if n == m:
            self.isSquare = True
        else:
            self.isSquare = False

        matrix_init = []
        self.rows = n
        self.cols = m
        k = 0
        for j in range(0, n):
            row_init = []
            for l in range(0, m):
                row_init.append(list_matrix_inputs[k])
                k += 1
            matrix_init.append(tuple(row_init))
        self.matrix = tuple(matrix_init)

    def __repr__(self):
        return str(self.matrix)

    def trace(self):
        if not self.isSquare:
            return -1
        tr = 0
        for i in range(0, self.rows):
            tr += self.matrix[i][i]
        return tr

    def det(self):
        if not self.isSquare:
            print("Matrix must be square to have a determinant.")
            return -1

        n = self.rows
        if self.rows == 2:
            determinant = self.matrix[0][0] * self.matrix[1][1] - self.matrix[0][1] * self.matrix[1][0]
            return sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(determinant))))
        else:
            first_row = []
            for i in range(0, n):
                first_row.append(self.matrix[0][i])

            determinant = 0

            for i in range(0, len(first_row)):
                sub_mat_ele = []
                for j in range(1, n):
                    for k in range(0, n):
                        if k == i:
                            continue
                        else:
                            sub_mat_ele.append(self.matrix[j][k])
                M = Matrix(n-1, n-1, sub_mat_ele)
                determinant += ((-1) ** i) * M.det() * first_row[i]
            return sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(determinant))))

    def transpose(self):
        if self.isSquare:
            matrix_init = []
            ele = []

            for i in range(0, self.rows):
                for j in range(0, self.cols):
                    ele.append(self.matrix[j][i])

            k = 0
            for j in range(0, self.cols):
                row_init = []
                for l in range(0, self.rows):
                    row_init.append(ele[k])
                    k += 1
                matrix_init.append(tuple(row_init))
            self.matrix = tuple(matrix_init)
        else:
            mat = list(self.matrix)
            t_mat = zip(*mat)
            self.matrix = tuple(t_mat)
            temp = self.rows
            self.rows = self.cols
            self.cols = temp



###
# Functions
###

# Create nxn Identity Matrix
def identity(n):
    ele = []
    temp = 0
    for i in range(0, n):
        for j in range(0, n):
            if j == temp:
                ele.append(1)
            else:
                ele.append(0)
        temp += 1

    Id = Matrix(n, n, ele)
    return Id

# Take dot product of two vectors
def dot(vec1, vec2):
    if (vec1.dim != vec2.dim):
        print("Vectors do not have same dimension.")
        return -1
    tot = 0
    for i in range(0, vec1.dim):
        tot += sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(vec1.components[i] * vec2.components[i]))))
    return sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(tot))))

# Take cross product of n-1 vectors in R^n
def cross(vec_arr):
    if len(vec_arr) == 1:
        print("Must have at least two vectors.")
        return Vector([-1])
    elif len(vec_arr) == 2:
        vec1 = vec_arr[0]
        vec2 = vec_arr[1]
        if vec1.dim != 3 or vec2.dim != 3:
            return Vector([-1])

        x = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(vec1.components[1] * vec2.components[2] - vec1.components[2] * vec2.components[1]))))
        y = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(vec1.components[2] * vec2.components[0] - vec1.components[0] * vec2.components[2]))))
        z = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(vec1.components[0] * vec2.components[1] - vec1.components[1] * vec2.components[0]))))
        cross = Vector([x, y, z])
        return cross
    else:
        dimension = vec_arr[0].dim
        if dimension == 3:
            print("You only need two vectors for a cross product.")
            return Vector([-1])

        for i in range(0, len(vec_arr)):
            if vec_arr[i].dim == dimension:
                pass
            else:
                print("Vectors must be of the same dimension.")
                return Vector([-1])

        basis = [sym.Symbol('b{}'.format(i), commutative=False) for i in range(0, dimension)]
        basisMatrix = Matrix(dimension, 1, basis)

        AM = augmented_matrix(basisMatrix, vec_arr[0])
        for i in range(1, len(vec_arr)):
            AM = augmented_matrix(AM, vec_arr[i])

        AM.transpose()
        determinant = AM.det()
        stringForm = str(sym.expand(determinant))
        string = stringForm.replace(" ", "")

        ele = []
        for i in range(0, dimension):
            ele.append("0")
        j = 0
        s = ""
        while j < len(string):
            if string[j] == "b" and s == "":
                s = s + "1"
                ele[int(string[j+1])] = s
                j += 2
                s = ""
            elif string[j] == "b" and s == "-":
                s = s + "1"
                ele[int(string[j+1])] = s
                j += 2
                s = ""
            elif string[j] == "b" and s == "+":
                s = s + "1"
                ele[int(string[j+1])] = s
                j += 2
                s = ""
            elif string[j] == "b" and s != "" and s[-1] != "*":
                ele[int(string[j+1])] = s
                j += 2
                s = ""
            elif string[j] == "b" and s != "" and s[-1] == "*":
                s = s[:-1]
                ele[int(string[j+1])] = s
                j += 2
                s = ""
            else:
                s = s + string[j]
                j += 1
        for i in range(0, len(ele)):
            ele[i] = sym.parse_expr(ele[i])
        cross = Vector(ele)
        return cross

# Take partial derivative of vector
def partial(vec, var):
    ele = []
    for i in range(0, vec.dim):
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(sym.diff(vec.components[i], var))))))

    partial = Vector(ele)
    return partial

# Calculate Magnitude of vector
def mag(vec):
    if not isinstance(vec, Vector):
        print("Input must be a vector.")
        return -1
    magnitude = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(sym.sqrt(dot(vec, vec))))))
    return magnitude

# Calculate Gauss Map for Parametrized Surface
def gauss_map(f, var_list):
    f_i = []
    n = len(var_list)

    for i in range(0, n):
        f_i.append(partial(f, var_list[i]))

    cprod = cross(f_i)
    magnitude = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(mag(cprod)))))
    ele = []
    for i in range(0, cprod.dim):
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(cprod.components[i] / magnitude)))))
    gauss = Vector(ele)
    return gauss

# Calculate First Fundamental Form for Parametrized Surface
def firstFF(f, var_list):
    g_elem = []
    f_i = []
    n = len(var_list)

    for i in range(0, n):
        f_i.append(partial(f, var_list[i]))

    for i in range(0, n):
        for j in range(0, n):
            g_elem.append(dot(f_i[i], f_i[j]))

    g = Matrix(n, n, g_elem)
    return g

# Calculate Second Fundamental Form for Parametrized Surface
def secondFF(f, var_list):
    f_i = []
    n = len(var_list)
    for i in range(0, n):
        f_i.append(partial(f, var_list[i]))

    f_ij = []

    for i in range(0, len(f_i)):
        for j in range(0, n):
            f_ij.append(partial(f_i[i], var_list[j]))

    nu = gauss_map(f, var_list)

    h_elem = []

    for i in range(0, len(f_ij)):
        h_elem.append(dot(nu, f_ij[i]))

    h = Matrix(n, n, h_elem)
    return h

def thirdFF(f, var_list):
    nu = gauss_map(f, var_list)

    nu_i = []
    for i in range(0, len(var_list)):
        nu_i.append(partial(nu, var_list[i]))

    ele = []
    for i in range(0, len(var_list)):
        for j in range(0, len(var_list)):
            ele.append(dot(nu_i[i], nu_i[j]))

    e = Matrix(len(var_list), len(var_list), ele)
    return e

# Calculate Gaussian Curvature
def gaussianCurvature(f, var_list):
    if f.dim != 2:
        print("Surface must be in R^3 to calculate Gaussian Curvature.")
        return -1
    g = firstFF(f, var_list)
    h = secondFF(f, var_list)

    K = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(h.det() / g.det()))))
    return K

# Calculate Christoffel Symbols of First Kind (for Surfaces in R^3) G_ij,k has the array indices [k][i][j]
def Christoffel_1(f, var_list):
    g = firstFF(f, var_list)
    gamma = []
    for i in range(0, len(var_list)):
        g_ele = []
        for j in range(0, len(var_list)):
            subele = []
            for k in range(0, len(var_list)):
                ele = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig((1/2) * (sym.diff(g.matrix[i][j], var_list[k]) + sym.diff(g.matrix[i][k], var_list[j]) - sym.diff(g.matrix[j][k], var_list[i]))))))
                subele.append(ele)
            g_ele.append(subele)
        gamma.append(g_ele)

    return gamma

###
# Calculate Christoffel Symbols of the First Kind Given the Metric Tensor g_ij
###

def Christoffel_1_metric(g, var_list):
    gamma = []
    for i in range(0, len(var_list)):
        g_ele = []
        for j in range(0, len(var_list)):
            subele = []
            for k in range(0, len(var_list)):
                ele = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig((1/2) * (sym.diff(g.matrix[i][j], var_list[k]) + sym.diff(g.matrix[i][k], var_list[j]) - sym.diff(g.matrix[j][k], var_list[i]))))))
                subele.append(ele)
            g_ele.append(subele)
        gamma.append(g_ele)
    
    return gamma

# Calculate Christoffel Symbols of Second Kind (for Surfaces in R^3) G^k_ij has the array indices [k][i][j]
def Christoffel_2(f, var_list):
    Gamma_1 = Christoffel_1(f, var_list)
    g = firstFF(f, var_list)
    g_inv = inverse(g)

    Gamma_2 = []
    for k in range(0, len(var_list)):
        g_ele = []
        for i in range(0, len(var_list)):
            sub_ele = []
            for j in range(0, len(var_list)):
                temp = 0
                for m in range(0, len(var_list)):
                    temp += g_inv.matrix[k][m] * Gamma_1[m][i][j]
                sub_ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(temp)))))
            g_ele.append(sub_ele)
        Gamma_2.append(g_ele)

    return Gamma_2

###
# Calculate Christoffel Symbols of the Second Kind Given the Metric Tensor g_ij
###

def Christoffel_2_metric(g, var_list):
    Gamma_1 = Christoffel_1_metric(g, var_list)
    g_inv = inverse(g)

    Gamma_2 = []
    for k in range(0, len(var_list)):
        g_ele = []
        for i in range(0, len(var_list)):
            sub_ele = []
            for j in range(0, len(var_list)):
                temp = 0
                for m in range(0, len(var_list)):
                    temp += g_inv.matrix[k][m] * Gamma_1[m][i][j]
                sub_ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(temp)))))
            g_ele.append(sub_ele)
        Gamma_2.append(g_ele)

    return Gamma_2


# Print Matrix
def printMatrix(M):
    if not isinstance(M, Matrix):
        print("Input must be a matrix.")

    print("Matrix = ")
    for row in M.matrix:
        print(row)

# Calculate the inverse of the Matrix A
def inverse(A):
    if not A.isSquare:
        print("Matrix must be a square matrix.")
        return Matrix(1, 1, [-1])

    if A.rows == 2:
        D = A.det()
        ele = []
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(A.matrix[1][1] / D)))))
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(-1 * A.matrix[0][1] / D)))))
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(-1 * A.matrix[1][0] / D)))))
        ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(A.matrix[0][0] / D)))))
        A_inv = Matrix(2, 2, ele)
        return A_inv
    else:
        Id = identity(A.rows)
        AM = augmented_matrix(A, Id)
        RREF = RRE_AM(AM)

        ele = []
        for i in range(0, A.rows):
            for j in range(A.cols, AM.cols):
                ele.append(RREF.matrix[i][j])

        A_inv = Matrix(A.rows, A.cols, ele)
        return A_inv

# Multiplies a matrix by a Constant
def MC_mult(M, const):
    if not isinstance(M, Matrix):
        print("Must have a matrix input.")
        return Matrix(1, 1, [-1])

    ele = []
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(M.matrix[i][j] * const)))))
    M_new = Matrix(M.rows, M.cols, ele)
    return M_new

# Matrix times matrix
def MM_mult(A, B):
    # A and B must be of type Matrix
    if not isinstance(A, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if not isinstance(B, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if A.cols != B.rows:
        print("Matrices are not compatible for multiplication.")
        return Matrix(1, 1, [-1])

    n = A.rows
    m = B.cols

    ele = []
    for i in range(0, n):
        for j in range(0, m):
            start = 0
            for k in range(0, A.cols):
                start += A.matrix[i][k] * B.matrix[k][j]
            ele.append(start)

    C = Matrix(n, m, ele)
    return C

# Matrix times vector
def MV_mult(M, v):
    print(isinstance(M, Matrix))
    if not isinstance(M, Matrix):
        print("Must have a matrix input.")
        return Vector([-1])
    if not isinstance(v, Vector):
        print("Must have a vector input.")
        return Vector([-1])
    if M.cols != v.dim:
        print("Matrix and vector are not compatible for multiplication")
        return Vector([-1])
    if v.isrow:
        print("Vector must be a column vector.")
        return Vector([-1])

    ele_list = []
    for i in range(0, M.rows):
        ele = 0
        for j in range(0, M.cols):
            ele += M.matrix[i][j] * v.components[j]
        ele_list.append(ele)

    w = Vector(ele_list)
    return w

# Transposed vector times matrix
def VTM_mult(vT, M):
    if not isinstance(M, Matrix):
        print("Must have a matrix input.")
        return Vector([-1])
    if not isinstance(vT, Vector):
        print("Must have a vector input.")
        return Vector([-1])
    if M.rows != vT.dim:
        print("Matrix and vector are not compatible for multiplication")
        return Vector([-1])
    if not v.isrow:
        print("Vector must be a row vector.")
        return Vector([-1])

    ele_list = []
    for i in range(0, M.cols):
        ele = 0
        for j in range(0, M.rows):
            ele += M.matrix[j][i] * v.components[j]
        ele_list.append(ele)

    w = Vector(ele_list)
    return w

# Add/Subtract Matrices
def M_add(A, B):
    if not isinstance(A, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if not isinstance(B, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if A.cols != B.cols or A.rows != B.rows:
        print("Matrices are not compatible.")
        return Matrix(1, 1, [-1])

    ele = []
    for i in range(0, A.rows):
        for j in range(0, A.cols):
            ele.append(A.matrix[i][j] + B.matrix[i][j])
    C = Matrix(A.rows, A.cols, ele)
    return C

def M_sub(A, B):
    if not isinstance(A, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if not isinstance(B, Matrix):
        print("One or both inputs are not matrices.")
        return Matrix(1, 1, [-1])
    if A.cols != B.cols or A.rows != B.rows:
        print("Matrices are not compatible.")
        return Matrix(1, 1, [-1])

    ele = []
    for i in range(0, A.rows):
        for j in range(0, A.cols):
            ele.append(A.matrix[i][j] - B.matrix[i][j])
    C = Matrix(A.rows, A.cols, ele)
    return C

# Add/Subtract Vectors:
def V_add(vec1, vec2):
    if not isinstance(vec1, Vector):
        print("One or both inputs are not vectors.")
        return Vector([-1])
    if not isinstance(vec2, Vector):
        print("One or both inputs are not vectors.")
        return Vector([-1])
    if not vec1.isrow == vec2.isrow:
        print("Vectors are not compatible.")
        return Vector([-1])

    ele = []
    for i in range(0, vec1.dim):
        ele.append(vec1.components[i] + vec2.components[i])
    v = Vector(ele, vec1.isrow)
    return v

def V_sub(vec1, vec2):
    if not isinstance(vec1, Vector):
        print("One or both inputs are not vectors.")
        return Vector([-1])
    if not isinstance(vec2, Vector):
        print("One or both inputs are not vectors.")
        return Vector([-1])
    if not vec1.isrow == vec2.isrow:
        print("Vectors are not compatible.")
        return Vector([-1])

    ele = []
    for i in range(0, vec1.dim):
        ele.append(vec1.components[i] - vec2.components[i])
    v = Vector(ele, vec1.isrow)
    return v

# Compare Matrices or Vectors
def isEqualV(vec1, vec2):
    if vec1.dim != vec2.dim:
        return False

    if not vec1.isrow == vec2.isrow:
        return False

    for i in range(0, vec1.dim):
        if vec1.components[i] != vec2.components[i]:
            return False
    return True

def isEqualM(A, B):
    if A.rows != B.rows:
        return False
    if A.cols != B.cols:
        return False

    for i in range(0, A.rows):
        for j in range(0, A.cols):
            if A.matrix[i][j] != B.matrix[i][j]:
                return False
    return True

def isZeroM(M):
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            if M.matrix[i][j] != 0:
                return False
    return True

def isZeroV(vec):
    for i in range(0, vec.dim):
        if vec.components[i] != 0:
            return False
    return True

# Differential Map of Parametrized Surface in R^3
def differential(f, var_list):
    f_i = []
    for i in range(0, len(var_list)):
        f_i.append(partial(f, var_list[i]))

    ele = []
    for i in range(0, f.dim):
        for j in range(0, len(var_list)):
            ele.append(f_i[j].components[i])

    Df = Matrix(f.dim, len(var_list), ele)
    return Df

###
# Row Operations
###

# Swap rows (row1 and row2 are the indices of the rows you want to swap)
def row_swap(M, row1, row2):
    ele = []
    for i in range(M.rows):
        for j in range(M.cols):
            if i == row1:
                ele.append(M.matrix[row2][j])
            elif i == row2:
                ele.append(M.matrix[row1][j])
            else:
                ele.append(M.matrix[i][j])

    M_swap = Matrix(M.rows, M.cols, ele)
    return M_swap

# Multiply a row by a constant
def row_const_mult(M, row, c):
    ele = []
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            if i == row:
                ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(c * M.matrix[row][j])))))
            else:
                ele.append(M.matrix[i][j])
    M_mult = Matrix(M.rows, M.cols, ele)
    return M_mult

# Divide a row by a constant
def row_const_div(M, row, c):
    ele = []
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            if i == row:
                ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(M.matrix[row][j] / c)))))
            else:
                ele.append(M.matrix[i][j])
    M_mult = Matrix(M.rows, M.cols, ele)
    return M_mult

# Add one row to another row
def row_add(M, row_orig, row_add):
    ele = []
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            if i == row_orig:
                ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(M.matrix[row_orig][j] + M.matrix[row_add][j])))))
            else:
                ele.append(M.matrix[i][j])

    M_add = Matrix(M.rows, M.cols, ele)
    return M_add

# Subtract one row from another row
def row_sub(M, row_orig, row_sub):
    ele = []
    for i in range(0, M.rows):
        for j in range(0, M.cols):
            if i == row_orig:
                ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(M.matrix[row_orig][j] - M.matrix[row_sub][j])))))
            else:
                ele.append(M.matrix[i][j])

    M_sub = Matrix(M.rows, M.cols, ele)
    return M_sub

# Create Augmented Matrix
def augmented_matrix(A, B):
    if isinstance(B, Vector):
        ele = []
        for i in range(0, A.rows):
            for j in range(0, A.cols + 1):
                if j == A.cols:
                    ele.append(B.components[i])
                else:
                    ele.append(A.matrix[i][j])

        AM = Matrix(A.rows, A.cols + 1, ele)
        return AM
    if isinstance(B, Matrix):
        ele = []
        for i in range(0, A.rows):
            for j in range(0, A.cols + B.cols):
                if j >= A.cols:
                    ele.append(B.matrix[i][j - A.cols])
                else:
                    ele.append(A.matrix[i][j])

        AM = Matrix(A.rows, A.cols + B.cols, ele)
        return AM


# Calculate Reduce Row Echelon Form of Augmented Matrix (i.e. solving Ax = b)
def RRE_AM(M):
    if isZeroM(M):
        return M
    lead = 0
    rowCount = M.rows
    columnCount = M.cols
    for r in range(0, rowCount):
        if lead >= columnCount:
            return M
        i = r
        while M.matrix[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return M
        M = row_swap(M, i, r)
        lv = M.matrix[r][lead]
        M = row_const_div(M, r, lv)
        for i in range(0, rowCount):
            if i != r:
                lv = M.matrix[i][lead]
                if lv != 0:
                    M = row_const_mult(M, r, lv)
                    M = row_sub(M, i, r)
                    M = row_const_div(M, r, lv)
        lead += 1

    return M

# Checks whether a vector field lies in the tangent plane
def isInTP(f, var_list, X):
    Df = differential(f, var_list)
    AM = augmented_matrix(Df, X)
    RREF = RRE_AM(AM)
    return RREF.matrix[-1][-1] == 0

# Change basis of X from standard basis to basis of the Tangent Plane
def TP_Basis(f, var_list, X):
    Df = differential(f, var_list)
    AM = augmented_matrix(Df, X)
    RREF = RRE_AM(AM)
    if not isInTP(f, var_list, X):
        print("Vector field does not lie in the tangent plane.")
        return Vector([-1])
    coeff = []
    for i in range(0, RREF.rows - 1):
        coeff.append(RREF.matrix[i][-1])
    X_tp = Vector(coeff)
    return X_tp

# Calculate first fundamental form acting on two vectors in the tangent plane
def I(f, var_list, X, Y):
    g = firstFF(f, var_list)
    v = MV_mult(g, Y)
    val = dot(X, v)
    return val

# Calculate second fundamental form acting on two vectors in tangent plane
def II(f, var_list, X, Y):
    h = secondFF(f, var_list)
    v = MV_mult(h, Y)
    val = dot(X, v)
    return val

# Calculate Covariant Derivative of vector field Y in direction of vector field X (both in basis of tangent space)
def Covariant_Deriv(f, var_list, X, Y):
    Gamma = Christoffel_2(f, var_list)
    ele = []
    for k in range(0, len(var_list)):
        component = 0
        for i in range(0, len(var_list)):
            component += sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(X.components[i] * sym.diff(Y.components[k], var_list[i])))))
            for j in range(0, len(var_list)):
                component += sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(X.components[i] * Y.components[j] * Gamma[k][i][j]))))
        ele.append(sym.simplify(symp.powsimp(sym.trigsimp(sym.expand_trig(component)))))

    DelXY = Vector(ele)
    return DelXY

# Calculates the Riemann Curvature Tensor R^rho_sigma,mu,nu
def Riemann_Curvature(f, var_list):
    Gamma = Christoffel_2(f, var_list)
    R = []
    for rho in range(0, len(var_list)):
        ele1 = []
        for sigma in range(0, len(var_list)):
            ele2 = []
            for mu in range(0, len(var_list)):
                ele3 = []
                for nu in range(0, len(var_list)):
                    val = 0
                    val += sym.diff(Gamma[rho][nu][sigma], var_list[mu])
                    val -= sym.diff(Gamma[rho][mu][sigma], var_list[nu])
                    for lamb in range(0, len(var_list)):
                        val += Gamma[rho][mu][lamb] * Gamma[lamb][nu][sigma]
                        val -= Gamma[rho][nu][lamb] * Gamma[lamb][mu][sigma]
                    ele3.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
                ele2.append(ele3)
            ele1.append(ele2)
        R.append(ele1)

    return R

###
# Calculates the Riemann Curvature Tensor given the Metric Tensor g_ij
###
def Riemann_Curvature_metric(g, var_list):
    Gamma = Christoffel_2_metric(g, var_list)
    R = []
    for rho in range(0, len(var_list)):
        ele1 = []
        for sigma in range(0, len(var_list)):
            ele2 = []
            for mu in range(0, len(var_list)):
                ele3 = []
                for nu in range(0, len(var_list)):
                    val = 0
                    val += sym.diff(Gamma[rho][nu][sigma], var_list[mu])
                    val -= sym.diff(Gamma[rho][mu][sigma], var_list[nu])
                    for lamb in range(0, len(var_list)):
                        val += Gamma[rho][mu][lamb] * Gamma[lamb][nu][sigma]
                        val -= Gamma[rho][nu][lamb] * Gamma[lamb][mu][sigma]
                    ele3.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
                ele2.append(ele3)
            ele1.append(ele2)
        R.append(ele1)

    return R

# Calculates the Ricci Curvature Tensor for a Parametrized Surface
def Ricci_Curvature(f, var_list):
    n = len(var_list)
    Riemann = Riemann_Curvature(f, var_list)

    Ric_ele = []
    for i in range(0, n):
        for j in range(0, n):
            val = 0
            for k in range(0, n):
                val += Riemann[k][i][k][j]
            Ric_ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
    Ric = Matrix(n, n, Ric_ele)
    return Ric

###
# Calculates the Ricci Curvature Tensor given the Metric Tensor g_ij
###
def Ricci_Curvature_metric(g, var_list):
    n = len(var_list)
    Riemann = Riemann_Curvature_metric(g, var_list)

    Ric_ele = []
    for i in range(0, n):
        for j in range(0, n):
            val = 0
            for k in range(0, n):
                val += Riemann[k][i][k][j]
            Ric_ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
    Ric = Matrix(n, n, Ric_ele)
    return Ric

# Calculate the Ricci Scalar for a Parametrized Surface
def Ricci_Scalar(f, var_list):
    Ric = Ricci_Curvature(f,var_list)
    g = firstFF(f, var_list)
    g_inv = inverse(g)
    R = 0
    for i in range(0, len(var_list)):
        for j in range(0, len(var_list)):
            R += g_inv.matrix[i][j] * Ric.matrix[i][j]
    R = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(R))))
    return R

###
# Calculates the Ricci Scalar given the Metric Tensor g_ij
###
def Ricci_Scalar_metric(g, var_list):
    Ric = Ricci_Curvature_metric(g,var_list)
    g_inv = inverse(g)
    R = 0
    for i in range(0, len(var_list)):
        for j in range(0, len(var_list)):
            R += g_inv.matrix[i][j] * Ric.matrix[i][j]
    R = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(R))))
    return R

# Returns the Minkowski Metric for n dimensional spacial components (-, +, +, +, +, ...)
def Minkowski(n):
    ele = []
    for i in range(0, n+1):
        for j in range(0, n+1):
            if i == j and i == 0:
                ele.append(-1)
            elif i == j and i != 0:
                ele.append(1)
            else:
                ele.append(0)
    eta = Matrix(n+1, n+1, ele)
    return eta

# Calculates the Einstein Tensor for a parametrized surface in R^n_1 with n-1 spacial parameters and 1 time parameter (assumes Minkowski metric is not added yet)
def Einstein_Tensor(f, var_list, Cos_const=True):
    if len(var_list) != f.dim:
        print("You must have the same number of parameters as the vector's dimension.")
        return Matrix(1, 1, [-1])
    if not Cos_const:
        g = M_add(Minkowski(f.dim - 1), firstFF(f, var_list))
        Ein = Einstein_Tensor_metric(g, var_list, Cos_const)
        return Ein
    else:
        g = M_add(Minkowski(f.dim - 1), firstFF(f, var_list))
        Ein = Einstein_Tensor_metric(g, var_list, Cos_const)
        return Ein

###
# Calculates the Einstein Tensor given the Metric Tensor g_ij (assumes already has Minkowski metric added to it)
###

def Einstein_Tensor_metric(g, var_list, Cos_const=True):
    if not Cos_const:
        Ric = Ricci_Curvature_metric(g, var_list)
        R = Ricci_Scalar_metric(g, var_list)
        Ein = M_sub(Ric, MC_mult(g, (R/2)))
        return Ein
    else:
        Ric = Ricci_Curvature_metric(g, var_list)
        R = Ricci_Scalar_metric(g, var_list)
        Ein = M_add(M_sub(Ric, MC_mult(g, (R/2))), MC_mult(g, LAMBDA))
        return Ein

# Calculates the Energy Momentum Tensor given a parametrized surface in R^n_1 with n-1 spacial parameters and 1 time parameter
def EMT(f, var_list, Cos_const=True):
    Ein = Einstein_Tensor(f, var_list, Cos_const)
    T = MC_mult(Ein, 1/KAPPA)
    return T

###
# Calculates the Energy Momentum Tensor given the Metric Tensor g_ij
###
def EMT_metric(g, var_list, Cos_const=True):
    Ein = Einstein_Tensor_metric(g, var_list, Cos_const)
    T = MC_mult(Ein, 1/KAPPA)
    return T

###
# Calculate the Kulkarni-Nomizu Product of two (0, 2) tensors A and B
###
def Kulkarni_Nomizu(A, B):
    n = A.rows
    KN = []
    for rho in range(0, n):
        ele1 = []
        for sigma in range(0, n):
            ele2 = []
            for mu in range(0, n):
                ele3 = []
                for nu in range(0, n):
                    val = A.matrix[rho][mu] * B.matrix[sigma][nu] + A.matrix[sigma][nu] * B.matrix[rho][mu] - A.matrix[rho][nu] * B.matrix[sigma][mu] - A.matrix[sigma][mu] * B.matrix[rho][nu]
                    ele3.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
                ele2.append(ele3)
            ele1.append(ele2)
        KN.append(ele1)

    return KN

###
# Calculate the Weyl Tensor for a Parametrized Surface
###
def Weyl_Tensor(f, var_list):
    n = len(var_list)
    Riemann = Riemann_Curvature(f, var_list)
    Ric = Ricci_Curvature(f, var_list)
    R = Ricci_Scalar(f, var_list)
    g = firstFF(f, var_list)

    term1 = Kulkarni_Nomizu(g, Ric)
    for  i in range(0, n):
        for j in range(0, n):
            for k in range(0, n):
                for l in range(0, n):
                    term1[i][j][k][l] = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(term1[i][j][k][l] * (1/(n - 2))))))
    term2 = Kulkarni_Nomizu(g, g)
    for  i in range(0, n):
        for j in range(0, n):
            for k in range(0, n):
                for l in range(0, n):
                    term2[i][j][k][l] = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(term2[i][j][k][l] * (R / (2 * (n - 1) * (n - 2)))))))

    Weyl = []
    for rho in range(0, n):
        ele1 = []
        for sigma in range(0, n):
            ele2 = []
            for mu in range(0, n):
                ele3 = []
                for nu in range(0, n):
                    val = Riemann[rho][sigma][mu][nu] - term1[rho][sigma][mu][nu] + term2[rho][sigma][mu][nu]
                    ele3.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
                ele2.append(ele3)
            ele1.append(ele2)
        Weyl.append(ele1)

    return Weyl

###
# Calculate the Weyl Tensor given the Metric Tensor g_ij
###
def Weyl_Tensor_metric(g, var_list):
    n = len(var_list)
    Riemann = Riemann_Curvature_metric(g, var_list)
    Ric = Ricci_Curvature_metric(g, var_list)
    R = Ricci_Scalar_metric(g, var_list)

    term1 = Kulkarni_Nomizu(g, Ric)
    for  i in range(0, n):
        for j in range(0, n):
            for k in range(0, n):
                for l in range(0, n):
                    term1[i][j][k][l] = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(term1[i][j][k][l] * (1/(n - 2))))))
    term2 = Kulkarni_Nomizu(g, g)
    for  i in range(0, n):
        for j in range(0, n):
            for k in range(0, n):
                for l in range(0, n):
                    term2[i][j][k][l] = sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(term2[i][j][k][l] * (R / (2 * (n - 1) * (n - 2)))))))

    Weyl = []
    for rho in range(0, n):
        ele1 = []
        for sigma in range(0, n):
            ele2 = []
            for mu in range(0, n):
                ele3 = []
                for nu in range(0, n):
                    val = Riemann[rho][sigma][mu][nu] - term1[rho][sigma][mu][nu] + term2[rho][sigma][mu][nu]
                    ele3.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
                ele2.append(ele3)
            ele1.append(ele2)
        Weyl.append(ele1)

    return Weyl

###
# Calculates the Schouten Tensor for a Parametrized Surface
###
def Schouten_Tensor(f, var_list):
    n = len(var_list)
    Ric = Ricci_Curvature(f, var_list)
    R = Ricci_Scalar(f, var_list)
    g = firstFF(f, var_list)

    ele = []
    for i in range(0, n):
        for j in range(0, n):
            val = (1 / (n - 2)) * (Ric.matrix[i][j] - (R / (2 * (n - 1))) * g.matrix[i][j])
            ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
    S = Matrix(n, n, ele)
    return S

###
# Calculates the Schouten Tensor given the Metric Tensor g_ij
###
def Schouten_Tensor_metric(g, var_list):
    n = len(var_list)
    Ric = Ricci_Curvature_metric(g, var_list)
    R = Ricci_Scalar_metric(g, var_list)

    ele = []
    for i in range(0, n):
        for j in range(0, n):
            val = (1 / (n - 2)) * (Ric.matrix[i][j] - (R / (2 * (n - 1))) * g.matrix[i][j])
            ele.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
    S = Matrix(n, n, ele)
    return S

###
# Calculates the Cotton Tensor for a Parametrized Surface
###
def Cotton_Tensor(f, var_list):
    n = len(var_list)
    S = Schouten_Tensor(f, var_list)

    C = []
    for rho in range(0, n):
        ele1 = []
        for sigma in range(0, n):
            ele2 = []
            for mu in range(0, n):
                val = sym.diff(S.matrix[sigma][rho], var_list[mu]) - sym.diff(S.matrix[mu][rho], var_list[sigma])
                ele2.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
            ele1.append(ele2)
        C.append(ele1)

    return C

###
# Calculates the Cotton Tensor given the Metric Tensor g_ij
###
def Cotton_Tensor_metric(g, var_list):
    n = len(var_list)
    S = Schouten_Tensor_metric(g, var_list)

    C = []
    for rho in range(0, n):
        ele1 = []
        for sigma in range(0, n):
            ele2 = []
            for mu in range(0, n):
                val = sym.diff(S.matrix[sigma][rho], var_list[mu]) - sym.diff(S.matrix[mu][rho], var_list[sigma])
                ele2.append(sym.simplify(sym.powsimp(sym.trigsimp(sym.expand_trig(val)))))
            ele1.append(ele2)
        C.append(ele1)

    return C
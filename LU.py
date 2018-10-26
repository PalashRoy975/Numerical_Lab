import numpy as np
from numpy.linalg import inv
from sympy import *

f = open("output_LU.txt", "w+")


def equation_to_matrix(filename):
    unique_symbols = set([])

    file = open(filename, "r")
    data = file.read()

    data = data.replace(" ", "")
    data = data.replace("=", "-")  # replace = with minus
   
    for i in range(len(data)):
        if (data[i] >= 'a' and data[i] <= 'z'):
            if ((i + 1) < len(data) and (data[i + 1]>='0' and data[i + 1] <= '9')):
                unique_symbols.add(str(data[i] + data[i + 1]))
            else:
                unique_symbols.add(data[i])

        elif data[i] >= 'A' and data[i] <= 'Z':
            unique_symbols.add(data[i]);

    temp = sorted(unique_symbols)  # sort symbols
    ret3 = temp

    notation = ""

    for i in range(len(ret3)):

        if i != 0:
            notation = notation + " "

        data = data.replace(ret3[i], "x" + str(i + 1))
        notation = notation + "x" + str(i + 1)

    i = 0
    while True:
        if (i + 1) < len(data) and (data[i] >= '0' and data[i] <= '9') and data[i + 1] == 'x':
            data = data[:i + 1] + "*" + data[i + 1:]
        if (i > len(data)):
            break

        i = i + 1

    variables = symbols(notation)

    data = data.splitlines()

    equations = data

    A, B = linear_eq_to_matrix(equations, variables)

    A = A.tolist()
    b = []

    for temp in B:
        b.append(temp)
    n = len(b)

    return A, b, ret3, n


def swap(a, b):
    temp = a
    a = b
    b = temp


def forward(L, U, pos, n, B):
    # partial pivoting not used
    process(U, pos, n, B, L)
    for i in range(pos + 1, n):
        dd = U[i][pos] / U[pos][pos]
        for j in range(0, n):
            U[i][j] = U[i][j] - dd * U[pos][j]
        U[i][pos] = 0
        L[i][pos] = dd
    print('Upper Triangular matrix creation step', file=f)
    print(U,file=f)
    print('Lower Triangular matrix creation step', file=f)
    print(L,file=f)

def process(U, pos, n, B, L):
    mx = pos
    for i in range(pos, n):
        if abs(U[i][pos]) > abs(U[mx][pos]):
            mx = i
    #swap(B[pos][0],B[mx][0])
    d = B[pos][0]
    B[pos][0] = B[mx][0]
    B[mx][0] = d
    for i in range(0, n):
        d = U[pos][i]
        U[pos][i] = U[mx][i]
        U[mx][i] = d
        d = L[pos][i]
        L[pos][i] = L[mx][i]
        L[mx][i] = d
    #print('pivoting')
    #print(U)
    #print(pos)
    #print(mx)

def main():
    matrix, val, variables, sz = equation_to_matrix("input_Gauss_Seidal.txt")

    print(matrix)
    print(val)
    print(variables)
    print(sz)
    for i in range(0,sz):
        for j in range(0,sz):
            print(matrix[i][j])

    mm = np.zeros(shape=(sz, sz))
    for i in range(0, sz):
        for j in range(0, sz):
            mm[i][j] = matrix[i][j]
    X = mm
    mm = np.zeros(shape=(sz, 1))
    for i in range(0, sz):
        mm[i][0] = val[i]
    n = sz
    B = mm
    U = X
    L = np.zeros(shape=(n, n))
    # print('Upper Triangular matrix creation steps', file=f)
    # print(U, file=f)
    # print('Lower Triangular matrix creation step', file=f)
    # print(L, file=f)
    for i in range(0, n):
        forward(L, U, i, n, B)
    for i in range(0, n):
        for j in range(0, n):
            L[i][i] = 1
    L_inv = inv(L)
    U_inv = inv(U)

    print('Final Upper Triangular matrix', file=f)
    print(U, file=f)
    print('Final Lower Triangular matrix', file=f)
    print(L, file=f)
    print('Final Upper Triangular inverse matrix', file=f)
    print(U_inv, file=f)
    print('Final Lower Triangular inverse matrix', file=f)
    print(L_inv, file=f)
    Z = np.matmul(L_inv, B)
    print('Final Matrix Z', file=f)
    print(Z, file=f)
    X = np.matmul(U_inv, Z)
    print('Final Matrix X', file=f)
    print(X, file=f)
    for i in range(0, n):
        print(variables[i] + ' = ', end='', file=f)
        print(round(X[i][0], 5), file=f)
    print(X)
    f.close()


main()

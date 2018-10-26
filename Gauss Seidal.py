import numpy as np
from numpy.linalg import inv
from sympy import *
import matplotlib.pyplot as plt

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

def Gauss_Seidal(mat,val):
    count =0
    X_Input,temp,temp2 = [],[],[]
    print("Enter Initial Guess:")
    for i in range(len(val)):
        X_Input.append(int(input()))
        temp.append(X_Input[i])
    X1_arr, Xerror_arr, X3_arr = [], [], []
    while True:
        for i in range(len(val)):
            a = val [i]
            for j in range(len(val)):
                if(i!=j):
                    a -=mat[i][j]*X_Input[j]
            a/=mat[i][i]
            X_Input[i]=a
            if(count!=0):
                err = abs(temp[i]-X_Input[i])/abs(X_Input[i])*100
                Xerror_arr.append(err)
            X1_arr.append(a)


        print("Iteration: " + str(count + 1))
        print(X1_arr)
        print("Errors In Iteration: " + str(count + 1)+" In Percentage")
        print(Xerror_arr)
        Xerror_arr.clear()
        temp.clear()
        for i in range(len(val)):
            temp.append(X_Input[i])

        X1_arr.clear()
        count=count+1
        if(count==10): break

def main():
    matrix, val, variables, sz = equation_to_matrix("input_Gauss_Seidal.txt")

    print(matrix)
    print(val)
    Gauss_Seidal(matrix,val)

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

    f.close()

main()

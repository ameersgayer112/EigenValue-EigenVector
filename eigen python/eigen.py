from numpy import array, identity, diagonal
from math import sqrt
import numpy as np
from timeit import timeit
import sys


# Find largest off-diagonal element a[k,l]
def maxElem(A):
    n = len(A)
    Amax = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if abs(A[i, j]) >= Amax:
                Amax = abs(A[i, j])
                k = i;
                l = j
    return Amax, k, l


def rotate(A, p, k, l):
    n = len(A)
    # print(n)
    Adiff = A[l, l] - A[k, k]
    if abs(A[k, l]) < abs(Adiff) * 1.0e-36:
        t = A[k, l] / Adiff
    else:
        phi = Adiff / (2.0 * A[k, l])
        t = 1.0 / (abs(phi) + sqrt(phi ** 2 + 1.0))
        if phi < 0.0:
            t = -t
    c = 1.0 / sqrt(t ** 2 + 1.0);
    s = t * c
    tau = s / (1.0 + c)
    temp = A[k, l]
    A[k, l] = 0.0
    A[k, k] = A[k, k] - t * temp
    A[l, l] = A[l, l] + t * temp
    for i in range(k):  # Case of i < k
        temp = A[i, k]
        A[i, k] = temp - s * (A[i, l] + tau * temp)
        A[i, l] = A[i, l] + s * (temp - tau * A[i, l])
    for i in range(k + 1, l):  # Case of k < i < l
        temp = A[k, i]
        A[k, i] = temp - s * (A[i, l] + tau * A[k, i])
        A[i, l] = A[i, l] + s * (temp - tau * A[i, l])
    for i in range(l + 1, n):  # Case of i > l
        temp = A[k, i]
        A[k, i] = temp - s * (A[l, i] + tau * temp)
        A[l, i] = A[l, i] + s * (temp - tau * A[l, i])

    for i in range(n):  # Update transformation matrix
        temp = p[i, k]
        p[i, k] = temp - s * (p[i, l] + tau * p[i, k])
        p[i, l] = p[i, l] + s * (temp - tau * p[i, l])
    p = p.T


def Jacobi(A, tol=0):
    n = len(A)
    maxRot = 5 * (n ** 2)  # Set limit on number of rotations
    p = identity(n) * 1.0  # Initialize transformation matrix
    for i in range(2000):  # Jacobi rotation loop
        print("itiration number :", i)
        Amax, k, l = maxElem(A)
        if Amax < tol: return diagonal(A), p
        rotate(A, p, k, l)
    return diagonal(A), p

def check_symmetric(mat):
    r,c= mat.shape
    print (r,c)
    if(r!=c): return False

    for i in range(r):
        for j in range(r):
            if (mat[i][j] != mat[j][i]):
                return False
    return True


def main():
    A = np.loadtxt("C:\data_8.txt", dtype=np.float, delimiter=',')
    np.set_printoptions(threshold=sys.maxsize)  # to print all of the matrix
    #A = np.array([[4, -14, -12], [-14, 10, 13], [-12, 13, 1]])
    if check_symmetric(A) is True:
        eigenvalues, eigenvectors = Jacobi(A, tol=1.0e-9)

        f = open("results.txt", 'w')

        f.write("------------------------------------------- The Eigenvectors  matrix ---------------------------\n")
        print("------------------------------------------- The Eigenvectors matrix  ---------------------------\n")
        # print(eigenvectors)
        f.write(str(eigenvectors) + "\n")
        print("-------------------------------------------The Eigenvalues   ---------------------------\n")
        f.write("-------------------------------------------The Eigenvalues  ---------------------------\n")
        # print(eigenvalues)
        f.write(str(eigenvalues))
        f.close()
    else:
        print("the matrix is not symmtrical")



print("the total time is : ", timeit(main, number=1))


def m2():
    A = np.loadtxt("C:\data_8.txt", dtype=np.float, delimiter=',')
    print (A.shape)


    print(np.linalg.eig(A))
#print("the total time is : ", timeit(m2, number=1))

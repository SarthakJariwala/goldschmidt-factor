import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import ABX_functions

# https://onlinelibrary-wiley-com.offcampus.lib.washington.edu/doi/epdf/10.1002/aenm.201902467
def goldschmidt(A, B, X, A_ratio=None, B_ratio=None, X_ratio=None):
    '''

    '''
    A_dictionary = {"MA" : 2.16, "FA" : 2.53, "EA" : 2.74, "Cs" : 1.67}
    B_dictionary = {"Pb" : 1.19, "Sn" : 1.15}
    X_dictionary = {"I" : 2.20, "Br" : 1.96, "Cl" : 1.84}

    if len(A) == len(B) == len(X) == 1:
        return ABX_functions.ABX(A, B, X, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 2 and len(B) == 1 and len(X) == 1:
        return ABX_functions.A2BX(A, B, X, A_ratio=A_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 1 and len(B) == 2 and len(X) == 1:
        return ABX_functions.AB2X(A, B, X, B_ratio=B_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 1 and len(B) == 1 and len(X) == 2:
        return ABX_functions.ABX2(A, B, X, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 1 and len(B) == 1 and len(X) == 3:
        return ABX_functions.ABX3(A, B, X, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 3 and len(B) == 1 and len(X) == 1:
        return ABX_functions.A3BX(A, B, X, A_ratio=A_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 2 and len(B) == 1 and len(X) == 2:
        return ABX_functions.A2BX2(A, B, X, A_ratio=A_ratio, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 3 and len(B) == 1 and len(X) == 2:
        return ABX_functions.A3BX2(A, B, X, A_ratio=A_ratio, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 2 and len(B) == 1 and len(X) == 3:
        return ABX_functions.A2BX3(A, B, X, A_ratio=A_ratio, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

    elif len(A) == 3 and len(B) == 1 and len(X) == 3:
        return ABX_functions.A3BX3(A, B, X, A_ratio=A_ratio, X_ratio=X_ratio, A_dictionary=A_dictionary, B_dictionary=B_dictionary, X_dictionary=X_dictionary)

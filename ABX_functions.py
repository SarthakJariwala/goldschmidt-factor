import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def ABX(A, B, X, A_dictionary, B_dictionary, X_dictionary):
    A_radii = A_dictionary.get(A[0])
    B_radii = B_dictionary.get(B[0])
    X_radii = X_dictionary.get(X[0])

    t_effective = (A_radii + X_radii) / (math.sqrt(2) * (B_radii + X_radii))
    return t_effective

def ABX2(A, B, X, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if X_ratio == None:
        X_radii = []
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        A_radii = A_dictionary.get(A[0])
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)
        r_effective = ratio * X_radii[0] + (1 - ratio) * X_radii[1]
        t_effective = (A_radii + r_effective) / (math.sqrt(2) * (B_radii + r_effective))

        ones = np.ones(ratio.shape)
        eights = ones*0.8
        nines = ones*0.9
        plt.plot(ratio, t_effective, color='red', lw=2)
        plt.plot(ratio, ones, c='black')
        plt.plot(ratio, eights, c='black')
        plt.ylim(0.6, 1.2)
        plt.xlim(0,1)
        title = plt.title("$%s%s_x%s_{3x}%s_{3(1-x)}$ tolerance factor as a function of %s molar fraction" % (A[0], B[0], X[0], X[1], X[0]))
        title.set_position([0.5, 1.05])
        plt.ylabel("Tolerance Factor t", fontsize=14)
        plt.xlabel("%s Molar Ratio" % X[0], fontsize=14)
        plt.fill_between(ratio, nines, eights, color='yellow', alpha=0.5)
        plt.fill_between(ratio, nines, ones, color='green', alpha=0.5)
        plt.show()

        df = pd.DataFrame(np.round(ratio, 2), columns=['%s Ratio' % X[0]])
        df['Tolerance Factor'] = t_effective
        return df

    else:
        if sum(X_ratio) == 1:
            X_radii = []
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            A_radii = A_dictionary.get(A[0])
            B_radii = B_dictionary.get(B[0])

            r_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1]
            t_effective = (A_radii + r_effective) / (math.sqrt(2) * (B_radii + r_effective))
            return t_effective
        else:
            print('Error: The sum of X_ratio is not equal to 1.')

def ABX3(A, B, X, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if X_ratio == None:
        X_radii = []
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        A_radii = A_dictionary.get(A[0])
        B_radii = B_dictionary.get(B[0])

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue
        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)
        r_effective = x_ratio * X_radii[0] + y_ratio * X_radii[1] + z_ratio * X_radii[2]
        t_effective = (A_radii + r_effective) / (math.sqrt(2) * (B_radii + r_effective))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        img = ax.scatter(x_ratio, y_ratio, z_ratio, c=t_effective, cmap=plt.hot())
        fig.colorbar(img)
        ax.set_xlabel("%s Molar Ratio" % X[0])
        ax.set_ylabel("%s Molar Ratio" % X[1])
        ax.set_zlabel("%s Molar Ratio" % X[2])
        title = plt.title("$%s%s%s_x%s_y%s_z$ tolerance factor as a function of halide composition" % (A[0], B[0], X[0], X[1], X[2]))
        title.set_position([0.5,1.05])
        plt.show()

        df = pd.DataFrame(x_ratio, columns =['%s Ratio' % X[0]])
        df['%s Ratio' % X[1]] = y_ratio
        df['%s Ratio' % X[2]] = z_ratio
        df['Tolerance Factor'] = t_effective
        return df

    else:
        if sum(X_ratio) == 1:
            X_radii = []
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            A_radii = A_dictionary.get(A[0])
            B_radii = B_dictionary.get(B[0])

            r_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1] + X_ratio[2] * X_radii[2]
            t_effective = (A_radii + r_effective) / (math.sqrt(2) * (B_radii + r_effective))
            return t_effective
        else:
            print('Error: The sum of X_ratio is not equal to 1.')

def AB2X(A, B, X, B_ratio, A_dictionary, B_dictionary, X_dictionary):
    if B_ratio == None:
        B_radii = []
        for ions in B:
            B_radii.append(B_dictionary.get(ions))
        A_radii = A_dictionary.get(A[0])
        X_radii = X_dictionary.get(X[0])

        ratio = np.linspace(0, 1, num=11)
        r_effective = ratio * B_radii[0] + (1 - ratio) * B_radii[1]
        t_effective = (A_radii + X_radii) / (math.sqrt(2) * (r_effective + X_radii))

        ones = np.ones(ratio.shape)
        eights = ones*0.8
        nines = ones*0.9
        plt.plot(ratio, t_effective, color='red', lw=2)
        plt.plot(ratio, ones, c='black')
        plt.plot(ratio, eights, c='black')
        plt.ylim(0.6, 1.2)
        plt.xlim(0,1)
        title = plt.title("$%s%s_x%s_{1-x}%s_3$ tolerance factor as a function of %s molar fraction" % (A[0], B[0], B[1], X[0], B[0]))
        title.set_position([0.5, 1.05])
        plt.ylabel("Tolerance Factor t", fontsize=14)
        plt.xlabel("%s Molar Ratio" % B[0], fontsize=14)
        plt.fill_between(ratio, nines, eights, color='yellow', alpha=0.5)
        plt.fill_between(ratio, nines, ones, color='green', alpha=0.5)
        plt.show()

        df = pd.DataFrame(np.round(ratio, 2), columns=['%s Ratio' % B[0]])
        df['Tolerance Factor'] = t_effective
        return df
    else:
        if sum(B_ratio) == 1:
            B_radii = []
            for ions in B:
                B_radii.append(B_dictionary.get(ions))
            A_radii = A_dictionary.get(A[0])
            X_radii = X_dictionary.get(X[0])

            r_effective = B_ratio[0] * B_radii[0] + B_ratio[1] * B_radii[1]
            t_effective = (A_radii + X_radii) / (math.sqrt(2) * (r_effective + X_radii))
            return t_effective
        else:
            print('Error: The sum of B_ratio is not equal to 1.')

def A2BX(A, B, X, A_ratio, A_dictionary, B_dictionary, X_dictionary):
    if A_ratio == None:
        A_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])
        X_radii = X_dictionary.get(X[0])

        ratio = np.linspace(0, 1, num=11)
        r_effective = ratio * A_radii[0] + (1 - ratio) * A_radii[1]
        t_effective = (r_effective + X_radii) / (math.sqrt(2) * (B_radii + X_radii))

        ones = np.ones(ratio.shape)
        eights = ones*0.8
        nines = ones*0.9
        plt.plot(ratio, t_effective, color='red', lw=2)
        plt.plot(ratio, ones, c='black')
        plt.plot(ratio, eights, c='black')
        plt.ylim(0.6, 1.2)
        plt.xlim(0,1)
        title = plt.title("$%s_x%s_{1-x}%s%s_3$ tolerance factor as a function of %s molar fraction" % (A[0], A[1], B[0], X[0], A[0]))
        title.set_position([0.5, 1.05])
        plt.ylabel("Tolerance Factor t", fontsize=14)
        plt.xlabel("%s Molar Ratio" % A[0], fontsize=14)
        plt.fill_between(ratio, nines, eights, color='yellow', alpha=0.5)
        plt.fill_between(ratio, nines, ones, color='green', alpha=0.5)
        plt.show()

        df = pd.DataFrame(np.round(ratio, 2), columns=['%s Ratio' % A[0]])
        df['Tolerance Factor'] = t_effective
        return df

    else:
        if sum(A_ratio == 1):
            A_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            B_radii = B_dictionary.get(B[0])
            X_radii = X_dictionary.get(X[0])

            r_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1]
            t_effective = (r_effective + X_radii) / (math.sqrt(2) * (B_radii + X_radii))
            return t_effective
        else:
            print('Error: The sum of A_ratio is not equal to 1.')

def A3BX(A, B, X, A_ratio, A_dictionary, B_dictionary, X_dictionary):
    if A_ratio == None:
        A_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        X_radii = X_dictionary.get(X[0])
        B_radii = B_dictionary.get(B[0])

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue
        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)
        r_effective = x_ratio * A_radii[0] + y_ratio * A_radii[1] + z_ratio * A_radii[2]
        t_effective = (r_effective + X_radii) / (math.sqrt(2) * (B_radii + X_radii))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        img = ax.scatter(x_ratio, y_ratio, z_ratio, c=t_effective, cmap=plt.hot())
        fig.colorbar(img)
        ax.set_xlabel("%s Molar Ratio" % A[0])
        ax.set_ylabel("%s Molar Ratio" % A[1])
        ax.set_zlabel("%s Molar Ratio" % A[2])
        title = plt.title("$%s%s%s_x%s_y%s_z$ tolerance factor as a function of A-site cation composition" % (A[0], A[1], A[2], B[0], X[0]))
        title.set_position([0.5,1.05])
        plt.show()

        df = pd.DataFrame(x_ratio, columns =['%s Ratio' % A[0]])
        df['%s Ratio' % A[1]] = y_ratio
        df['%s Ratio' % A[2]] = z_ratio
        df['Tolerance Factor'] = t_effective
        return df

    else:
        if (sum(A_ratio) == 1):
            A_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            X_radii = X_dictionary.get(X[0])
            B_radii = B_dictionary.get(B[0])

            r_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1] + A_ratio[2] * A_radii[2]
            t_effective = (r_effective + X_radii) / (math.sqrt(2) * (B_radii + X_radii))
            return t_effective
        else:
            print('Error: The sum of A_ratio is not equal to 1.')

def A2BX2(A, B, X, A_ratio, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if (A_ratio == None and X_ratio == None):
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        A_effective = ratio * A_radii[0] + (1-ratio) * A_radii[1]
        X_effective = ratio * X_radii[0] + (1-ratio) * X_radii[1]

        t_effective = []
        for i in A_effective:
                t_effective.append((i + X_effective) / (math.sqrt(2) * (B_radii + X_effective)))

        df = pd.DataFrame(ratio, columns =['%s Ratio' % A[0]])
        df['%s Ratio' % A[1]] = 1-ratio
        #df['Tolerance Factor'] = t_effective
        i_count = 0
        ratio = np.round(ratio, decimals=2)
        for i in ratio:
            df['%s' % i] = t_effective[i_count]
            i_count += 1
        df = df.rename(columns = {'0.0' : '%s Ratio : 0.0' % X[0]})
        return df

    elif ((A_ratio == None and X_ratio != None) or (A_ratio != None and X_ratio == None)):
        print('Warning: Insert a list of ratios for both A_ratio and X_ratio to calculate a specifice tolerance factor')
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        A_effective = ratio * A_radii[0] + (1-ratio) * A_radii[1]
        X_effective = ratio * X_radii[0] + (1-ratio) * X_radii[1]

        t_effective = []
        for i in A_effective:
                t_effective.append((i + X_effective) / (math.sqrt(2) * (B_radii + X_effective)))

        df = pd.DataFrame(ratio, columns =['%s Ratio' % A[0]])
        df['%s Ratio' % A[1]] = 1-ratio
        #df['Tolerance Factor'] = t_effective
        i_count = 0
        ratio = np.round(ratio, decimals=2)
        for i in ratio:
            df['%s' % i] = t_effective[i_count]
            i_count += 1
        df = df.rename(columns = {'0.0' : '%s Ratio : 0.0' % X[0]})
        return df

    elif (A_ratio != None and X_ratio != None):
        if (sum(A_ratio) == 1 and sum(X_ratio_ ==1)):
            A_radii = []
            X_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            B_radii = B_dictionary.get(B[0])

            A_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1]
            X_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1]
            t_effective = (A_effective + X_effective) / (math.sqrt(2) * (B_radii + X_effective))
            return t_effective
        else:
            print('Error: Either the sum of A_ratio or X_ratio is not equal to 1.')

def A3BX2(A, B, X, A_ratio, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if (A_ratio == None and X_ratio == None):
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        A_effective = x_ratio * A_radii[0] + y_ratio * A_radii[1] + z_ratio * A_radii[2]
        X_effective = ratio * X_radii[0] + (1-ratio) * X_radii[1]

        t_effective = []
        for i in A_effective:
                t_effective.append((i + X_effective) / (math.sqrt(2) * (B_radii + X_effective)))
        df = pd.DataFrame(x_ratio, columns=['%s Ratio' % A[0]])
        df['%s Ratio' % A[1]] = y_ratio
        df['%s Ratio' % A[2]] = z_ratio
        df['A_effective'] = A_effective

        df2 = pd.DataFrame(t_effective, columns = np.round(ratio,2))

        df_merged = pd.merge(df,df2,left_index=True,right_index=True)
        df_merged = df_merged.rename(columns = {0.0 : '%s Ratio : 0.0' % X[0]})
        return df_merged

    elif ((A_ratio == None and X_ratio != None) or (A_ratio != None and X_ratio == None)):
        print('Warning: Insert a list of ratios for both A_ratio and X_ratio to calculate a specifice tolerance factor')
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        A_effective = x_ratio * A_radii[0] + y_ratio * A_radii[1] + z_ratio * A_radii[2]
        X_effective = ratio * X_radii[0] + (1-ratio) * X_radii[1]

        t_effective = []
        for i in A_effective:
                t_effective.append((i + X_effective) / (math.sqrt(2) * (B_radii + X_effective)))
        df = pd.DataFrame(x_ratio, columns=['%s Ratio' % A[0]])
        df['%s Ratio' % A[1]] = y_ratio
        df['%s Ratio' % A[2]] = z_ratio
        df['A_effective'] = A_effective

        df2 = pd.DataFrame(t_effective, columns = np.round(ratio,2))

        df_merged = pd.merge(df,df2,left_index=True,right_index=True)
        df_merged = df_merged.rename(columns = {0.0 : '%s Ratio : 0.0' % X[0]})
        return df_merged

    elif (A_ratio != None and X_ratio != None):
        if (sum(A_ratio) == 1 and sum(X_ratio) == 1):
            A_radii = []
            X_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            B_radii = B_dictionary.get(B[0])

            A_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1] + A_ratio[2] * A_radii[2]
            X_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1]
            t_effective = (A_effective + X_effective) / (math.sqrt(2) * (B_radii + X_effective))
            return t_effective
        else:
            print('Error: Either the sum of A_ratio or X_ratio is not equal to 1.')

def A2BX3(A, B, X, A_ratio, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if (A_ratio == None and X_ratio == None):
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        X_effective = x_ratio * X_radii[0] + y_ratio * X_radii[1] + z_ratio * X_radii[2]
        A_effective = ratio * A_radii[0] + (1-ratio) * A_radii[1]

        t_effective = []
        for i in X_effective:
                t_effective.append((A_effective + i) / (math.sqrt(2) * (B_radii + i)))
        df = pd.DataFrame(x_ratio, columns=['%s Ratio' % X[0]])
        df['%s Ratio' % X[1]] = y_ratio
        df['%s Ratio' % X[2]] = z_ratio
        df['X_effective'] = X_effective

        df2 = pd.DataFrame(t_effective, columns = np.round(ratio,2))

        df_merged = pd.merge(df,df2,left_index=True,right_index=True)
        df_merged = df_merged.rename(columns = {0.0 : '%s Ratio : 0.0' % A[0]})
        return df_merged

    elif ((A_ratio == None and X_ratio != None) or (A_ratio != None and X_ratio == None)):
        print('Warning: Insert a list of ratios for both A_ratio and X_ratio to calculate a specifice tolerance factor')
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        X_effective = x_ratio * X_radii[0] + y_ratio * X_radii[1] + z_ratio * X_radii[2]
        A_effective = ratio * A_radii[0] + (1-ratio) * A_radii[1]

        t_effective = []
        for i in X_effective:
                t_effective.append((A_effective + i) / (math.sqrt(2) * (B_radii + i)))
        df = pd.DataFrame(x_ratio, columns=['%s Ratio' % X[0]])
        df['%s Ratio' % X[1]] = y_ratio
        df['%s Ratio' % X[2]] = z_ratio
        df['X_effective'] = X_effective

        df2 = pd.DataFrame(t_effective, columns = np.round(ratio,2))

        df_merged = pd.merge(df,df2,left_index=True,right_index=True)
        df_merged = df_merged.rename(columns = {0.0 : '%s Ratio : 0.0' % A[0]})
        return df_merged

    elif (A_ratio != None and X_ratio != None):
        if sum(A_ratio) == 1 and sum(X_ratio) == 1:
            A_radii = []
            X_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            B_radii = B_dictionary.get(B[0])

            A_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1]
            X_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1] + X_ratio[2] * X_radii[2]
            t_effective = (A_effective + X_effective) / (math.sqrt(2) * (B_radii + X_effective))
            return t_effective
        else:
            print('Error: Either the sum of A_ratio or X_ratio is not equal to 1.')

def A3BX3(A, B, X, A_ratio, X_ratio, A_dictionary, B_dictionary, X_dictionary):
    if (A_ratio == None and X_ratio == None):
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        X_effective = x_ratio * X_radii[0] + y_ratio * X_radii[1] + z_ratio * X_radii[2]
        A_effective = x_ratio * A_radii[0] + y_ratio * A_radii[1] + z_ratio * A_radii[2]

        t_effective = np.zeros(shape=[len(X_effective), len(A_effective)])
        i_count = 0
        for i in X_effective:
            j_count = 0
            for j in A_effective:
                t_effective[i_count][j_count] = (j + i) / (math.sqrt(2) * (B_radii + i))
                j_count += 1
            i_count += 1

        X_labels = []
        A_labels = []
        for i in range(len(x_ratio)):
            X_labels.append("%s: %s, %s: %s, %s: %s" % (X[0], np.round(x_ratio[i],2), X[1], np.round(y_ratio[i],2), X[2], np.round(z_ratio[i],2)))
            A_labels.append("%s: %s, %s: %s, %s: %s" % (A[0], np.round(x_ratio[i],2), A[1], np.round(y_ratio[i],2), A[2], np.round(z_ratio[i],2)))

        df = pd.DataFrame(t_effective, columns=A_labels)
        df['X Ratio Index'] = X_labels
        df = df.set_index('X Ratio Index')
        return df

    elif ((A_ratio == None and X_ratio != None) or (A_ratio != None and X_ratio == None)):
        print('Warning: Insert a list of ratios for both A_ratio and X_ratio to calculate a specifice tolerance factor')
        A_radii = []
        X_radii = []
        for ions in A:
            A_radii.append(A_dictionary.get(ions))
        for ions in X:
            X_radii.append(X_dictionary.get(ions))
        B_radii = B_dictionary.get(B[0])

        ratio = np.linspace(0, 1, num=11)

        x_ratio = []
        y_ratio = []
        z_ratio = []

        x = np.linspace(0,1,11)
        y = np.linspace(0,1,11)
        xx, yy = np.meshgrid(x, y)
        z = -xx -yy +1
        for i in range(len(x)):
            for j in range(len(y)):

                if z[i][j] >= 0:
                    x_ratio.append(x[i])
                    y_ratio.append(y[j])
                    z_ratio.append(z[i][j])
                else:
                    continue

        x_ratio = np.array(x_ratio)
        y_ratio = np.array(y_ratio)
        z_ratio = np.array(z_ratio)

        X_effective = x_ratio * X_radii[0] + y_ratio * X_radii[1] + z_ratio * X_radii[2]
        A_effective = x_ratio * A_radii[0] + y_ratio * A_radii[1] + z_ratio * A_radii[2]

        t_effective = np.zeros(shape=[len(X_effective), len(A_effective)])
        i_count = 0
        for i in X_effective:
            j_count = 0
            for j in A_effective:
                t_effective[i_count][j_count] = (j + i) / (math.sqrt(2) * (B_radii + i))
                j_count += 1
            i_count += 1

        X_labels = []
        A_labels = []
        for i in range(len(x_ratio)):
            X_labels.append("%s: %s, %s: %s, %s: %s" % (X[0], np.round(x_ratio[i],2), X[1], np.round(y_ratio[i],2), X[2], np.round(z_ratio[i],2)))
            A_labels.append("%s: %s, %s: %s, %s: %s" % (A[0], np.round(x_ratio[i],2), A[1], np.round(y_ratio[i],2), A[2], np.round(z_ratio[i],2)))

        df = pd.DataFrame(t_effective, columns=A_labels)
        df['X Ratio Index'] = X_labels
        df = df.set_index('X Ratio Index')
        return df
    elif (A_ratio != None and X_ratio != None):
        if sum(A_ratio) == 1 and sum(X_ratio) == 1:
            A_radii = []
            X_radii = []
            for ions in A:
                A_radii.append(A_dictionary.get(ions))
            for ions in X:
                X_radii.append(X_dictionary.get(ions))
            B_radii = B_dictionary.get(B[0])

            A_effective = A_ratio[0] * A_radii[0] + A_ratio[1] * A_radii[1] + A_ratio[2] * A_radii[2]
            X_effective = X_ratio[0] * X_radii[0] + X_ratio[1] * X_radii[1] + X_ratio[2] * X_radii[2]
            t_effective = (A_effective + X_effective) / (math.sqrt(2) * (B_radii + X_effective))
            return t_effective
        else:
            print('Error: Either the sum of A_ratio or X_ratio is not equal to 1.')

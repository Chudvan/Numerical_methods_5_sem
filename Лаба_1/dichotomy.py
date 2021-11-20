from sympy import *
import numpy as np
import matplotlib.pyplot as plt


class DichotomysMethod:
    @staticmethod
    def dichotomy(a, b, F, var_, E=1e-3, i=0):
        if a >= b:
            print('ERROR: a >= b.')
            return
        i += 1
        f_a = F.subs(var_, a)
        f_b = F.subs(var_, b)
        c = (a + b) / 2
        if f_a * f_b >= 0:
            print("ERROR: Can't solve by dichotomy.")
            return
        if abs(f_a - f_b) <= 2 * E:
            return c
        print(f'x{i} = {c}')
        f_c = F.subs(var_, c)
        if f_c * f_b < 0:
            return __class__.dichotomy(c, b, F, var_, E, i)
        return __class__.dichotomy(a, c, F, var_, E, i)

    @staticmethod
    def plot(a, b, F, var):
        a_b = np.linspace(a, b, 100)
        f = lambdify(var, F, 'numpy')(a_b)
        title = '$' + 'y = ' + latex(F, mode='inline')[1:]
        plt.plot(a_b, f, label=title)
        plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0)

    @staticmethod
    def residual(F, var_, solution):
        return F.subs(var_, solution)

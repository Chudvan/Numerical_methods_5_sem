# -*- coding: utf-8 -*-
from sympy import *
import matplotlib.pyplot as plt
import numpy as np

init_printing()


class LeastSquares:
    @staticmethod
    def create_polinom(x_, x_list, y_list, n):
        # n - порядок полинома
        if len(x_list) != len(y_list):
            return
        x_buffer = [1 for _ in range(len(x_list))]
        y_buffer = [i for i in y_list]

        A = zeros(n + 1)
        b = Matrix([0 for _ in range(n + 1)])

        # последняя строка матрицы A и матрица b
        for j in range(n, -1, -1):
            b[j] = sum(y_buffer)
            A[n, j] = sum(x_buffer)
            x_buffer = [x_buffer[i] * x_list[i] for i in range(len(x_list))]
            y_buffer = [y_buffer[i] * x_list[i] for i in range(len(x_list))]

        # n первых строк матрицы A
        for i in range(n - 1, -1, -1):
            for j in range(n, 0, -1):
                A[i, j] = A[i + 1, j - 1]
            A[i, 0] = sum(x_buffer)
            x_buffer = [x_buffer[j] * x_list[j] for j in range(len(x_list))]

        print("Matrix A:")
        display(A)
        print("Matrix b:")
        display(b)

        # Решение СЛАУ
        solution = tuple(linsolve((A, b)))[0]
        coeffs = [N(i, 5) for i in solution]

        # Формирование многочлена n-ой стпени
        x_current = 1
        res = 0
        for i in range(len(coeffs) - 1, -1, -1):
            res += x_current * coeffs[i]
            x_current *= x_
        s = 0
        for i, x_i in enumerate(x_list):
            s += (y_list[i] - res.subs(x_, x_i)) ** 2
        print("Сумма квадратов отклонений между узловыми значениями функции и аппроксимирующей кривой:", s)
        return res, coeffs

    @staticmethod
    def split_intervals(x_list, accuracy):
        x_values = []
        for i in range(len(x_list) - 1):
            x_i = x_list[i]
            x_values.append(x_i)
            step = (x_list[i + 1] - x_i) / (accuracy + 1)
            for j in range(1, accuracy + 1):
                x_values.append(x_i + step * j)
        x_values.append(x_list[-1])
        return x_values

    @staticmethod
    def plot(x_list, coeffs):
        y_list = [np.polyval(coeffs, x_i) for x_i in x_list]
        plt.plot(x_list, y_list)

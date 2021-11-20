# -*- coding: utf-8 -*-
from sympy import *
init_printing()


class Interpolate:
    @staticmethod
    def build_Lagrange(x_, x_list, y_list):
        if len(x_list) != len(y_list):
            return
        res = 0
        for j, y_j in enumerate(y_list):
            r = 1
            x_j = x_list[j]
            for i, x_i in enumerate(x_list):
                if i != j:
                    r *= (x_ - x_i) / (x_j - x_i)
            res += r * y_j
        return res

    @staticmethod
    def get_y_values(x_, x_list, y_list, n, accuracy):
        # n - количество точек для интерполяции
        # для полинома 5-ой степени n = 6
        # accuracy - количество промежуточных точек между
        # узлами (степень гладкости)
        if len(x_list) != len(y_list):
            return
        number_intervals = len(x_list) - 1
        # Интерполяция происходит следующим образом:
        # несколько левых(left) интервалов первым полиномом
        # интервалы в центре - по одному на каждый новый полином
        # правые(right) интервалы - последним полиномом
        left = n // 2
        right = n - left
        x_values = []
        y_values = []
        # Параметры для запуска 3-х циклов:
        # (number_intervals, l_from, l_to, disp)
        # disp = число пройденных интервалов
        interpolate_parts = [(left, 0, n, 0),
                             *[(1, i + 1, n + 1 + i, left + i) for i in range(number_intervals - n)],
                             (right, -n, len(x_list), number_intervals - right)]
        for param in interpolate_parts:
            disp = param[3]
            l_from = param[1]
            l_to = param[2]
            lagrange = __class__.build_Lagrange(x_, x_list[l_from:l_to], y_list[l_from:l_to])
            # Итерируемся по интервалам
            for i in range(disp, param[0] + disp):
                x_values.append(x_list[i])
                y_values.append(y_list[i])
                step = (x_list[i + 1] - x_list[i]) / (accuracy + 1)
                for j in range(1, accuracy + 1):
                    x_values.append(x_list[i] + step * j)
                    y_values.append(lagrange.subs(x_, x_list[i] + step * j))
        x_values.append(x_list[-1])
        y_values.append(y_list[-1])
        return x_values, y_values

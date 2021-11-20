from sympy import *

ROUNDING = 10

def _stop(delta, E):
    m = max(map(lambda x: abs(x), delta))
    return m <= E


def _create_w(F, list_var):
    m = []
    for eq in F:
        empty = []
        for var_ in list_var:
            empty.append(eq.diff(var_))
        m.append(empty)
    return Matrix(m)


def _solve_linear_system(system, list_var):
    solution = solve_linear_system(system, *list_var)
    return Matrix([N(solution[var_], ROUNDING) for var_ in list_var])


class NewtonsMethod:
    @staticmethod
    def solve(x_0, F, *var_list, E=1e-3):
        x_k = Matrix(x_0)
        delta = []
        for i in range(1, len(var_list) + 1):
            var(f'dx{i}')
            delta.append(eval(f'dx{i}'))
        print('Матрица системы нелинейных уравнений F:')
        display(F)
        w = _create_w(F, var_list)
        print('Матрица Якоби W:')
        display(w)
        i = 0
        print(f'x{i} = {x_0}')
        w_k = w
        f_k = F
        for var_ in zip(var_list, x_k):
            w_k = w_k.subs(var_[0], var_[1])
            f_k = f_k.subs(var_[0], var_[1])
        delta_k = _solve_linear_system(w_k.row_join(-f_k), delta)
        x_k_1 = x_k + delta_k
        i += 1
        print(f'x{i} = {list(x_k_1)}')
        while not _stop(delta_k, E):
            x_k = x_k_1
            w_k = w
            f_k = F
            for var_ in zip(var_list, x_k):
                w_k = w_k.subs(var_[0], var_[1])
                f_k = f_k.subs(var_[0], var_[1])
            delta_k = _solve_linear_system(w_k.row_join(-f_k), delta)
            x_k_1 = x_k + delta_k
            i += 1
            print(f'x{i} = {list(x_k_1)}')
        return list(x_k_1)

    @staticmethod
    def residual(F, *var_list, solution=None):
        if not solution:
            print('ERROR: Give Solution!')
        residual = F
        for var_ in zip(var_list, solution):
            residual = residual.subs(var_[0], var_[1])
        return residual

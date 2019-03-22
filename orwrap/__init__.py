from ortools.linear_solver import pywraplp
try:
    from ortools_ext.matrix_constraints import MakeMatrixConstraint
except ImportError:
    print("Could not import cythonized MakeMatrixConstraint!")

    def MakeMatrixConstraint(solver, coefficients, lin_expr, lb, ub):
        for i, row in enumerate(coefficients):
            expr = solver.Sum((c * expr for c, expr in zip(row, lin_expr)))
            solver.Add(pywraplp.LinearConstraint(expr, lb[i], ub[i]))


# provide a dictionary with status strings
status_txt = {
    pywraplp.Solver.OPTIMAL: 'OPTIMAL',
    pywraplp.Solver.FEASIBLE: 'FEASIBLE',
    pywraplp.Solver.INFEASIBLE: 'INFEASIBLE',
    pywraplp.Solver.UNBOUNDED: 'UNBOUNDED',
    pywraplp.Solver.ABNORMAL: 'ABNORMAL',
    pywraplp.Solver.NOT_SOLVED: 'NOT_SOLVED'
}


class Solver(pywraplp.Solver):
    def MakeMatrixConstraint(self, coef, lin_expr, lb, ub):
        """
        Allows fast insertion of many constraints at once.
        lb <= coef * lin_expr <= ub

        :param coef: np.ndarray[rows, cols]
        :param lin_expr: list if <cols> entries containing linear expressions
        :param lb: np.ndarray[rows] with lower bound
        :param ub: np.ndarray[rows] with upper bound
        :return:
        """
        return MakeMatrixConstraint(self, coef, lin_expr, lb, ub)

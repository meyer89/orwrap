from ortools.linear_solver import pywraplp
from ortools.linear_solver.linear_solver_natural_api import Constant
try:
    from orwrap.matrix_constraints import MakeMatrixConstraint
except ImportError:
    from scipy.sparse import isspmatrix
    print("Could not import cythonized MakeMatrixConstraint!")

    def MakeMatrixConstraint(solver, coefficients, lin_expr, lb, ub):
        if isspmatrix(coefficients):
            coefficients = coefficients.toarray()

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
    # define an expression which does nothing
    ZeroExpr = Constant(0)

    def MakeMatrixConstraint(self, coef, lin_expr, lb, ub):
        """
        Allows fast insertion of many constraints at once.
        lb <= coef * lin_expr <= ub

        :param coef: np.ndarray[rows, cols]
        :param lin_expr: list if <cols> entries containing linear expressions
        :param lb: np.ndarray[rows] with lower bound
        :param ub: np.ndarray[rows] with upper bound
        """
        return MakeMatrixConstraint(self, coef, lin_expr, lb, ub)

    def NumVarDict(self, index, lb, ub, name):
        return {i: self.NumVar(lb, ub, name+"_"+str(i)) for i in index}

    def BoolVarDict(self, index, name):
        return {i: self.BoolVar(name+"_"+str(i)) for i in index}

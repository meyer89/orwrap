from time import time
import pytest
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix, issparse, hstack, eye
from orwrap import Solver, pywraplp


@pytest.mark.parametrize("coef, lb, ub", [
    (np.tri(100), np.zeros(100), np.ones(100)),
    (csr_matrix(np.tri(100)), np.zeros(100), np.ones(100)),
    (coo_matrix(np.tri(200)), np.zeros(200), np.ones(200)),
    (hstack((coo_matrix(np.tri(200)), eye(200))), np.zeros(200), np.ones(200))
])
def test_speedup_lpformat(coef, lb, ub):
    # create the model the fast way with matrix constraints
    t0 = time()
    res_fast = create_solver_fast(coef, lb, ub)
    t_fast = time() - t0

    # create the same model slow with loop over rows
    t0 = time()
    res_slow = create_solver_slow(coef, lb, ub)
    t_slow = time() - t0

    # make sure there is a significant speedup
    assert t_fast < 0.05 * t_slow
    # make sure the exported model is still identical
    assert res_fast.ExportModelAsLpFormat(obfuscated=False) == res_slow.ExportModelAsLpFormat(obfuscated=False)


def create_solver_fast(coef, lb, ub):
    solver, lin_expr = init_solver(Solver, coef.shape[1])
    solver.MakeMatrixConstraint(coef, lin_expr, lb, ub)
    return solver


def create_solver_slow(coef, lb, ub):
    if issparse(coef):
        coef = coef.toarray()

    solver, lin_expr = init_solver(pywraplp.Solver, coef.shape[1])
    for i, row in enumerate(coef):
        expr = solver.Sum((c * expr for c, expr in zip(row, lin_expr)))
        solver.Add(pywraplp.LinearConstraint(expr, lb[i], ub[i]))
    return solver


def init_solver(fcn, n_vars):
    solver = fcn('Test', Solver.GLOP_LINEAR_PROGRAMMING)
    lin_expr = [(i+1) * solver.NumVar(lb=0, ub=10, name='x' + str(i)) + i for i in range(n_vars)]
    return solver, lin_expr

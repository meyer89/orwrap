# cython: language_level=3
# distutils: language=c++
from cython cimport floating
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc
from ortools.linear_solver.linear_solver_natural_api import OFFSET_KEY, CastToLinExp
from ortools.linear_solver.pywraplp import Solver
from scipy.sparse import isspmatrix, coo_matrix
from numpy import argsort, float64


# Define the Swig struct which has the real pointer to MPVariable/MPSolver hidden inside
cdef extern from "swigpyobject.h":
    ctypedef struct SwigPyObject: 
        void *ptr


# Define the few function from linear_solver.h which are needed for constraint construction
cdef extern from "ortools/linear_solver/linear_solver.h" namespace "operations_research":        
    cdef cppclass MPVariable:
        double value_solution()

    cdef cppclass MPConstraint:
        void SetCoefficient(MPVariable*, double)
        void SetBounds(double, double)

    cdef cppclass MPSolver:
        MPConstraint* MakeRowConstraint()


# define a struct used in linear expression
ctypedef struct LinExprEntry:
    MPVariable *var
    double weight


# define simple C-structure for linear expression with variables and offset
ctypedef struct LinExpr:
    vector[LinExprEntry] entries
    double offset


# Wrapper for the internal function which converts the data
def MakeMatrixConstraint(solver, coefficients, lin_expr, double[:] lb, double[:] ub):
    # Make sure it is the right SwigPyObject
    assert isinstance(solver, Solver)
    assert coefficients.shape[0] == lb.shape[0]
    assert coefficients.shape[0] == ub.shape[0]
    assert coefficients.shape[1] == len(lin_expr)
    
    # extract MPSolver from SWIG
    cdef SwigPyObject* swig_obj = <SwigPyObject*?>solver.this
    cdef MPSolver* solver_ptr = <MPSolver*>swig_obj.ptr

    # Extract the linear expression
    cdef vector[LinExpr] vec_lin_expr
    for x in lin_expr:
        vec_lin_expr.push_back(pycoefs_to_linexpr(CastToLinExp(x).GetCoeffs()))

    if isspmatrix(coefficients):
        A = coefficients if coefficients.format == 'coo' else coo_matrix(coefficients)
        A.sum_duplicates()
        idx = argsort(A.row)
        if A.data.dtype == float64:
            return MakeMatrixConstraintSparse64(solver_ptr, A.data[idx], A.col[idx], A.row[idx], vec_lin_expr, lb, ub)
        else:
            return MakeMatrixConstraintSparse32(solver_ptr, A.data[idx], A.col[idx], A.row[idx], vec_lin_expr, lb, ub)
    else:
        # call internal function
        if coefficients.dtype == float64:
            return MakeMatrixConstraintDense64(solver_ptr, coefficients, vec_lin_expr, lb, ub)
        else:
            return MakeMatrixConstraintDense32(solver_ptr, coefficients, vec_lin_expr, lb, ub)


cdef LinExpr pycoefs_to_linexpr(coefs_in):
    cdef SwigPyObject* swig_obj
    cdef MPVariable* current_variable
    cdef double term_factor

    # initialize variables for the linear expression struct
    cdef vector[LinExprEntry] entries
    cdef double offset = coefs_in.pop(OFFSET_KEY, 0.0)

    # iterate over the variables
    for py_var, term_factor in coefs_in.items():
        swig_obj = <SwigPyObject*?>py_var.this
        current_variable = <MPVariable*?>swig_obj.ptr
        # Now we have a pointer of MPVariable and can build or small struct with the data
        if term_factor != 0:
            entries.push_back(LinExprEntry(var=current_variable, weight=term_factor))

    # return the struct with the entries
    return LinExpr(entries=entries, offset=offset)


cdef MakeMatrixConstraintSparse64(MPSolver* solver_ptr, double[:] coef_data, int[:] coef_col, int[:] coef_row, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    return MakeMatrixConstraintSparse(solver_ptr, coef_data, coef_col, coef_row, lin_expr, lb, ub)


cdef MakeMatrixConstraintSparse32(MPSolver* solver_ptr, float[:] coef_data, int[:] coef_col, int[:] coef_row, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    return MakeMatrixConstraintSparse(solver_ptr, coef_data, coef_col, coef_row, lin_expr, lb, ub)


cdef MakeMatrixConstraintSparse(MPSolver* solver_ptr, floating[:] coef_data, int[:] coef_col, int[:] coef_row, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    # Extract the sizes of the matrix
    cdef Py_ssize_t n_coef = coef_data.shape[0]

    # Variable definition with static type information
    cdef MPConstraint* current_constraint
    cdef double constraint_offset
    cdef int i_constraint, col

    # create a new constraint without any information
    current_constraint = solver_ptr.MakeRowConstraint()
    # initialize the offset to zero
    constraint_offset = 0

    # define the index of current constraint
    i_constraint = 0

    for i in range(n_coef):
        # check if there is a new constraint
        if coef_row[i] != i_constraint:
            # finalize previous constraint
            current_constraint.SetBounds(lb[i_constraint] - constraint_offset, ub[i_constraint] - constraint_offset)

            # add new constraint
            current_constraint = solver_ptr.MakeRowConstraint()
            constraint_offset = 0
            i_constraint = coef_row[i]

        # handle this in the sparse matrix
        col = coef_col[i]
        constraint_offset += coef_data[i] * lin_expr[col].offset
        AddLinExpr(current_constraint, coef_data[i], lin_expr[col])

    # finalize the last constraint
    current_constraint.SetBounds(lb[i_constraint] - constraint_offset, ub[i_constraint] - constraint_offset)

    # Assume we do not need to change anything in these constraints, so we just give a status
    return True


# Define a function which takes a dense matrix with coefficients and linear expression lists for each column
cdef MakeMatrixConstraintDense64(MPSolver* solver_ptr, double[:, :] coefficients, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    return MakeMatrixConstraintDense(solver_ptr, coefficients, lin_expr, lb, ub)


# Define a function which takes a dense matrix with coefficients and linear expression lists for each column
cdef MakeMatrixConstraintDense32(MPSolver* solver_ptr, float[:, :] coefficients, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    return MakeMatrixConstraintDense(solver_ptr, coefficients, lin_expr, lb, ub)


# Define a function which takes a dense matrix with coefficients and linear expression lists for each column
cdef MakeMatrixConstraintDense(MPSolver* solver_ptr, floating[:, :] coefficients, vector[LinExpr] lin_expr, double[:] lb, double[:] ub):
    # Extract the sizes of the matrix
    cdef Py_ssize_t rows = coefficients.shape[0]
    cdef Py_ssize_t cols = coefficients.shape[1]

    # Variable definition with static type information
    cdef MPConstraint* current_constraint
    cdef double constraint_offset, coefficient_factor

    # create one constraints for each row in the matrix
    for row in range(rows):
        # create a new constraint without any information
        current_constraint = solver_ptr.MakeRowConstraint()
        # initialize the offset to zero
        constraint_offset = 0

        # multiply the linear expression with the coefficient factor
        for col in range(cols):
            coefficient_factor = coefficients[row, col]
            # check if the calculation is actually needed
            if coefficient_factor != 0:
                # add the offset to the constraint offset
                constraint_offset += coefficient_factor * lin_expr[col].offset
                AddLinExpr(current_constraint, coefficient_factor, lin_expr[col])

        # Set the bounds with the constant offset in the expression
        current_constraint.SetBounds(lb[row] - constraint_offset, ub[row] - constraint_offset)

    # Assume we do not need to change anything in these constraints, so we just give a status
    return True


# add one single linear expression with weight to the specified constraint
cdef void AddLinExpr(MPConstraint* constraint, double weight, LinExpr expr_coef):
    cdef vector[LinExprEntry].iterator it = expr_coef.entries.begin()
    cdef LinExprEntry current_entry
    while it != expr_coef.entries.end():
        current_entry = deref(it)
        constraint.SetCoefficient(current_entry.var, current_entry.weight * weight)
        inc(it)

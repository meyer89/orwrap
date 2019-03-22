# cython: language_level=3
# distutils: language=c++
from ortools.linear_solver.linear_solver_natural_api import OFFSET_KEY
from ortools.linear_solver.pywraplp import Solver


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


# Wrapper for the internal function which converts the data
def MakeMatrixConstraint(solver, double[:,:] coefficients, lin_expr, double[:] lb, double[:] ub):
    # Make sure it is the right SwigPyObject
    assert isinstance(solver, Solver)
    assert coefficients.shape[0] == lb.shape[0]
    assert coefficients.shape[0] == ub.shape[0]
    assert coefficients.shape[1] == len(lin_expr)
    
    # extract MPSolver from SWIG
    cdef SwigPyObject* swig_obj = <SwigPyObject*?>solver.this
    cdef MPSolver* solver_ptr = <MPSolver*>swig_obj.ptr
    
    # Extract the linear expression.
    lin_coef = [x.GetCoeffs() for x in lin_expr]
    offset = [x.pop(OFFSET_KEY, 0.0) for x in lin_coef]
    
    # call internal function
    return MakeMatrixConstraintDense(solver_ptr, coefficients, lin_coef, offset, lb, ub)


# Define a function which takes a dense matrix with coefficients and linear expression lists for each column
cdef MakeMatrixConstraintDense(MPSolver* solver_ptr, double[:,:] coefficients, expr_coef, expr_offset, double[:] lb, double[:] ub):
    # Extract the sizes of the matrix
    cdef Py_ssize_t rows = coefficients.shape[0]
    cdef Py_ssize_t cols = coefficients.shape[1]

    # Variable definition with static type information
    cdef SwigPyObject* swig_obj
    cdef MPConstraint* current_constraint
    cdef MPVariable* current_variable
    cdef double term_factor, constraint_offset, coefficient_factor

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
                constraint_offset += coefficient_factor * expr_offset[col]
                
                # iterate over the variables. This still handles with Python objects and could be done faster
                for py_var, term_factor in expr_coef[col].items():
                    swig_obj = <SwigPyObject*?>py_var.this
                    current_variable = <MPVariable*?>swig_obj.ptr
                    # Now we have a pointer of MPVariable and can use the C-function for the coefficient
                    current_constraint.SetCoefficient(current_variable, term_factor * coefficient_factor)
                    
        # Set the bounds with the constant offset in the expression
        current_constraint.SetBounds(lb[row] - constraint_offset, ub[row] - constraint_offset)

    # Assume we do not need to change anything in these constraints, so we just give a status
    return True
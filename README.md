# orwrap for ortools

Implements faster insert of matrix constraints for ortools (**can be 100x faster**)!
Can be used as a drop-in replacement for the Solver in ortools.linear_solver.pywraplp

## Installation

You need to have ortools and Cython already installed and you **must** use Linux at this time.
The installation was tested on Ubuntu and works fine here.

The headers in the _include_ directory are from the ubuntu 18.10 package. 

## Usage
 
The most important added function to the solver is MakeMatrixConstraint.

    def MakeMatrixConstraint(coef, lin_expr, lb, ub):
        """
        Allows fast insertion of many constraints at once.
        lb <= coef * lin_expr <= ub
    
        :param coef: np.ndarray[rows, cols]
        :param lin_expr: list if <cols> entries containing linear expressions
        :param lb: np.ndarray[rows] with lower bound
        :param ub: np.ndarray[rows] with upper bound
        """

## Example
    
    from orwrap import Solver
    from numpy import zeros, ones, tril
    model = Solver("Test", Solver.GLOP_LINEAR_PROGRAMMING)
    model.MakeMatrixConstraint(tril(ones((10, 10))), model.NumVarDict(range(10), 0, 10, "test").values(), zeros(10), ones(10))
    print(model.ExportModelAsLpFormat(obfuscated=False))
    
using Pardiso
using SparseArrays

verbose = true

n = 4  # The number of equations.
m = 3  # The number of right-hand sides.

A = sparse([ 0. -2  3 0
            -2  4 -4 1
            -3  5  1 1
             1 -3  0 2 ])

# Generate a random collection of right-hand sides.
B = ones(n,m)

# Initialize the PARDISO internal data structures.
ps = MKLPardisoSolver()
set_nprocs!(ps,4)
if verbose
    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
end

# If we want, we could just solve the system right now.
# Pardiso.jl will automatically detect the correct matrix type,
# solve the system and free the data
X1 = solve(ps, A, B)

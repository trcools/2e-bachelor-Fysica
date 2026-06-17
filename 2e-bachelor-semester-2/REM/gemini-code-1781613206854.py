import sympy as sp

# Define symbols
gamma, beta, c = sp.symbols('gamma beta c')
Ex, Ey, Ez = sp.symbols('E_x E_y E_z')
Bx, By, Bz = sp.symbols('B_x B_y B_z')

# Define Lambda for y-boost
Lambda = sp.Matrix([
    [gamma, 0, -gamma*beta, 0],
    [0, 1, 0, 0],
    [-gamma*beta, 0, gamma, 0],
    [0, 0, 0, 1]
])

# Define F (Griffiths convention)
F = sp.Matrix([
    [0, Ex/c, Ey/c, Ez/c],
    [-Ex/c, 0, Bz, -By],
    [-Ey/c, -Bz, 0, Bx],
    [-Ez/c, By, -Bx, 0]
])

# Compute F' = Lambda * F * Lambda^T
# Since Lambda is symmetric, Lambda^T = Lambda
F_prime = Lambda * F * Lambda

# Display the elements of interest to compare with Eq (1)
print("F'_01 (should be E'_x / c):", sp.simplify(F_prime[0,1]))
print("F'_02 (should be E'_y / c):", sp.simplify(F_prime[0,2]))
print("F'_03 (should be E'_z / c):", sp.simplify(F_prime[0,3]))
print("F'_12 (should be B'_z):", sp.simplify(F_prime[1,2]))
print("F'_13 (should be -B'_y):", sp.simplify(F_prime[1,3]))
print("F'_23 (should be B'_x):", sp.simplify(F_prime[2,3]))
# Name: Roshan Yadav
# Roll No: 2311144
# Function lib Assignment 15: Shooting method and PDE.

import math


# For Question 1: Solving Boundary value problem using Shooting With RK4
def f(T):
    """
    Returns d²T/dx² = -α(Ta - T)
    where α = 0.01, Ta = 20
    """
    alpha = 0.01
    Ta = 20
    return -alpha * (Ta - T)


def RK4_2(a=0, b=10, T_a=40, T_b=200, h=0.01, v0=0):
    """
    RK4 integrator for the coupled ODEs:
    dT/dx = v
    dv/dx = f(T)
    """
    N = int((b - a) / h)
    V = []
    X = []
    T = []
    V.append(v0)
    T.append(T_a)
    X.append(a)
    dx = h
    
    for i in range(0, N):
        v = V[i]
        temp = T[i]
        x = X[i]

        k_1T = dx * v
        k_1v = dx * f(temp)

        k_2T = dx * (v + k_1v / 2)
        k_2v = dx * f(temp + k_1T / 2)

        k_3T = dx * (v + k_2v / 2)
        k_3v = dx * f(temp + k_2T / 2)

        k_4T = dx * (v + k_3v)
        k_4v = dx * f(temp + k_3T)

        new_temp = temp + (k_1T + 2*k_2T + 2*k_3T + k_4T) / 6
        new_v = v + (k_1v + 2*k_2v + 2*k_3v + k_4v) / 6
        new_x = x + dx

        T.append(new_temp)
        V.append(new_v)
        X.append(new_x)

    return X, V, T


def RK4_shooting(a=0, b=10, T_a=40, T_b=200, h=0.01, max_iter=100, tol=0.01):
    """
    Shooting method to solve boundary value problem.
    Uses Lagrange interpolation to refine the initial slope guess.
    """
    
    # Initial guesses for the slope at x=a
    zeta_l = -1.5  
    zeta_h = -0.5  
    
    # Get solutions with initial guesses
    X, V_l, T_l = RK4_2(a, b, T_a, T_b, h, zeta_l)
    X, V_h, T_h = RK4_2(a, b, T_a, T_b, h, zeta_h)
    
    # Ensure the solution is with in the bracket
    if (T_l[-1] - T_b) * (T_h[-1] - T_b) > 0:
        if T_l[-1] < T_b and T_h[-1] < T_b:
            zeta_h = 0.5
            X, V_h, T_h = RK4_2(a, b, T_a, T_b, h, zeta_h)
        elif T_l[-1] > T_b and T_h[-1] > T_b:
            zeta_l = -3.0
            X, V_l, T_l = RK4_2(a, b, T_a, T_b, h, zeta_l)
    
    # Iterate using Lagrange interpolation
    for iteration in range(max_iter):
        # Check convergence
        if abs(T_h[-1] - T_b) < tol:
            return X, V_h, T_h
        
        if abs(T_l[-1] - T_b) < tol:
            return X, V_l, T_l
        
        # Lagrange interpolation to get new guess
        zeta_new = zeta_l + (zeta_h - zeta_l) / (T_h[-1] - T_l[-1]) * (T_b - T_l[-1])
        
        # Get solution with new guess
        X, V_new, T_new = RK4_2(a, b, T_a, T_b, h, zeta_new)
        
        # Update brackets
        if T_new[-1] < T_b:
            zeta_l, T_l, V_l = zeta_new, T_new, V_new
        else:
            zeta_h, T_h, V_h = zeta_new, T_new, V_new
    
    # Return best solution
    if abs(T_h[-1] - T_b) < abs(T_l[-1] - T_b):
        return X, V_h, T_h
    else:
        return X, V_l, T_l


# For question 2: Solving Heat Equations
def g(x):
    """
    Initial condition: heated to 300°C at center, 0°C elsewhere
    """
    if abs(x - 1.0) < 0.05:  # At center (x=1 for L=2)
        return 300
    else:
        return 0


def create_matrix_A(alpha, n):
    """
    Creates the evolution matrix for explicit heat equation scheme
    """
    A = []
    for i in range(n):
        row = [0] * n
        A.append(row)
    
    # Fill the tridiagonal matrix
    for i in range(n):
        A[i][i] = 1 - 2*alpha
        if i > 0:
            A[i][i-1] = alpha
        if i < n - 1:
            A[i][i+1] = alpha
    
    return A


def matmul(A, V):
    """
    Matrix-vector multiplication
    """
    C = []
    for i in range(len(V)):
        sum_val = 0
        for j in range(len(V)):
            sum_val += A[i][j] * V[j]
        C.append(sum_val)
    
    return C


def PDE_H(L=2, dx=0.1, dt=0.005, T=2000):
    """
    Solves 1D heat equation using explicit scheme
    u_xx = u_t
    Default dt=0.005 ensures stability (α = 0.5)
    """
    N = int(L / dx) + 1  
    alpha = dt / (dx**2)
    
    # Create spatial grid
    X = []
    for i in range(N):
        X.append(i * dx)
    
    # Initial condition
    u0 = []
    for x in X:
        u0.append(g(x))
    
    U = []
    U.append(u0)
    
    A = create_matrix_A(alpha, N - 2)  
    
    # Time evolution
    for i in range(1, T):
        V_prev = U[i-1]
        
        V_interior = V_prev[1:-1]
        
        V_new_interior = matmul(A, V_interior)
        
        V_new = [0] + V_new_interior + [0]
        
        U.append(V_new)
    
    return U, X
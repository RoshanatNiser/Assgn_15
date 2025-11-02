# Name: Roshan Yadav
# Roll No: 2311144
# Assignment 15: Shooting method and PDE.

from Func_lib_for_assgn_15 import *
import matplotlib.pyplot as plt

# Question 1: Shooting Method
print("Question 1: Solving boundary value problem using Shooting Method")


X, V, T = RK4_shooting(a=0, b=10, T_a=40, T_b=200, h=0.01)



# Find where temperature is 100°C
print("Finding x where T = 100°C:")
for i in range(len(T)):
    if abs(T[i] - 100) < 0.1:
        print(f"At index {i}: x = {X[i]:.4f} m, T = {T[i]:.4f}°C\n")
        break



# Question 2: Heat Equation
print("Question 2: Solving 1D Heat Equation\n")


U, X = PDE_H(L=2, dx=0.1, dt=0.005, T=2000)

# Plot temperature profiles at different time steps
plt.figure(figsize=(10, 6))

Time= [0,5,10,50,100,200,600]
for i in Time:
    if i < len(U):
        plt.plot(X, U[i], marker='o', label=f't = {i}', markersize=4)

plt.title('Question 2: Heat Equation - Temperature Evolution')
plt.xlabel('Position x (m)')
plt.ylabel('Temperature T(x) (°C)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("Assgn_15_Question_2.png", dpi=300, bbox_inches='tight')
plt.close()

print("Plot saved as 'Assgn_15_Question_2.png\n")


#Results
"""
Question 1: Solving boundary value problem using Shooting Method
Finding x where T = 100°C:
At index 442: x = 4.4200 m, T = 99.9218°C

Question 2: Solving 1D Heat Equation

Plot saved as 'Assgn_15_Question_2.png

"""
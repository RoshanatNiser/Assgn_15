# Name: Roshan Yadav
# Roll No: 2311144
# Assignment 15: Shooting method and PDE.

from Func_lib_for_assgn_15 import *

import matplotlib.pyplot as plt

# Question 1

X,V,T=BVP(a=0,b=10,T_a=40,T_b=200,h=0.01,x0=1,p=-1.5,q=0.1)

for i in range(0,len(T)):
    if abs(T[i] - 100) < 0.1:
        print(i)

# Question 2
U,X=PDE_H(L=2,dx=0.1,dt=0.01,T=1000)

for i in range (0,1000,200):
    plt.plot(X,U[i], label=f'T ={i}'markersize=4)

plt.title('For Question 2: Heat Equation')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"Assgn_12_Question_2.png", dpi=300, bbox_inches='tight')
plt.close()
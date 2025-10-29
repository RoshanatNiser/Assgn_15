# Name: Roshan Yadav
# Roll No: 2311144
# Function lib Assignment 15: Shooting method and PDE.


import math


# For Question 1: Solving Boundary value problem using Shooting With RK4
def f(x):
    a=0.01
    t=20
    return -a(t-x)


def RK4_2(a=0,b=10,T_a=40,T_b=200,h=0.01,x0=1,v0=0):
    d=20
    c=0.01

    N=int((b-a)/h)
    V=[]
    X=[]
    T=[]
    V.append(v0)
    T.append(T_a)
    X.append(a)
    dx=h
    
    for i in range(0,N):
        v=V[i]
        T=T[i]
        x=X[i]

        k_1T=dx*v
        k_1v=dx*f(T)

        k_2T=dx*(v + (k_1v)/2)
        k_2v=dx*(f(T + (k_1T)/2))

        k_3T=dx*(v + (k_2v)/2)
        k_3v=dx*(f(T + (k_2T)/2))

        k_4T=dx*(v+k_3v)
        k_4v=dx*(f(T + k_3T))

        t= t + (k_1T + 2*k_2T + 2*k_3T + k_4T)/6
        v= v + (k_1v + 2*k_2v + 2*k_3v + k_4v)/6
        x= x + dx

        T.append(t)
        V.append(v)
        X.append(x)

    return X,V,T

def BVP(a=0,b=10,T_a=40,T_b=200,h=0.01,x0=1,p=-1.5,q=0.1,max_iter=100):
    for i in range(0,max_iter):
        X,V_h,T=RK4_2(v0=p)
        X,V_l,T=RK4_2(V0=q)

        if abs(T[-1] - T_b) < 0.001:
            return X,V_h,T
        
        if T[-1]>T_b:
            p_l=p
            
        if T[-1]<T_b:
            p_h=q
            
        else: 
            p=p + ((q-p)/(V_h[0] - V_l[0]))*(T_b - T[-1])


    



# For question 2: Solving Heat Equations
def g(x):
    if x==1:
        return 300
    else:
        return 0

def A(a,n):
    A=[]
    I=[]
    for i in range(0,n):
        for j in range(0,n):
            I.append(0)
        A.append(I)
    
    A[0][0]=1-a

    for i in range(1,n):
        A[i][i]=1-a
        A[i][i-1]=a
        A[i-1][i]=a
    
    return A

def matmul(A,V):
    C=[]
    for i in range(0,len(V)):
        a=0
        for j in range(0,len(V)):
            b=(A[i][j])*V[j]
            a=a+b
        C.append(a)
    
    return C

def PDE_H(L=2,dx=0.1,dt=0.01,T=1000):
    N=int(L/dx)
    t0=0
    X=[]
    for i in range(0,N):
        X.append(i)
    
    u0=g(X)
    U=[]
    U.append(u0)
    a=dt/dx**2

    for i in range(1,T):
        V_0=U[i]
        A=A(a,N)
        V= matmul(A,V_0)
        U.append(V)
    
    return U,X







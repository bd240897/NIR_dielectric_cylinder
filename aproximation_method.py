# -*- coding: utf-8 -*-
import scipy.special as sc
import math
import cmath
import matplotlib.pyplot as plt 
import numpy as np

# ========================================================================
##  1 - входные данные 

# Параметры среды-1
e1 = 20
mu1 = 1
k1 = 2*math.pi*math.sqrt(mu1*e1)
eta1 = 120*math.pi*math.sqrt(mu1/e1)

# Параметры среды-2
k2 = 2*math.pi
eta2 = 120*math.pi

# Параметры программы
a = 1  # радиус круга
N = 10 # количество точек на круга 
# phi_i = 90 # угол падения [град]
phi_i = math.radians(90) # угол падения [рад]

# ========================================================================
#  Построение круга
dphi = 360/N # угловой шаг на круга [град]
# dphi = 2*math.pi/N_circl # угловой шаг на круга [рад]

phi_circl = []
x = []
z = []
for i in range(N):
    # текущий угол для круга
    cur_phi_grad = i*dphi
    cur_phi_rad = math.radians(cur_phi_grad)
    
    # cохраним угол в градусах
    phi_circl.append(cur_phi_grad)
    
    # посчитаем координаты через угол в радианах 
    x.append(math.cos(cur_phi_rad))
    z.append(math.sin(cur_phi_rad))
    
# построим график круга 
def circl_graf():
    fig = plt.figure(figsize=(8., 6.)) 
    ax = fig.add_subplot(111)   
    ax.plot(x, z, label='Круг')
    ax.legend()    

circl_graf()

# 3 - cам ходя расчета 

# замены для удобства 
pi = math.pi
sqrt = math.sqrt
atan = math.atan
ln =  math.log
atan = math.atan
J = sc.jv
H = sc.hankel1
def dH(n, x):
    return n/x*H(n,x) - H(n+1,x)    
def dJ(n, x):
    return n/x*J(n,x) - J(n+1,x)

gamma = 1 #TODO
wmn = 1
umn = 1
rmn = 10
dx = 1

# списик для хранения элементов системы 
# Z1 = np.zeros((1,N),dtype=complex)
# Z2 = np.zeros((1,N),dtype=complex)
# Y1 = np.zeros((1,N),dtype=complex)
# Y2 = np.zeros((1,N),dtype=complex)


Z = np.zeros((2*N,2*N),dtype=complex)
Y = np.zeros((2*N,2*N),dtype=complex)
    
tm = tn = np.array([1,2])
nm = nn = np.array([2,3])



for m in range(N):
    xm = x[m]
    zm = z[m]
    # tn = t[n]
    # nn = n[n]
    for n in range(N):
        xn = x[n]
        zn = z[n]
        # tn = t[n]
        # nn = n[n]
        
        # локальные замены 
        umn_plus = umn+dx/2  
        umn_minus = umn-dx/2  
        r_minus = sqrt(umn_plus**2 + wmn**2)
        r_plus = sqrt(umn_minus**2 + wmn**2)
        rmn = sqrt(umn**2 + wmn**2)
            
        if n == m:
            # электрическое поле
            E1y = k1**2*dx*1j/4*(1 + 2*1j/pi*ln(gamma*k1*rmn/2))
            E2y = k2**2*dx*1j/4*(1 + 2*1j/pi*ln(gamma*k2*rmn/2))
            
            # магнитное поле
            H1u = -1/2
            H2u = 1/2
            H1w = H2w = 0
        elif abs(n-m) <= 1:
        
            # wmn = 
            # umn = 

                
            # электрчиеское поле
            E1y = k1**2*1j/4*dx - k1**2/(2*pi) *(2*wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            + umn_plus*ln(gamma*k1/2*sqrt(umn_plus**2+wmn**2)) 
            - umn_minus*ln(gamma*k1/2*sqrt(umn_minus**2+wmn**2)) - dx)
            
            E2y = k2**2*1j/4*dx - k2**2/(2*pi) *(2*wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            + umn_plus*ln(gamma*k2/2*sqrt(umn_plus**2+wmn**2)) 
            - umn_minus*ln(gamma*k2/2*sqrt(umn_minus**2+wmn**2)) - dx)  
            
            # магнитное поле 
            H1u = 1j/8*wmn*k1**2*dx + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            H2u = 1j/8*wmn*k2**2*dx + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            
            H1w = 1j/4*(H(0,k1*r_minus) - H(0,k1*r_plus))
            H2w = 1j/4*(H(0,k2*r_minus) - H(0,k2*r_plus))
        else:
            # электрчиеское поле
            E1y = k1**2*1j/4*dx*H(0,k1*rmn)
            E2y = k2**2*1j/4*dx*H(0,k2*rmn)
               
            # магнитное поле     
            H1u = 1j/4*k1*dx*(wmn/rmn*H(1,k1*rmn))
            H2u = 1j/4*k2*dx*(wmn/rmn*H(1,k2*rmn))
            
            H1w = 1j/4*(H(0,k1*r_plus) - H(0,k1*r_minus))
            H2w = 1j/4*(H(0,k2*r_plus) - H(0,k2*r_minus))
        
        # элементы системы       
        E1mn = E1y 
        E2mn = -E2y
          
        H1mn = np.dot(tm,tn)*H1u + np.dot(tm,nn)*H1w
        H2mn = -np.dot(tm,tn)*H2u - np.dot(tm,nn)*H2w
            
        
        Z[m][n] = E1mn
        Z[m][n+N] = E2mn

        Z[m+N][n] = H1mn
        Z[m+N][n+N] = H2mn
        
    #     Z1[0,n] = E1mn 
    #     Z2[0,n] = E2mn
       
    #     Y1[0,n] = H1mn
    #     Y2[0,n] = H2mn
       
    # Z[m][0] = Z1 + Z2  
    # Y[m][0] = Y1 + Y2


# # падающее поле 
# for i in range(N):
#     xn = x[n]
#     zn = z[n]
    
#     Eiy = exp(-1j*k2*(xn*cos(phi_i) + zn*sin(phi_i)))
    
#     Hix = 1/eta2*sin(phi_i) * Eiy
#     Hiy = 1/eta2*cod(phi_i) * Eiy

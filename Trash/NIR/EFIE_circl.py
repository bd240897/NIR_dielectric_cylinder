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
a = 2  # радиус круга
N = 200 # количество точек на круга
# phi_i = 90 # угол падения [град]
phi_i = math.radians(0) # угол падения [рад]
gamma = 1.781072417990

# замены для удобства 
pi = math.pi
sqrt = math.sqrt
atan = math.atan
ln =  cmath.log
atan = cmath.atan
exp = cmath.exp
cos = math.cos
sin = math.sin

J = sc.jv
H = sc.hankel1
def dH(n, x):
    return n/x*H(n,x) - H(n+1,x)    
def dJ(n, x):
    return n/x*J(n,x) - J(n+1,x)


# ========================================================================
#  Построение круга
dphi_grad = 360/N # угловой шаг на круга [град]
dphi_rad =  math.radians(dphi_grad) # угловой шаг на круга [рад]
dx = 2*a*math.sin(dphi_rad/2)
delta = dx/2 

phi_circl = []
x = []
z = []

for i in range(N):
    # текущий угол для круга
    cur_phi_grad = i*dphi_grad
    cur_phi_rad = i*dphi_rad
    
    # cохраним угол в градусах в массив (для графика)
    phi_circl.append(cur_phi_grad)
    
    # посчитаем координаты через угол в радианах 
    zi = a*sin(cur_phi_rad)
    xi = a*cos(cur_phi_rad)
    
    # упакуем их в массив
    x.append(xi)
    z.append(zi)
    
# построим график круга 
def circl_graf():
    fig = plt.figure(figsize=(8., 6.)) 
    ax = fig.add_subplot(111)   
    ax.plot(x, z, label='Круг')
    ax.legend()    

circl_graf()

# Генерация касательных
t_ = [0]*N
for i in range(N):
    if i == N-1:
        tx = (x[0] - x[i])/dx
        tz = (z[0] - z[i])/dx
    else:
        tx = (x[i+1] - x[i])/dx
        tz = (z[i+1] - z[i])/dx
    t_[i] = (tx, tz)
t_ = np.array(t_) # переведем в тип np

# Генерация нормалей
alph = 90
alhp_r = math.radians(alph)
matrix = np.array([[cos(alhp_r), -sin(alhp_r)],[sin(alhp_r), cos(alhp_r)]])
n_ = [0]*N
for i in range(N):
    n_[i] = np.dot(matrix, t_[i])
n_ = np.array(n_) # переведем в тип np    

# координаты центров отрезков
x_midle = np.zeros(N)
z_midle = np.zeros(N)
for i in range(N):
    if i == N-1:
        x_midle[i] = (x[0] + x[i])/2
        z_midle[i] = (z[0] + z[i])/2
    else:
        x_midle[i] = (x[i+1] + x[i])/2
        z_midle[i] = (z[i+1] + z[i])/2       

# ========================================================================
# 3 - cам ходя расчета 
# списик для хранения элементов системы 
Z = np.zeros((2*N,2*N),dtype=complex)
Y = np.zeros((2*N,2*N),dtype=complex)

Es = np.zeros((N,N),dtype=complex)

for m in range(N):
    xm = x_midle[m]
    zm = z_midle[m]
    tm = t_[m]
    nm = n_[m]
    for n in range(N):
        xn = x_midle[n]
        zn = z_midle[n]
        tn = t_[n]
        nn = n_[n]
              
        # Растояние между началом координат О2 и точкой наблюдений (в ГСК)        
        x_mO2 = xm - xn
        z_mO2 = zm - zn
        
        # Положим этот вектор на новую СК
        # Координаты этой точки (в ЛСК)
        umn = x_mO2*tn[0] + z_mO2*tn[1]
        wmn = x_mO2*nn[0] + z_mO2*nn[1]
        
        
        # локальные замены 
        xmn = xm - xn
        zmn = zm - zn
        r = sqrt(xmn**2 + zmn**2)
        
        umn_plus = umn + delta  
        umn_minus = umn - delta 
        r_minus = sqrt(umn_minus**2 + wmn**2)
        r_plus = sqrt(umn_plus**2 + wmn**2)
        rmn = sqrt(umn**2 + wmn**2)
            
        if n == m:
            # электрическое поле
            # E1y = k1**2*dx*1j/4*(1 + 2*1j/pi*ln(gamma*k1*delta/2))
            E2y = k2**2*dx*1j/4*(1 + 1j*2/pi*(ln(gamma*k2*delta/2)-1))

            # # магнитное поле
            # H1u = -1/2
            # H2u = 1/2
            # H1w = H2w = 0

        elif abs(n-m) <= 1:
        
            # электрчиеское поле
            # E1y = k1**2*1j/4*delta - k1**2/(2*pi) *(2*wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            # + umn_plus*ln(gamma*k1/2*sqrt(umn_plus**2+wmn**2))
            # - umn_minus*ln(gamma*k1/2*sqrt(umn_minus**2+wmn**2)) - delta)
            
            # E2y = k2**2*1j/4*delta - k2**2/(2*pi) * (2*wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            # + umn_plus*ln(gamma*k2/2*sqrt(umn_plus**2+wmn**2))
            # - umn_minus*ln(gamma*k2/2*sqrt(umn_minus**2+wmn**2)) - delta)

            E2y = k2*dx + k2*2*1j/pi * (-wmn*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            + umn_plus*ln(gamma*k2/2*sqrt(umn_plus**2+wmn**2))
            - umn_minus*ln(gamma*k2/2*sqrt(umn_minus**2+wmn**2)) - dx)
            E2y = E2y * k2 * 1j/4



            # # магнитное поле
            # H1u = 1j/8*wmn*k1**2*delta + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            # H2u = 1j/8*wmn*k2**2*delta + 1/(2*pi)*(atan(umn_plus/wmn) - atan(umn_minus/wmn))
            #
            # H1w = 1j/4*(H(0,k1*r_minus) - H(0,k1*r_plus))
            # H2w = 1j/4*(H(0,k2*r_minus) - H(0,k2*r_plus))
        else:
            # электрчиеское поле
            # E1y = k1**2*1j/4*delta*H(0,k1*rmn)
            E2y = k2**2*1j/4*2*delta*H(0,k2*rmn)
               
            # магнитное поле     
            # H1u = 1j/4*k1*delta*(wmn/rmn*H(1,k1*rmn))
            # H2u = 1j/4*k2*delta*(wmn/rmn*H(1,k2*rmn))
            #
            # H1w = 1j/4*(H(0,k1*r_plus) - H(0,k1*r_minus))
            # H2w = 1j/4*(H(0,k2*r_plus) - H(0,k2*r_minus))
        
        # элементы системы       
        # E1mn = E1y
        # E2mn = -E2y
        #
        # H1mn = np.dot(tm,tn)*H1u + np.dot(tm,nn)*H1w
        # H2mn = -np.dot(tm,tn)*H2u - np.dot(tm,nn)*H2w
        #
        # # запись в массив для дальнешего решения
        # Z[m][n] = E1mn
        # Z[m][n+N] = E2mn
        #
        # Z[m+N][n] = H1mn
        # Z[m+N][n+N] = H2mn

        Es[m,n] = E2y


E = np.zeros(2*N, dtype=complex)
Ei = np.zeros(N, dtype=complex)

# падающее поле 
for m in range(N):
    xm = x_midle[m]
    zm = z_midle[m]
    # tm = t[m]
    # nm = n[m]
    
    # cоставляющие полей
    Eiy = exp(-1j*k2*(xm*cos(phi_i) + zm*sin(phi_i)))
    # Hix = 1/eta2*sin(phi_i) * Eiy
    # Hiz = 1/eta2*cos(phi_i) * Eiy
    
    # элементы системы  
    # Ei = Eiy
    # Hi = tm[0]*Hix + tm[1]*Hiz

    # запись в массив для дальнешего решения 
    # E[m] = Ei
    # E[m+N] = Hi
    #
    Ei[m] = Eiy
   
# Расчитаем токи
# тут зашито сразу 2 тока
# I = np.linalg.solve(Z, E)
I = np.linalg.solve(Es, Ei)

# отделим один ток от другого
# j1 = I[:N]
# j2 = I[N:]

# построим эти графики 
fig = plt.figure(figsize=(8., 6.)) 
ax = fig.add_subplot(111)
ax.plot(phi_circl, abs(I), color = 'r', label='j1')
ax.legend()  
# ax2 = fig.add_subplot(212)
# ax2.plot(phi_circl, abs(j2), color = 'g', label='j2')
# ax2.legend()
plt.show()
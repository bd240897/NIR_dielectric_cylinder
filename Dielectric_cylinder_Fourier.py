# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 22:40:12 2021

@author: Дмитрий
"""

import scipy.special as sc
import math
import cmath
import matplotlib.pyplot as plt 


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
N_harm = 100 # количество гармоник
N_circl = 100 # количество точек на круга 
# phi_i = 0 # угол падения [град]
phi_i = 2*math.pi/360*0 # угол падения [рад]

# ========================================================================
#  Построение круга
# dphi = 360/N_circl # угловой шаг на круга [град]
dphi = 2*math.pi/N_circl # угловой шаг на круга [рад]

phi_circl = []
for i in range(N_circl):
    phi_circl.append(i*dphi)
 
def graf_circl():
    """график круга в полярных координатах"""
    fig = plt.figure(figsize=(8., 6.))
    ax = fig.add_subplot(111, projection='polar')
    ax.plot(phi_circl, [a]*len(phi_circl))

# graf_circl()
# ========================================================================
# 3 - вычисление самих коэфицентов 

# для удобства обозначим цилиндрические функции как 
J = sc.jv
H = sc.hankel1
def dH(n, x):
    return n/x*H(n,x) - H(n+1,x)    
def dJ(n, x):
    return n/x*J(n,x) - J(n+1,x)

# харнение значения тока в хависимости от угла наблюдеения на сечении цилиндра

# реальный токи - сумма коэф * exp
I1_real = []
I2_real = []

# коэффициенты разложения токо
j1 = []
j2 = []

for phi in phi_circl:
    
    I1 = 0
    I2 = 0
    for n in range(N_harm):
        
        An = ((-1)**n)*cmath.exp(1j*n*phi_i)       
        Z = -eta1*J(n,k1*a)*dH(n,k2*a) + eta2*H(n,k2*a)*dJ(n,k1*a)
        j1n = 4*An/(k1*H(n,k1*a)) * (dH(n,k2*a)*J(n,k2*a) - H(n,k2*a)*dJ(n,k2*a))/Z
        j2n = 4*An/(k2*J(n,k2*a)) * (-eta1/eta2*J(n,k1*a)*dJ(n,k2*a) + dJ(n,k1*a)*J(n,k2*a))/Z
    
        # запишем коэф разложения токов чтоб потом использовать
        j1.append(j1n)
        j2.append(j2n)
    
        I1 += j1n*cmath.exp(1j*n*phi)  
        I2 += j2n*cmath.exp(1j*n*phi) 
        
    I1_real.append(I1)
    I2_real.append(I2)
    
    
I1_real_abs = list(map(abs, I1_real))
I2_real_abs = list(map(abs, I2_real))

angle_for_grad = list(map(lambda x: x * 360/(2*math.pi), phi_circl))

def graf_I_liner():
    fig = plt.figure(figsize=(8., 6.)) 
    ax = fig.add_subplot(111)   
    ax.plot(angle_for_grad, I1_real_abs, label='I1')
    ax.plot(angle_for_grad, I2_real_abs, label='I2')
    ax.legend()                         

# graf_I_liner()  
    
def graf_I_liner():
    fig = plt.figure(figsize=(8., 6.)) 
    ax = fig.add_subplot(111, projection='polar')
    ax = fig.add_subplot(111)   
    ax.plot(phi_circl, I1_real_abs, label='I1')
    ax.plot(phi_circl, I2_real_abs, label='I2')
    ax.legend()
    
graf_I_liner()   

Es_OUR_list = []
Hs_OUR_list = []
Ep_OUR_list = []
Hp_OUR_list = []

# 4 - посчитаем поле методом нашего разложения на растоянии ro = 2*a

# растояние на котром будем смотреть поле
ro = 2*a

for phi in phi_circl:    
    Es = 0
    Hs = 0
    Ep = 0
    Hp= 0
       
    for n in range(N_harm):
        
        j1n = j1[n]
        j2n = j2[n]
            
        # отраженое поле
        Es += -k2*eta2*j2n*H(n,k2*ro)*J(n,k2*a)*cmath.exp(1j*n*phi)
        Hs += -k2*1j/4*j2n*dH(n,k2*ro)*J(n,k2*a)*cmath.exp(1j*n*phi)
        
        # прошедшее поле
        Ep += -k1*eta1*j1n*H(n,k1*a)*J(n,k1*ro)*cmath.exp(1j*n*phi)
        Hp += -k1*1j/4*j1n*H(n,k1*a)*dJ(n,k1*ro)*cmath.exp(1j*n*phi)

    Es_OUR_list.append(Es)
    Hs_OUR_list.append(Hs)  
    Ep_OUR_list.append(Ep) 
    Hp_OUR_list.append(Hp) 
 
Es_OUR_list_abs = list(map(abs, Es_OUR_list))
Hs_OUR_list_abs = list(map(abs, Hs_OUR_list))
Ep_OUR_list_abs = list(map(abs, Ep_OUR_list))
Hp_OUR_list_abs = list(map(abs, Hp_OUR_list))

def graf_E_polar(phi_graf: list, 
                     Es,
                     Hs,
                     Ep,
                     Hp,
                     name: str,
                     ):
    """Функция для построение графика полей 
    в ПОЛЯРНЫХ координатах"""
    
    fig = plt.figure(figsize=(8., 6.))
    ax1 = fig.add_subplot(221, projection='polar')
    ax1.plot(phi_graf, Es, color = 'r', label='Es'+name)
    ax1.legend()
    
    ax2 = fig.add_subplot(222, projection='polar')
    ax2.plot(phi_graf, Hs, color = 'g',  label='Hs'+name)
    ax2.legend()
    
    ax3 = fig.add_subplot(223, projection='polar')
    ax3.plot(phi_graf, Ep, color = 'b',  label='Ep'+name)
    ax3.legend()
    
    ax4 = fig.add_subplot(224, projection='polar')
    ax4.plot(phi_graf, Hp, color = 'y',  label='Hp'+name)
    ax4.legend()
    
# график поля для нашего разложения токов в ряд 
graf_E_polar(phi_circl,
            Es_OUR_list_abs,
            Hs_OUR_list_abs,
            Ep_OUR_list_abs,
            Hp_OUR_list_abs,
            '_OUR'
            )
# ========================================================================                                                 # 
# 4 - посчитаем поле методом никольского на растоянии ro = 2*a
    
# растояние на котром будем смотреть поле
ro = 2*a

Es_list = []
Hs_list = []
Ep_list = []
Hp_list = []

for phi in phi_circl:
    
    Es = 0
    Hs = 0
    Ep = 0
    Hp= 0
       
    for n in range(N_harm):
        
        def dH(n, x):
            return n/x*H(n,x) - H(n+1,x)
        
        def dJ(n, x):
            return n/x*J(n,x) - J(n+1,x)        
        
        An = (-1)**n
        
        Z2 = J(n,k1*a)*dH(n,k2*a) - eta2/eta1*H(n,k2*a)*dJ(n,k1*a)     
        bn = An * (dH(n,k2*a)*J(n,k2*a) - H(n,k2*a)*dJ(n,k2*a))/Z2
        cn = An * (-J(n,k1*a)*dJ(n,k2*a) + eta2/eta1*dJ(n,k1*a)*J(n,k2*a))/Z2

        # отраженое поле
        Es += cn*H(n, k2*ro)*cmath.exp(1j*n*phi)
        Hs += 1j/eta2*cn*dH(n, k2*ro)*cmath.exp(1j*n*phi)
        
        # прошедшее поле
        Ep += bn*H(n, k1*ro)*cmath.exp(1j*n*phi)
        Hp += 1j/eta1*bn*dH(n, k1*ro)*cmath.exp(1j*n*phi)

    Es_list.append(Es)
    Hs_list.append(Hs)  
    Ep_list.append(Ep) 
    Hp_list.append(Hp) 
 
Es_list_abs = list(map(abs, Es_list))
Hs_list_abs = list(map(abs, Hs_list))
Ep_list_abs = list(map(abs, Ep_list))
Hp_list_abs = list(map(abs, Hp_list))
    
# графтк для поля из Никольского
graf_E_polar(phi_circl,
            Es_list_abs,
            Hs_list_abs,
            Ep_list_abs,
            Hp_list_abs,
            '_NIC'
            )    
    
def graf_E_liner(phi_graf: list, 
                Es,
                Hs,
                Ep,
                Hp,
                name: str,
                ):
    """Функция для построение графика полей 
    в ДЕКАРТОВЫХ координатах"""
    # обычный график 
    fig = plt.figure(figsize=(8., 6.)) 
    ax = fig.add_subplot(221)   
    ax.plot(phi_graf, Es, color = 'r', label='Es'+name)
    ax.legend()  

    ax2 = fig.add_subplot(222)   
    ax2.plot(phi_graf, Hs, color = 'g', label='Hs'+name)
    ax2.legend()
    
    ax3 = fig.add_subplot(223)   
    ax3.plot(phi_graf, Ep, color = 'b', label='Ep'+name)
    ax3.legend()
    
    ax4 = fig.add_subplot(224)   
    ax4.plot(phi_graf, Hp, color = 'y', label='Hp'+name)
    ax4.legend()
    
graf_E_liner(phi_circl,
            Es_list_abs,
            Hs_list_abs,
            Ep_list_abs,
            Hp_list_abs,
            '_NIC'
            ) 
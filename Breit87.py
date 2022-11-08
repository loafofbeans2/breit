# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 23:16:37 2022

@author: benmc
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

s = 1/2
L = 0
J = L+s
I = 3/2
mu_b = 1.3996
h = 6.626*10**(-34)
g_j = 1 + (J*(J+1)+s*(s+1)-L*(L+1))/(2*J*(J+1))
g_I = 0.000999
W_hfs = 3417*(I+1/2)
mu_I = -mu_b *g_I
N = 50000
f_0 = 4
B_max = 250

def x(B):
    return(g_j-g_I)*B*mu_b/W_hfs

def first(B, m_f):
    return -mu_b*g_I*m_f*B
def second(B, m_f):
    return -mu_I * B * m_f
def bracket(B, m_f, x):
    return W_hfs*0.5*((1+(4*m_f*x))/(2*I+1))**0.5

def breit(B, m_f):
    w_pos = first(B, m_f) + second(B, m_f) + bracket(B, m_f, x(B))
    w_neg = first(B, m_f) +second(B, m_f) - bracket(B, m_f, x(B))
    return w_pos, w_neg

B = np.linspace(0, B_max, N)


fig, ax = plt.subplots(1, 2)
ax[1].set_box_aspect(1)
ax[0].set_box_aspect(1)
#ax[2].set_box_aspect(1)

def gap(B, m_f):
    return abs(0.5*(breit(B, m_f+1)[0] - breit(B, m_f)[0]))
def gap_neg(B, m_f):
    return abs(0.5*(breit(B, m_f+1)[1] - breit(B, m_f)[1]))

F = 2
for m_f in range(-F, F,1):
    print(m_f)
    ax[1].plot(B, gap(B, m_f))
F = 1
for m_f in range( -F,F):
    print(m_f)
    ax[1].plot(B, gap_neg(B, m_f))
#Defining functions to use
def rabi1(B):
    return gap_neg(B, -1)-f_0
def rabi2(B):
    return gap_neg(B, 0)-f_0
def rabi3(B):
    return gap(B, -2)-f_0
def rabi4 (B):
    return gap(B, 1)-f_0
def rabi5 (B):
    return gap(B, 0)-f_0
def rabi6 (B):
    return gap(B, 1)-f_0

gap1 = fsolve(rabi1, 1)[0]
gap2 = fsolve(rabi2, 1)[0]
gap3 = fsolve(rabi3, 1)[0]
gap4 = fsolve(rabi4, 1)[0]
gap5 = fsolve(rabi5, 1)[0]
gap6 = fsolve(rabi6, 1)[0]
gaplist = np.array([gap1, gap2, gap3, gap4, gap5, gap6])
print('The Calculated Values for 87Rb gaps is: ', gaplist)
sweepfield =np.array([-0.09948,-0.0763,-0.05312,-0.02628,-0.00066,0.0274,])

mainfield = gaplist+sweepfield

test_data = np.array([8.94, 8.96318,8.98636,9.0132,9.03882, 9.06688])
noisey = np.random.normal(gaplist, 0.25)
print('THe Observed Values Are: ', noisey)
diff = noisey - gaplist

print('The Difference is: ', diff)
error = (diff/test_data)
print('With Percentage Difference: ', error*100)
test = gap(test_data, 2)
ax[0].scatter([1,2,3,4,5,6], diff, label = 'difference')
ax[0].plot([1,2,3,4,5,6], test_data, label = 'Actual')
ax[0].plot([1,2,3,4,5,6], gaplist, label = 'Theoretical')
#print(test)
print(x(B))

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "Anastasia Lyupa"
__date__ = "$Oct 17, 2015 9:15:50 PM$"

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg  as la

PLOTS_DIR = 'plots'
DATA_DIR = 'data'

LX = 5.0
LY = LZ = 1.0
NX = 100
P_ATM = 1.0
P0 = 300 * P_ATM
SW0 = 0.35
DH = LX / NX
TAU = 2.5e-2

G = 9.8e-5
PHI = 0.4
K = 6.64e-11
NK = 3.25
RO_W=1000.0
RO_N = 850.0
BETA_W = 4.4e-5
BETA_N = 1.0e-4
MU_W = 1.0e-8
MU_N = 1.0e-7
SWR = SNR = 0.05
LAMBDA_W = K / MU_W
LAMBDA_N = K / MU_N

Pw = np.zeros(NX, dtype=np.float64)
Pw_new = np.zeros(NX, dtype=np.float64)
Pw_old = np.zeros(NX, dtype=np.float64)
Sw = np.zeros(NX, dtype=np.float64)
Sw_old = np.zeros(NX, dtype=np.float64)
Sw_tmp = np.zeros(NX, dtype=np.float64)
Sn = np.zeros(NX, dtype=np.float64)
Rw = np.zeros(NX, dtype=np.float64)
Rn = np.zeros(NX, dtype=np.float64)
Lw = np.zeros(NX - 1, dtype=np.float64)
Ln = np.zeros(NX - 1, dtype=np.float64)
Qw = np.zeros(NX, dtype=np.float64)
Qn = np.zeros(NX, dtype=np.float64)
Cw = np.zeros(NX, dtype=np.float64)
Cn = np.zeros(NX, dtype=np.float64)

def init_p_s():
    for i in range(NX):
        Pw[i] = P0 + (1.1 * P0 - P0) * i / (NX - 1)
    Pw_new[:] = Pw
    Pw_old[:] = Pw
    Sw[:] = SW0
    Sw_old[:] = Sw
    Sn[:] = 1.0 - Sw
    Qn[:] = 0.0
    Qw[:] = 0.0

def kw(sw):
    k = 0.0
    if (sw >= 1.0 - SNR):
        k = 1.0
    elif (sw <= SWR):
        k = 0.0
    else:
        sw = (sw - SWR) / (1.0 - SWR - SNR)
        k = pow(sw, 0.5) * pow((1.0 - pow((1.0 - pow(sw, NK / (NK - 1))),
        (NK - 1) / NK)), 2)
    return k

def kn(sw):
    k = 0.0
    if (sw >= 1.0 - SNR):
        k = 0.0
    elif (sw <= SWR):
        k = 1.0
    else:
        sw = (sw - SWR) / (1.0 - SWR - SNR)
        k = pow(1.0 - sw, 0.5) * pow((1.0 - pow(sw, NK / (NK - 1))),
        2 * (NK - 1) / NK)
    return k

def set_lambda():
    for i in range(NX - 1):
        if (Pw[i + 1] - Pw[i] - Rw[i] * G * DH <= 0):
            Lw[i] = 0.5 * (Rw[i + 1] + Rw[i]) * LAMBDA_W * kw(Sw[i])
        else:
            Lw[i] = 0.5 * (Rw[i + 1] + Rw[i]) * LAMBDA_W * kw(Sw[i + 1])
        if (Pw[i + 1] - Pw[i] - Rn[i] * G * DH <= 0):
            Ln[i] = 0.5 * (Rn[i + 1] + Rn[i]) * LAMBDA_N * kn(Sw[i])
        else:
            Ln[i] = 0.5 * (Rn[i + 1] + Rn[i]) * LAMBDA_N * kn(Sw[i + 1])
        
def set_rho():
        Rw[:] = RO_W * (1.0 + BETA_W * (Pw - P_ATM))
        Rn[:] = RO_N * (1.0 + BETA_N * (Pw - P_ATM))
        
def set_c():
        Cw[:] = RO_W * BETA_W / Rw
        Cn[:] = RO_N * BETA_N / Rn

def set_s(dt):
    Sw[0] = 0.7
    Sw[NX - 1] = SW0
    for i in range(1, NX - 1):
        Sw[i] = Sw[i] + (dt / (PHI * Rw[i])) * \
        ((Lw[i] * (Pw_new[i + 1] - Pw_new[i]) - \
        Lw[i - 1] * (Pw_new[i] - Pw_new[i - 1])) / (DH * DH) - \
        (G / DH) * (0.5 * (Rw[i + 1] + Rw[i]) * Lw[i] - \
        0.5 * (Rw[i] + Rw[i - 1]) * Lw[i - 1]) + Qw[i]) - \
        Sw[i] * Cw[i] * (Pw_new[i] - Pw[i])
    Sn[:] = 1.0 - Sw

def explicit3_set_s(dt):
    Sw_tmp[:] = Sw
    Sw[0] = 0.7
    Sw[NX - 1] = SW0
    for i in range(1, NX - 1):
        Sw[i] = (dt / (PHI * Rw[i] * (1.0 + TAU / dt + \
        TAU * Cw[i] * (Pw_new[i] - Pw[i]) / dt))) * \
        (((Lw[i] * (Pw_new[i + 1] - Pw_new[i]) - \
        Lw[i - 1] * (Pw_new[i] - Pw_new[i - 1])) / (DH * DH) - \
        (G / DH) * (0.5 * (Rw[i + 1] + Rw[i]) * Lw[i] - \
        0.5 * (Rw[i] + Rw[i - 1]) * Lw[i - 1]) + Qw[i]) - \
        PHI * Rw[i] * Sw[i] * Cw[i] * (Pw_new[i] - Pw[i]) / dt + \
        PHI * Rw[i] * Sw[i] / dt + TAU * PHI * Rw[i] * (2 * Sw[i] - \
        Sw_old[i] - Sw[i] * Cw[i] * (Pw_new[i] - 2 * Pw[i] + Pw_old[i]) + \
        Cw[i] * (Pw_new[i] - Pw[i]) * Sw_old[i]) / (dt * dt))
    Sn[:] = 1.0 - Sw
    Sw_old[:] = Sw_tmp

def update_p():
    global Pw
    global Pw_new
    Pw, Pw_new = Pw_new, Pw
    
def explicit3_update_p():
    global Pw_old
    global Pw
    global Pw_new
    Pw, Pw_old = Pw_old, Pw
    Pw, Pw_new = Pw_new, Pw

def impes_set_p(dt):
    A = sp.lil_matrix((NX, NX), dtype=np.float64)
    b = np.zeros(NX)
    for i in range(1, NX - 1):
        A[i, i - 1] = (Ln[i - 1] / Rn[i] + Lw[i - 1] / Rw[i]) / (DH * DH)
        A[i, i] = (-1.0) * ((Ln[i] + Ln[i - 1]) / Rn[i] + \
        (Lw[i] + Lw[i - 1]) / Rw[i]) / (DH * DH) - (1.0 / dt) * PHI * \
        (Sn[i] * Cn[i] + Sw[i] * Cw[i])
        A[i, i + 1] = (Ln[i] / Rn[i] + Lw[i] / Rw[i]) / (DH * DH)
        b[i] = (G / DH) * ((Ln[i] * 0.5 * (Rn[i + 1] + Rn[i]) - \
        Ln[i - 1] * 0.5 * (Rn[i - 1] + Rn[i])) / Rn[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) \
        - (1.0 / dt) * Pw[i] * PHI * (Sn[i] * Cn[i] + Sw[i] * Cw[i]) - \
        (Qn[i] / Rn[i] - Qw[i] / Rw[i])
    A[0, 0] = 1.0
    A[NX - 1, NX - 1] = 1.0
    b[0] = P0 * 1.1
    b[NX - 1] = P0

    SpA = sp.csr_matrix(A)
    #print SpA.todense()
    Pw_new[:] = la.spsolve(SpA, b)

def explicit2_set_p(dt):
    Pw_new[0] = P0 * 1.1
    Pw_new[NX - 1] = P0
    for i in range(1, NX - 1):
        Pw_new[i] = Pw[i] + (dt / (PHI * (Sn[i] * Cn[i] + Sw[i] * Cw[i]))) * \
        ((-1) / DH * G * ((Ln[i] * 0.5 * (Rn[i + 1] + Rn[i]) - \
        Ln[i - 1] * 0.5 * (Rn[i - 1] + Rn[i])) / Rn[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) + \
        (Qn[i] / Rn[i] - Qw[i] / Rw[i]) + \
        ((Ln[i] * (Pw[i + 1] - Pw[i]) - Ln[i - 1] * (Pw[i] - Pw[i - 1])) / Rn[i] + \
        (Lw[i] * (Pw[i + 1] - Pw[i]) - Lw[i - 1] * (Pw[i] - Pw[i - 1])) / Rw[i]) / \
        (DH * DH))
        
def explicit3_set_p(dt):
    Pw_new[0] = P0 * 1.1
    Pw_new[NX - 1] = P0    
    for i in range(1, NX - 1):
        Pw_new[i] = (dt / (PHI * ((Sn[i] * Cn[i] + Sw[i] * Cw[i]) * \
        (1.0 + TAU / dt) + TAU * (Cw[i] - Cn[i]) * (Sw[i] - Sw_old[i])))) * \
        ((-1) / DH * G * ((Ln[i] * 0.5 * (Rn[i + 1] + Rn[i]) - \
        Ln[i - 1] * 0.5 * (Rn[i - 1] + Rn[i])) / Rn[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) + \
        (Qn[i] / Rn[i] - Qw[i] / Rw[i]) + \
        ((Ln[i] * (Pw[i + 1] - Pw[i]) - Ln[i - 1] * (Pw[i] - Pw[i - 1])) / Rn[i] + \
        (Lw[i] * (Pw[i + 1] - Pw[i]) - Lw[i - 1] * (Pw[i] - Pw[i - 1])) / Rw[i]) / \
        (DH * DH) + (PHI / dt) * (Sn[i] * Cn[i] + Sw[i] * Cw[i]) * (Pw[i] + \
        (TAU / dt) * (2 * Pw[i] - Pw_old[i])) + (2 * PHI * TAU / (dt * dt)) * \
        (Cw[i] - Cn[i]) * (Sw[i] - Sw_old[i]) * Pw[i])

def print_plots(test_name, time):
    h = np.linspace(0, LX, NX)
    fig = plt.figure()
    fig.canvas.set_window_title('Test ' + test_name + ' at ' + str(time) + 's')
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel("S")
    ax1.set_title("Saturations")
    line1, = ax1.plot(h, Sw, color='blue', lw=2)
    line2, = ax1.plot(h, Sn, color='red', lw=2)
    ax1.set_ylim([0.0, 1.0])
    
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel("P, bar")
    ax2.set_title("Pressure")
    line, = ax2.plot(h, Pw, color='green', lw=2)
    ax2.set_ylim([P0 * 0.85, P0 * 1.15])

    plt.subplots_adjust(wspace=0.2,hspace=.4)
    #plt.show()
    saveto = PLOTS_DIR + '/' + test_name + '_' + str(time) + 's.png'
    fig.savefig(saveto, format='png')
    plt.close(fig)

def save_data(test_name, time):
    datafile = open(DATA_DIR + '/' + test_name + '_' + str(time) +'s.out', 'w')
    np.savetxt(datafile, Pw, header='Pw')
    np.savetxt(datafile, Sw, header='Sw')
    datafile.close()

def get_eps(method1, method2, time):
    Pw_tmp1 = np.zeros(NX)
    Pw_tmp2 = np.zeros(NX)
    Sw_tmp1 = np.zeros(NX)
    Sw_tmp2 = np.zeros(NX)
    
    datafile = open(DATA_DIR + '/' + method1 + '_' + str(time) +'s.out', 'r')
    Data_tmp = np.loadtxt(datafile, dtype=np.float64)
    Pw_tmp1[:] = Data_tmp[:NX]
    Sw_tmp1[:] = Data_tmp[NX:]
    datafile.close()
    
    datafile = open(DATA_DIR + '/' + method2 + '_' + str(time) +'s.out', 'r')
    Data_tmp = np.loadtxt(datafile, dtype=np.float64)
    Pw_tmp2[:] = Data_tmp[:NX]
    Sw_tmp2[:] = Data_tmp[NX:]
    datafile.close()
    
    Pw_diff = np.linalg.norm(Pw_tmp1 - Pw_tmp2)
    Sw_diff = np.linalg.norm(Sw_tmp1 - Sw_tmp2)

    Pw_eps = Pw_diff / np.mean(Pw_tmp1) / NX
    Sw_eps = Sw_diff / np.mean(Sw_tmp1) / NX
    
    print 'Pressure error: ' + str(Pw_eps * 100) + '%'
    print 'Saturation error: ' + str(Sw_eps * 100) + '%'
    
def time_step(method, dt):
    set_rho()
    set_c()
    set_lambda()
    if (method == 'impes'):
        impes_set_p(dt)
        set_s(dt)
        update_p()
    elif (method == 'explicit2'):
        explicit2_set_p(dt)
        set_s(dt)
        update_p()
    elif (method == 'explicit3'):
        explicit3_set_p(dt)
        explicit3_set_s(dt)
        explicit3_update_p()

def solve_system(method, timesteps, print_t, save_plots_t, save_data_t, dt):
    init_p_s()
    # the first step
    if (method == 'impes' or method == 'explicit2'):
        time_step(method, dt)
    elif (method == 'explicit3'):
        time_step('explicit2', dt)
    # other steps
    for t in range(2, timesteps + 1):
        time_step(method, dt)
        if (t % print_t == 0):
            print('time = ' + str(t * dt) + 's')
        if (t % save_plots_t == 0):
            print_plots(method, t * dt)
        if (t % save_data_t == 0):
            save_data(method, t * dt)
    
def hello():
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    
    #solve_system('impes', 1000000, 1000, 50000, 50000, 1.0e-3)
    #solve_system('explicit2', 10000, 1000, 100000, 10000, 1.0e-5)
    #solve_system('explicit3', 100000, 1000, 100000, 10000, 1.0e-3)
    #get_eps('impes', 'explicit2', 0.1)
    #get_eps('explicit3', 'explicit2', 0.1)
    get_eps('explicit3', 'impes', 100.0)
    
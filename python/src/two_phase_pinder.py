# To change this license header, choose License Headers ig Project Properties.
# To change this template file, choose Tools | Templates
# and opeg the template ig the editor.

__author__ = "Anastasia Lyupa"
__date__ = "$Oct 17, 2015 9:15:50 PM$"

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg  as la

PLOTS_DIR = 'plots'
DATA_DIR = 'data'

LX = 67.0
LY = LZ = 1.0
NX = 100
P_ATM = 1000000.0
P0 = P_ATM * 1.3
#P0 = P_ATM + RO_W * G * 31.5
SW0 = 0.985
DH = LX / NX
TAU = 2.5e-2

G = 981.0
PHI = 0.37
K = 5.0e-7
NK = 10.0
RO_W=0.9982
RO_G = 0.000129
BETA_W = 4.4e-8
BETA_G = 1.0e-5
MU_W = 0.01
MU_G = 0.0015
SWR = 0.12
SGR = 0.0
LAMBDA_W = K / MU_W
LAMBDA_G = K / MU_G
ALPHA = 1.0e-5

P = np.zeros(NX, dtype=np.float64)
P_new = np.zeros(NX, dtype=np.float64)
P_old = np.zeros(NX, dtype=np.float64)
Pw = np.zeros(NX, dtype=np.float64)
Pw_new = np.zeros(NX, dtype=np.float64)
Pw_old = np.zeros(NX, dtype=np.float64)
Pg = np.zeros(NX, dtype=np.float64)
Pc = np.zeros(NX, dtype=np.float64)
Pc_old = np.zeros(NX, dtype=np.float64)
Sw = np.zeros(NX, dtype=np.float64)
Sw_old = np.zeros(NX, dtype=np.float64)
Sw_tmp = np.zeros(NX, dtype=np.float64)
Sg = np.zeros(NX, dtype=np.float64)
Rw = np.zeros(NX, dtype=np.float64)
Rg = np.zeros(NX, dtype=np.float64)
Lw = np.zeros(NX - 1, dtype=np.float64)
Lg = np.zeros(NX - 1, dtype=np.float64)
Qw = np.zeros(NX, dtype=np.float64)
Qg = np.zeros(NX, dtype=np.float64)
Cw = np.zeros(NX, dtype=np.float64)
Cg = np.zeros(NX, dtype=np.float64)

def init_p_s():
    for i in range(NX):
        P[i] = P_ATM + 0.3 * P_ATM * (1.0 - i / (NX - 1))
    P_new[:] = P
    P_old[:] = P
    Sw[:] = SW0
    Sw_old[:] = Sw
    Sg[:] = 1.0 - Sw
    Qg[:] = 0.0
    Qw[:] = 0.0
    for i in range(NX):
        Pc[i] = pc(Sw[i])
        Pw[i] = P[i] - Pc[i] * 0.5
        Pg[i] = P[i] + Pc[i] * 0.5
    Pw_new[:] = Pw
    Pw_old[:] = Pw

def kw(sw):
    k = 0.0
    if (sw >= 1.0 - SGR):
        k = 1.0
    elif (sw <= SWR):
        k = 0.0
    else:
        sw = (sw - SWR) / (1.0 - SWR - SGR)
        k = pow(sw, 0.5) * pow((1.0 - pow((1.0 - pow(sw, NK / (NK - 1))),
        (NK - 1) / NK)), 2)
    return k

def kg(sw):
    k = 0.0
    if (sw >= 1.0 - SGR):
        k = 0.0
    elif (sw <= SWR):
        k = 1.0
    else:
        sw = (sw - SWR) / (1.0 - SWR - SGR)
        k = pow(1.0 - sw, 0.5) * pow((1.0 - pow(sw, NK / (NK - 1))),
        2 * (NK - 1) / NK)
    return k

def pcgw(sw):
    p = 1.0 / ALPHA * pow(pow(sw, NK / (1.0 - NK)) - 1.0,1.0 / NK)
    return p

def pcgws(sw):
    ps = 1.0 / ALPHA * 1.0 / (1.0 - NK) * pow(pow(sw, NK / (1.0 - NK)) - 1.0,
    (1.0 / NK - 1.0) * pow(sw, (NK / (1.0 - NK) - 1.0)) / (1.0 - SWR - SGR))
    return ps

def pc(sw):
    sw1 = 0.1
    sw2 = 0.9
    p = 0.0
    sw = (sw - SWR) / (1.0 - SWR - SGR)
    if (sw <= sw1):
        p = pcs(sw1)* (sw - sw1) + pcgw(sw1)
    elif (sw2 < sw and sw < 1.0):
        p = pcs(sw2) * (sw - sw2) + pcgw(sw2)
    elif (sw >= 1.0):
        p = 0.0
    else:
        p = pcgw(sw)
    return p

def pcs(sw):
    sw1 = 0.1
    sw2 = 0.9
    ps = 0
    sw = (sw - SWR) / (1.0 - SWR - SGR)
    if (sw <= sw1):
        ps = pcgws(sw1)
    elif (sw2 <= sw):
        ps = pcgws(sw2)
    else:
        ps = pcgws(sw)
    return ps


def set_lambda():
    for i in range(NX - 1):
        if (Pw[i + 1] - Pw[i] - Rw[i] * G * DH <= 0):
            Lw[i] = 0.5 * (Rw[i + 1] + Rw[i]) * LAMBDA_W * kw(Sw[i])
        else:
            Lw[i] = 0.5 * (Rw[i + 1] + Rw[i]) * LAMBDA_W * kw(Sw[i + 1])
        if (Pg[i + 1] - Pg[i] - Rg[i] * G * DH <= 0):
            Lg[i] = 0.5 * (Rg[i + 1] + Rg[i]) * LAMBDA_G * kg(Sw[i])
        else:
            Lg[i] = 0.5 * (Rg[i + 1] + Rg[i]) * LAMBDA_G * kg(Sw[i + 1])
        
def set_rho():
    Rw[:] = RO_W * (1.0 + BETA_W * (Pw - P_ATM))
    Rg[:] = RO_G * (1.0 + BETA_G * (Pg - P_ATM))

def hw(sw):
    h = LAMBDA_G * LAMBDA_W * kg(sw) * kw(sw) / \
    (LAMBDA_G * kg(sw) + LAMBDA_W * kw(sw)) * pcs(sw)
    return h
def fw(sw):
    f = LAMBDA_W * kw(sw) / (LAMBDA_G * kg(sw) + LAMBDA_W * kw(sw))
    return f

def set_c():
    Cw[:] = RO_W * BETA_W / Rw
    Cg[:] = RO_G * BETA_G / Rg

def set_s(dt):
    for i in range(1, NX - 1):
        Sw[i] = Sw[i] + (dt / (PHI * Rw[i])) * \
        ((Lw[i] * (Pw_new[i + 1] - Pw_new[i]) - \
        Lw[i - 1] * (Pw_new[i] - Pw_new[i - 1])) / (DH * DH) - \
        (G / DH) * (0.5 * (Rw[i + 1] + Rw[i]) * Lw[i] - \
        0.5 * (Rw[i] + Rw[i - 1]) * Lw[i - 1]) + Qw[i]) - \
        Sw[i] * Cw[i] * (Pw_new[i] - Pw[i])
    Sw[0] = Sw[1]
    Sw[NX - 1] = Sw[NX - 2]
    Sg[:] = 1.0 - Sw
    
def set_pc():
    Pc_old[:] = Pc
    for i in range(NX):
        Pc[i] = pc(Sw[i])

def set_pw_pg():
    for i in range(0, NX):   
        Pw[i] = P[i] - Pc[i] * 0.5
        Pg[i] = P[i] + Pc[i] * 0.5

def explicit3_set_s(dt):
    Sw_tmp[:] = Sw
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
    Sw[0] = Sw[1]
    Sw[NX - 1] = Sw[NX - 2]
    Sg[:] = 1.0 - Sw
    Sw_old[:] = Sw_tmp

def update_p():
    global P, P_new, Pw, Pw_new
    P, P_new = P_new, P
    Pw, Pw_new = Pw_new, Pw

def explicit3_update_p():
    global P
    global P_old
    global P_new
    P, P_old = P_old, P
    P, P_new = P_new, P

def impes_set_p(dt):
    A = sp.lil_matrix((NX, NX), dtype=np.float64)
    b = np.zeros(NX)
    for i in range(1, NX - 1):
        A[i, i - 1] = (Lg[i - 1] / Rg[i] + Lw[i - 1] / Rw[i]) / (DH * DH)
        A[i, i] = (-1.0) * ((Lg[i] + Lg[i - 1]) / Rg[i] + \
        (Lw[i] + Lw[i - 1]) / Rw[i]) / (DH * DH) - (1.0 / dt) * PHI * \
        (Sg[i] * Cg[i] + Sw[i] * Cw[i])
        A[i, i + 1] = (Lg[i] / Rg[i] + Lw[i] / Rw[i]) / (DH * DH)
        b[i] = (G / DH) * ((Lg[i] * 0.5 * (Rg[i + 1] + Rg[i]) - \
        Lg[i - 1] * 0.5 * (Rg[i - 1] + Rg[i])) / Rg[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) \
        - (1.0 / dt) * P[i] * PHI * (Sg[i] * Cg[i] + Sw[i] * Cw[i]) - \
        (Qg[i] / Rg[i] + Qw[i] / Rw[i]) - \
        0.5 * ((Lg[i] * (Pc[i + 1] - Pc[i]) - \
        Lg[i - 1] * (Pc[i] - Pc[i - 1])) / Rg[i] + \
        (Lw[i] * (Pc[i + 1] - Pc[i]) - \
        Lw[i - 1] * (Pc[i] - Pc[i - 1])) / Rw[i]) / (DH * DH) + \
        (0.5 / dt) * PHI * (Sg[i] * Cg[i] - Sw[i] * Cw[i]) * (Pc[i] - Pc_old[i])
    A[0, 0] = 1.0
    A[NX - 1, NX - 1] = 1.0
    b[0] = 1.3 * P_ATM
    b[NX - 1] = P_ATM

    SpA = sp.csr_matrix(A)
    P_new[:] = la.spsolve(SpA, b)

def explicit2_set_p(dt):
    Pw_new[0] = P0
    Pw_new[NX - 1] = P_ATM
    for i in range(1, NX - 1):
        Pw_new[i] = Pw[i] + (dt / (PHI * (Sg[i] * Cg[i] + Sw[i] * Cw[i]))) * \
        ((-1) / DH * G * ((Lg[i] * 0.5 * (Rg[i + 1] + Rg[i]) - \
        Lg[i - 1] * 0.5 * (Rg[i - 1] + Rg[i])) / Rg[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) + \
        (Qn[i] / Rg[i] + Qw[i] / Rw[i]) + \
        ((Lg[i] * (Pw[i + 1] - Pw[i]) - Lg[i - 1] * (Pw[i] - Pw[i - 1])) / Rg[i] + \
        (Lw[i] * (Pw[i + 1] - Pw[i]) - Lw[i - 1] * (Pw[i] - Pw[i - 1])) / Rw[i]) / \
        (DH * DH))
        
def explicit3_set_p(dt):
    Pw_new[0] = P0
    Pw_new[NX - 1] = P_ATM    
    for i in range(1, NX - 1):
        Pw_new[i] = (dt / (PHI * ((Sg[i] * Cg[i] + Sw[i] * Cw[i]) * \
        (1.0 + TAU / dt) + TAU * (Cw[i] - Cg[i]) * (Sw[i] - Sw_old[i])))) * \
        ((-1) / DH * G * ((Lg[i] * 0.5 * (Rg[i + 1] + Rg[i]) - \
        Lg[i - 1] * 0.5 * (Rg[i - 1] + Rg[i])) / Rg[i] + \
        ((Lw[i] * 0.5 * (Rw[i + 1] + Rw[i]) - \
        Lw[i - 1] * 0.5 * (Rw[i - 1] + Rw[i])) / Rw[i])) + \
        (Qn[i] / Rg[i] + Qw[i] / Rw[i]) + \
        ((Lg[i] * (Pw[i + 1] - Pw[i]) - Lg[i - 1] * (Pw[i] - Pw[i - 1])) / Rg[i] + \
        (Lw[i] * (Pw[i + 1] - Pw[i]) - Lw[i - 1] * (Pw[i] - Pw[i - 1])) / Rw[i]) / \
        (DH * DH) + (PHI / dt) * (Sg[i] * Cg[i] + Sw[i] * Cw[i]) * (Pw[i] + \
        (TAU / dt) * (2 * Pw[i] - Pw_old[i])) + (2 * PHI * TAU / (dt * dt)) * \
        (Cw[i] - Cg[i]) * (Sw[i] - Sw_old[i]) * Pw[i])

def print_plots(test_name, time):
    h = np.linspace(0, LX, NX)
    fig = plt.figure()
    fig.canvas.set_window_title('Test ' + test_name + ' at ' + str(time) + 's')
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel("S")
    ax1.set_title("Saturations")
    line1, = ax1.plot(h, Sw, color='blue', lw=2)
    line2, = ax1.plot(h, Sg, color='red', lw=2)
    ax1.set_ylim([0.0, 1.0])
    
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel("P, bar")
    ax2.set_title("Pressure")
    line1, = ax2.plot(h, Pw, color='green', lw=2)
    line2, = ax2.plot(h, Pg, color='yellow', lw=2)
    line3, = ax2.plot(h, P, color='brown', lw=2)
    ax2.set_ylim([P_ATM * 0.85, P_ATM * 1.5])

    plt.subplots_adjust(wspace=0.2,hspace=.4)
    #plt.show()
    saveto = PLOTS_DIR + '/' + test_name + '_' + str(time) + 's.png'
    fig.savefig(saveto, format='png')
    plt.close(fig)

def save_data(test_name, time):
    datafile = open(DATA_DIR + '/' + test_name + '_' + str(time) +'s.out', 'w')
    np.savetxt(datafile, P, header='P')
    np.savetxt(datafile, Sw, header='Sw')
    datafile.close()

def get_eps(method1, method2, time):
    P_tmp1 = np.zeros(NX)
    P_tmp2 = np.zeros(NX)
    Sw_tmp1 = np.zeros(NX)
    Sw_tmp2 = np.zeros(NX)

    datafile = open(DATA_DIR + '/' + method1 + '_' + str(time) +'s.out', 'r')
    Data_tmp = np.loadtxt(datafile, dtype=np.float64)
    P_tmp1[:] = Data_tmp[:NX]
    Sw_tmp1[:] = Data_tmp[NX:]
    datafile.close()

    datafile = open(DATA_DIR + '/' + method2 + '_' + str(time) +'s.out', 'r')
    Data_tmp = np.loadtxt(datafile, dtype=np.float64)
    P_tmp2[:] = Data_tmp[:NX]
    Sw_tmp2[:] = Data_tmp[NX:]
    datafile.close()
    
    P_diff = np.linalg.norm((P_tmp1 - P_tmp2) / P_tmp1)
    Sw_diff = np.linalg.norm((Sw_tmp1 - Sw_tmp2) / Sw_tmp1)

    P_eps = P_diff / NX
    Sw_eps = Sw_diff / NX
    
    print 'Pressure error: ' + str(P_eps * 100) + '%'
    print 'Saturation error: ' + str(Sw_eps * 100) + '%'

def time_step(method, dt):
    set_pc()
    set_rho()
    set_c()
    set_lambda()
    if (method == 'impes'):
        impes_set_p(dt)
        set_pw_pg()
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
    else:
        print('error: wrong method')
        return

def solve_system(method, timesteps, print_t, save_plots_t, save_data_t, dt):
    init_p_s()
    # the first step
    if (method == 'impes' or method == 'explicit2'):
        time_step(method, dt)
    elif (method == 'explicit3'):
        time_step('explicit2', dt)
    else:
        print('error: wrong method')
        return
    # other steps
    for t in range(2, timesteps + 1):
        time_step(method, dt)
        if (t % print_t == 0):
            print('time = ' + str(t * dt) + 's')
        if (t % save_plots_t == 0):
            print_plots(method, t * dt)
        if (t % save_data_t == 0):
            save_data(method, t * dt)

def print_pc():
    PC = np.zeros(NX, dtype=np.float64)
    x = 0.0
    for i in range(0, NX):
        x += 1.0 / NX
        PC[i] = pc(x)
        
    h = np.linspace(0, 1, NX)
    fig = plt.figure()
    fig.canvas.set_window_title('pc()')
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel("pc")
    ax1.set_title("pc")
    line, = ax1.plot(h, PC, color='green', lw=2)
    ax1.set_ylim([0.0, 0.2 * P_ATM])
    plt.subplots_adjust(wspace=0.2,hspace=.4)
    #plt.show()
    saveto = PLOTS_DIR + '/' + 'pc.png'
    fig.savefig(saveto, format='png')
    plt.close(fig)

def solve_two_phase_pinder_test():
    if not os.path.exists(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
        
    print_pc()
    solve_system('impes', 10000000, 1000, 1000, 1000, 1.0e-3)

    
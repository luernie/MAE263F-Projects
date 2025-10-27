import numpy as np
import matplotlib.pyplot as plt
from utils.plotting import plot_rod, plot_results
from physics.forces import getFb, getFs, objfun

def run_simulation():
    # === Simulation parameters ===
    nv = 11
    ndof = 2 * nv
    midNode = nv // 2 + 1
    dt = 0.01
    RodLength = 0.1
    deltaL = RodLength / (nv - 1)
    totalTime = 50
    plotStep = 250
    saveImage = 0
    maximum_iter = 1000

    # === Material properties ===
    r0 = 1e-3
    Y = 1e9
    visc = 1000.0
    rho_metal = 7000
    rho_gl = 1000
    rho = rho_metal - rho_gl

    EI = Y * np.pi * r0**4 / 4
    EA = Y * np.pi * r0**2
    tol = EI / RodLength ** 2 * 1e-3

    # === Geometry setup ===
    nodes = np.zeros((nv, 2))
    for c in range(nv):
        nodes[c, 0] = c * deltaL
        nodes[c, 1] = 0.0

    # === Radii ===
    R = np.full(nv, deltaL / 10)
    R[midNode - 1] = 0.025

    # === Mass and damping matrices ===
    m = np.zeros(2 * nv)
    for k in range(nv):
        m_node = 4/3 * np.pi * R[k]**3 * rho_metal
        m[2*k:2*k+2] = m_node
    mMat = np.diag(m)

    C = np.zeros((2*nv, 2*nv))
    for k in range(nv):
        C[2*k, 2*k] = 6*np.pi*visc*R[k]
        C[2*k+1, 2*k+1] = 6*np.pi*visc*R[k]

    # === Gravity ===
    g = np.array([0, -9.8])
    W = np.zeros(2 * nv)
    for k in range(nv):
        W[2*k:2*k+2] = 4/3 * np.pi * R[k]**3 * rho * g

    # === Initial Conditions ===
    q0 = nodes.flatten()
    u0 = np.zeros(2 * nv)
    all_DOFs = np.arange(ndof)
    fixed_index = np.array([0, 1, 2, 3])
    free_index = np.setdiff1d(all_DOFs, fixed_index)

    Nsteps = round(totalTime / dt)
    ctime = 0
    all_pos = np.zeros(Nsteps)
    all_vel = np.zeros(Nsteps)
    mid_angle = np.zeros(Nsteps)

    for timeStep in range(1, Nsteps):
        q_new, error = objfun(q0, u0, dt, tol, maximum_iter, m, mMat, EI, EA, W, C, deltaL, free_index)
        if error < 0:
            print("Could not converge.")
            break

        u_new = (q_new - q0) / dt
        ctime += dt

        all_pos[timeStep] = q_new[2*midNode-1]
        all_vel[timeStep] = u_new[2*midNode-1]

        q0 = q_new.copy()
        u0 = u_new.copy()

        if timeStep % plotStep == 0:
            plot_rod(q_new, ctime)

    plot_results(totalTime, all_pos, all_vel, mid_angle, saveImage)
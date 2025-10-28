import numpy as np
import matplotlib.pyplot as plt
from utils.plotting import plot_rod, plot_results
from physics.forces import getFb, getFs, objfun

def run_simulation():
    # === Simulation parameters ===
    nv = 50 #11
    ndof = 2 * nv
    midNode = nv // 2 + 1
    dt = 0.01
    RodLength = 1 #0.1 before
    deltaL = RodLength / (nv - 1)
    totalTime = 50
    plotStep = 250
    saveImage = 0
    maximum_iter = 1000

    # === Material properties ===
    R0 = 0.013 # outer rod radius
    r0 = 0.011 # inner rod radius 1e-3
    Y = 70e9 # elastic modulus
    # visc = 1000.0 #No viscous force
    rho_metal = 2700 # density of aluminum
    # rho_metal = 7000 # No need for these if not under a liquid
    # rho_gl = 1000
    # rho = rho_metal - rho_gl

    EI = Y * np.pi * (R0**4 - r0**4) / 12 # bending stiffness shelled sphere
    print(EI)
    EA = Y * np.pi * (R0**2 - r0**2) # axial stiffness
    tol = EI / RodLength ** 2 * 1e-3 # tolerance number

    # === Geometry setup ===
    nodes = np.zeros((nv, 2))
    for c in range(nv):
        nodes[c, 0] = c * deltaL #Node x location
        nodes[c, 1] = 0.0 #Node Y location

    # # === Radii ===
    R = np.full(nv, R0) #Array for all the outer radius
    r = np.full(nv, r0) #inner radius
    # R = np.full(nv, deltaL / 10)
    # R[midNode - 1] = 0.025 #specify radius of middle node

    # === Mass and damping matrices ===
    m = np.zeros(2 * nv) #initialize mass matrix
    for k in range(nv):
        m_node = (np.pi * (R[k]**2 - r[k]**2) * RodLength * rho_metal) / (nv - 1)
        # m_node = 4/3 * np.pi * R[k]**3 * rho_metal # rho_metal # used to separate out 
        m[2*k:2*k+2] = m_node
    mMat = np.diag(m)

    # No Viscosity in this case
    C=0
    # C = np.zeros((2*nv, 2*nv))
    # for k in range(nv):
    #     C[2*k, 2*k] = 6*np.pi*visc*R[k]
    #     C[2*k+1, 2*k+1] = 6*np.pi*visc*R[k]

    # === Gravity ===
    g = np.array([0, -9.8])
    W = np.zeros(2 * nv)
    # for k in range(nv):
    #     W[2*k:2*k+2] = 4/3 * np.pi * R[k]**3 * rho_metal * g

    # === External Load === %TODO: Add in the P value, 
    P = 20000 # Applied Force
    Ploc = 0.75 # Location of the Force

    # === Initial Conditions ===
    q0 = nodes.flatten()
    u0 = np.zeros(2 * nv)
    all_DOFs = np.arange(ndof)
    # fixed_index = np.array([0, 1, 2, 3])
    fixed_index = np.array([0, 1, nv*2-1]) # There is a pin on the first node and a roller on the last node
    free_index = np.setdiff1d(all_DOFs, fixed_index)

    Nsteps = round(totalTime / dt)
    ctime = 0
    all_pos = np.zeros(Nsteps)
    all_vel = np.zeros(Nsteps)
    mid_angle = np.zeros(Nsteps)

    for timeStep in range(1, Nsteps):
        q_new, error = objfun(q0, u0, dt, tol, maximum_iter, m, mMat, EI, EA, W, C, deltaL, free_index, P, Ploc, nodes)
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

    #TODO: Work on the P force value 2000N 0.75 away and take away gravity

    # Get it so that 700N on the 75% percent nodes

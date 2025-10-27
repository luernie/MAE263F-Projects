import numpy as np
from physics.bending import gradEb, hessEb
from physics.stretching import gradEs, hessEs

def getFs(q, EA, deltaL):
    # Compute stretching force
    pass

def getFb(q, EI, deltaL):
    # Compute bending force
    pass

def objfun(q_old, u_old, dt, tol, maximum_iter, m, mMat, EI, EA, W, C, deltaL, free_index):
    # Nonlinear solver for equilibrium update
    pass

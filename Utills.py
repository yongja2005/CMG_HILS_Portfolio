import numpy as np
import math
from scipy.spatial.transform import Rotation
from scipy.integrate import solve_ivp


def degtorad(*args):
    result = np.zeros( shape=(1,len(args)) )
    j = 0
    for i in args:
        result[0,j] = math.radians(i)
        j += 1
    return np.array(result)

def transf(elem):
    # transfomation matrix
    if len(elem) == 3: # len(elem) == 3  => euler to quaternion
        q_cur_rot = Rotation.from_euler('xyz', [elem], degrees=False)
        return q_cur_rot.as_quat()[:, ::-1] # reverse order
    else: # len(elem) == 4  => quaternion to euler
        euler_rot = Rotation.from_quat([elem])
        return euler_rot[:,::-1].as_euler('xyz', degrees=True)

def skew(elem):
    # skew-symmetric matrix
    if len(elem) == 3: # len(elem) == 3 => skew-symmetric matrix of w, Kinematics Equation 에서 사용됨
        return np.array([
            [ 0, elem[2], -elem[1], elem[0] ],
            [ -elem[2], 0, elem[0], elem[1] ],
            [ elem[1], -elem[0], 0, elem[2] ],
            [ -elem[0], -elem[1], -elem[2], 0 ]
            ])
    else: # len(elem) == 4 => skew-symmetric matrix of q_cmd, 쿼터니언 오차 구할 때 사용됨
        return np.array([
            [elem[3], elem[2], -elem[1], -elem[0]],
            [-elem[2], elem[3], elem[0], -elem[1]],
            [elem[1], -elem[0], elem[3], -elem[2]],
            [elem[0], elem[1], elem[2], elem[3]]
            ])

def AngMomCMG(gamma, beta, h_a): # Angular Momentum of cmg
    # A : matrix of configuration of cmg
    A = np.array([
    [-math.cos(beta)*math.sin(gamma[0]), -math.cos(gamma[1]), math.cos(beta)*math.sin(gamma[2]), math.cos(gamma[3])],
    [math.cos(gamma[0]), -math.cos(beta)*math.sin(gamma[1]), -math.cos(gamma[2]), math.cos(beta)*math.sin(gamma[3])],
    [math.sin(beta)*math.sin(gamma[0]), math.sin(beta)*math.sin(gamma[1]), math.sin(beta)*math.sin(gamma[2]), math.sin(beta)*math.sin(gamma[3])]
    ], np.float)
    return np.array([
    [
        A[0,0] + A[0,1] + A[0,2] + A[0,3] ,
        A[1,0] + A[1,1] + A[1,2] + A[1,3] ,
        A[2,0] + A[2,1] + A[2,2] + A[2,3]
    ]
    ], np.float)*h_a

def dotconfMAT(gamma, beta):
    # diff. matrix A
    return np.array([
    [-math.cos(beta)*math.cos(gamma[0]), math.sin(gamma[1]), math.cos(beta)*math.cos(gamma[2]), -math.sin(gamma[3])],
    [-math.sin(gamma[0]), -math.cos(beta)*math.cos(gamma[1]), math.sin(gamma[2]), math.cos(beta)*math.cos(gamma[3])],
    [math.sin(beta)*math.cos(gamma[0]), math.sin(beta)*math.cos(gamma[1]), math.sin(beta)*math.cos(gamma[2]), math.sin(beta)*math.cos(gamma[3])]
    ], np.float)

def dither(tt, e_0):
    # 디더제어항
    return np.array([
    [ 1, e_0*math.sin(tt + math.pi), e_0*math.sin(tt + math.pi/2) ],
    [ e_0*math.sin(tt + math.pi), 1, e_0*math.sin(tt) ],
    [ e_0*math.sin(tt + math.pi/2), e_0*math.sin(tt), 1 ]
    ], np.float)

def w_rk45(w, u_real, k, dt, I_sc):
    dotw = lambda t, W: np.array([ np.linalg.inv(I_sc).dot( u_real - np.cross(w,I_sc.dot(w)) ) ])
    t_eval=np.linspace(dt*(k), dt*(k+1), 3)
    W = solve_ivp(dotw, [dt*(k), dt*(k+1)], w, method='RK45', t_eval=t_eval)
    return W.y[:,2]

def q_rk45(w, q_cur, k, dt):
    dotq = lambda t,Q: np.array([ 0.5*skew(w).dot(q_cur) ])
    t_eval=np.linspace(dt*(k), dt*(k+1), 3)
    Q = solve_ivp(dotq, [dt*(k), dt*(k+1)], q_cur, method='RK45', t_eval=t_eval)
    return Q.y[:,2]

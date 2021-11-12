# 파일 실행하려면 별도의 method들을 정리한 Utills.py파일이 필요함. 하단 링크로 들어가면 됨.
# https://github.com/yongja2005/CMG_HILS_Portfolio/blob/main/Utills.py

import numpy as np
import pandas as pd
import math
import Utills

# stepping motor conf.
step_angle = Utills.degtorad(0.05625) # min. step per 0.1sec

# Initial condition
dt = 0.1 #sec
time = np.arange(0, 50+dt, dt)
I_sc = np.array([
    [0.028, 0, 0], [0, 0.028, 0], [0, 0, 0.034]
    ], float) # kgm^2
beta = Utills.degtorad(53.13) #deg to rad
h_a = 1.312e-6*(100*math.pi) #Nmsec (100pi rad/s = 3000 RPM)
zeta = 0.707 #damping ratio
t_s = 15 #settling time(sec)
w_n = 4/zeta/t_s #natural frequency
k_D = 2*zeta*w_n*I_sc #D gain
k_P = 2*w_n*w_n*I_sc #P gain
limit_GV = math.pi #rad/s
e_0 = 0.01

# create empty array
w = np.zeros(shape=(len(time),3))
q_cur = np.zeros(shape=(len(time),4))
u_cmd = np.zeros(shape=(len(time),3))
gamma = np.zeros(shape=(len(time),4))
dotgamma = np.zeros(shape=(len(time),4))
dotgamma_real = np.zeros(shape=(len(time),4))
h_cmg = np.zeros(shape=(len(time),3))
u_cmg = np.zeros(shape=(len(time),3))
u_real = np.zeros(shape=(len(time),3))
RPY = np.zeros(shape=(len(time),3))
command_step = np.zeros(shape=(len(time),4))
m = np.zeros(shape=(len(time), 1))

#initial condition
w[0,:] = Utills.degtorad(0, 0, 0) #deg/s to rad/s (max = 0.1rad/s)
euler_0 = Utills.degtorad(0, 0, 0) #deg to rad
euler_cmd = Utills.degtorad(-36, 20, 30)
q_cur[0,:] = Utills.transf(euler_0[0,:])
q_cmd = Utills.transf(euler_cmd[0,:])
gamma[0,:] = Utills.degtorad(0, 0, 0, 0) # deg to rad
RPY[0,:] = Utills.transf(q_cur[0,:])

for k in range(len(time)-1):
    # Command torque
    q_e = Utills.skew(q_cmd[0,:]).dot(q_cur[k,:])
    u_cmd[k,:] = np.array([
        - ( k_P.dot(q_e[0:3].T) + k_D.dot(w[k,:].T) )
    ], float)

    # Angular Momentum about 3-axis of cmg : h_cmg
    h_cmg[k,:] = Utills.AngMomCMG(gamma[k,:],beta, h_a)

    # Torque about 3-axis of cmg : u_cmg
    u_cmg[k,:] = np.array([
        -( u_cmd[k,:].T + np.cross(w[k,:],h_cmg[k,:]) )
    ], float)

    # Gimbal velocity
    tt = dt*k
    E = Utills.dither(tt, e_0)
    dotA = Utills.dotconfMAT(gamma[k,:], beta)
    m[k,0] = np.linalg.det(dotA.dot(dotA.T))
    a = e_0*np.exp(-10*m[k,0])
    dotAA = dotA.T.dot(np.linalg.inv(dotA.dot(dotA.T) + a*E))
    dotgamma[k,:] = np.array([
        dotAA.dot(u_cmg[k,:].T)/h_a
    ])

    for i in range(0, 4): #dotgamma를 step(0.05625deg)별로 쪼개기
        if dotgamma[k,i] > limit_GV:
            temp = (limit_GV*dt)//step_angle
            dotgamma_real[k,i] = temp*step_angle/dt
        elif dotgamma[k,i] < -limit_GV:
            temp = -(limit_GV*dt)//step_angle
            dotgamma_real[k,i] = temp*step_angle/dt
        elif -step_angle < dotgamma[k,i]*dt < step_angle:
            dotgamma_real[k,i] = 0
        else:
            if dotgamma[k,i] > 0:
                temp = dotgamma[k,i]*dt//step_angle
                dotgamma_real[k,i] = temp*step_angle/dt
            else:
                temp = dotgamma[k,i]*dt//step_angle
                dotgamma_real[k,i] = (temp+1)*step_angle/dt

    gamma[k+1,:] = np.array([gamma[k,:]+ dotgamma_real[k,:]*dt], float)
    # step size of step motor
    command_step[k, :] = np.round(dotgamma_real[k, :]/step_angle*dt)

    # Real torque of cmg : u_real
    u_real[k,:] = np.array([
        -(h_a*dotA.dot(dotgamma_real[k,:]) + np.cross(w[k,:],h_cmg[k,:]))
    ], float)

    # Dynamic Eqauation
    w[k+1,:] = Utills.w_rk45(w[k,:], u_real[k,:], -1+k, dt, I_sc)

    # Kinematics Equation
    q_cur[k+1,:] = Utills.q_rk45(w[k+0,:], q_cur[k+0,:], -1+k, dt)

    #Quaternion -> euler angle(deg)
    RPY[k+1,:] = Utills.transf(q_cur[k+1,:])

print(gamma*180/math.pi)

# pandas DataFrame으로 변경 후 csv확장자로 저장
RPY_ = pd.DataFrame(RPY)
RPY_.to_csv('RPY_.csv', index=False, header=None)
w_ = pd.DataFrame(w)
w_.to_csv('w_.csv', index=False, header=None)
dotgamma_real_ = pd.DataFrame(dotgamma_real)
dotgamma_real_.to_csv('dotgamma_real_.csv', index=False, header=None)
gamma_ = pd.DataFrame(gamma)
gamma_.to_csv('gamma_.csv', index=False, header=None)
u_real_ = pd.DataFrame(u_real)
u_real_.to_csv('u_real_.csv', index=False, header=None)
m_ = pd.DataFrame(m)
m_.to_csv('m_.csv', index=False, header=None)


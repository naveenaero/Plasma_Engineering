import numpy as np
from math import sin,cos,exp,sqrt,radians
import matplotlib.pyplot as plt

g = 9.81

actual_x = lambda v,t,theta: 1 + v*t*cos(radians(theta))
actual_y = lambda v,t,theta: v*t*sin(radians(theta)) - 0.5*g*t**2 

dt = 0.01
time_horizon = 0.50
nt = int(time_horizon/dt)
drag_coefficient = 0

class Particle:
    def __init__(self, x, y, v, drag_coefficient, particle_type, direction):
        self.x = x
        self.y = y
        self.vx = v*cos(radians(direction))
        self.vy = v*sin(radians(direction))
        self.b = drag_coefficient
        if particle_type=='electron':
            self.charge_mass_ratio = -1.6e-19/9.1e-31
        elif particle_type=='ion':
            self.charge_mass_ratio = 1.6e-19*1835/9.1e-31
        else:
            self.charge_mass_ratio = 0

    def update(self, integrator):
        if integrator=='euler':
            self.euler()
            
        elif integrator=='rk2':
            self.rk2()

        elif integrator=='boris':
            self.boris_pusher()

        else:
            print 'Error!'

    def get_gravity_force(self, v, b):
        Fx = -b*v[0]
        Fy = -g - b*v[1]
        return np.array([Fx, Fy])          

    def euler(self):
        self.x += self.vx*dt
        self.y += self.vy*dt
        F = self.get_force([self.vx, self.vy], self.b)
        self.vx += F[0]*dt
        self.vy += F[1]*dt


    def rk2(self):
        kx = np.zeros(2)
        kv = np.zeros(2)

        xn = [self.x, self.y]
        vn = [self.vx, self.vy]

        kx = vn + self.get_force(vn, self.b)*dt/2
        kv = vn + self.get_force(kx, self.b)*dt

        self.x += kx[0]*dt
        self.y += kx[1]*dt
        self.vx = kv[0]
        self.vy = kv[1]

    def boris_pusher(self):
        v = np.zeros(2)
        x = np.zeros(2)
        vn = [self.vx, self.vy]
        xn = [self.x, self.y]
        

        v = self.v 
# class ChargedParticle:


def error(true, computed):
    error = np.linalg.norm(true-computed)/np.linalg.norm(true)
    print 'L2-error-norm = ', error

particle1 = Particle(1, 0, 4, drag_coefficient, 45)
particle2 = Particle(1, 0, 4, drag_coefficient, 45)

particle_position1 = np.zeros((nt+1,2))
particle_position2 = np.zeros((nt+1,2))

particle_position1[0,0] = particle1.x
particle_position1[0,1] = particle1.y
particle_position2[0,0] = particle2.x
particle_position2[0,1] = particle2.y

# analytical position parameters
true_position = np.zeros((nt+1,2))
time = np.arange(0,time_horizon + dt , dt)
true_position[:,0] = actual_x(4,time,45)
true_position[:,1] = actual_y(4,time,45)
r_true = (true_position[:,0]**2+true_position[:,1]**2)**0.5

for t in range(nt):
    particle1.update('euler')
    particle2.update('rk2')
    particle_position1[t+1,0] = particle1.x
    particle_position1[t+1,1] = particle1.y
    particle_position2[t+1,0] = particle2.x
    particle_position2[t+1,1] = particle2.y
    
    

rk2_pos = (particle_position2[:,0]**2+particle_position2[:,1]**2)**0.5
euler_pos = (particle_position1[:,0]**2+particle_position1[:,1]**2)**0.5

# error
error(r_true[-1], rk2_pos[-1])
error(r_true[-1], euler_pos[-1])
error(true_position[-1,1], particle_position2[-1,1])
error(true_position[-1,1], particle_position1[-1,1])


plt.plot(particle_position1[:,0],particle_position1[:,1],linestyle='-')
plt.plot(particle_position1[-1,0],particle_position1[-1,1],'ro')
plt.plot(particle_position2[:,0],particle_position2[:,1],linestyle='-')
plt.plot(particle_position2[-1,0],particle_position2[-1,1],'bo')
plt.plot(true_position[:,0], true_position[:,1])

plt.show()

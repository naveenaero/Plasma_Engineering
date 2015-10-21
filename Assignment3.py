import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import sin,cos,exp,sqrt,radians
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10

g = 9.81

actual_x = lambda v,t,theta: 1 + v*t*cos(radians(theta))
actual_y = lambda v,t,theta: v*t*sin(radians(theta)) - 0.5*g*t**2 

dt = 0.01
time_horizon = 0.01
nt = int(time_horizon/dt)
drag_coefficient = 0

class Particle:
    def __init__(self, x, y, v, drag_coefficient, direction):
        self.x = x
        self.y = y
        self.vx = v*cos(radians(direction))
        self.vy = v*sin(radians(direction))
        self.b = drag_coefficient
        
    def update(self, integrator):
        if integrator=='euler':
            self.euler()
            
        elif integrator=='rk2':
            self.rk2()
        
        else:
            print 'Error!'

    def get_force(self, v, b):
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

    
class ChargedParticle:
    def __init__(self, position, velocity, particle_type):
        self.E = np.array([0,0,0])
        self.B = np.array([0,0,0.01])
        
        self.velocity = velocity
        
        if particle_type == 'electron':
            self.charge_mass_ratio = -1*1.6e-19/9.1e-31
            larmor_radius = 9.1e-31*velocity[1]/(1.6e-19*self.B[2])
        elif particle_type == 'ion':
            self.charge_mass_ratio = 1835*1.6e-19/9.1e-31
            larmor_radius = 1835*9.1e-31*velocity[1]/(1.6e-19*self.B[2])
        else:
            print "Error, enter correct particle type"

        self.position = np.zeros(3)
        self.position[0] = larmor_radius
        self.position[1] = position[1]
        self.position[2] = position[2]
        

    def boris_pusher(self):
        vn = self.velocity
        xn = self.position

        v_minus = np.zeros(3)
        v_plus = np.zeros(3)
        t = np.zeros(3)
        s = np.zeros(3)

        t = self.charge_mass_ratio*self.B*dt/2
        s = 2*t/(1+np.linalg.norm(t)**2)

        v_minus = vn + self.charge_mass_ratio*self.E*dt/2
        
        v_plus = v_minus + np.cross(v_minus, s) + np.cross(np.cross(v_minus,t),s)
        print v_plus
        
        self.velocity = v_plus + self.charge_mass_ratio*self.E*dt/2
        print self.velocity
        
        self.position = xn + self.velocity*dt





def error(true, computed):
    error = np.linalg.norm(true-computed)/np.linalg.norm(true)
    print 'L2-error-norm = ', error

particle1 = Particle(1, 0, 4, drag_coefficient, 45)
particle2 = Particle(1, 0, 4, drag_coefficient, 45)

init_position = np.array([0,0,0])
init_velocity = np.array([0,1e5,0])
particle3 = ChargedParticle(init_position, init_velocity, 'electron')

particle_position1 = np.zeros((nt+1,2))
particle_position2 = np.zeros((nt+1,2))
particle_position3 = np.zeros((nt+1,3))

particle_position1[0,0] = particle1.x
particle_position1[0,1] = particle1.y
particle_position2[0,0] = particle2.x
particle_position2[0,1] = particle2.y
particle_position3[0,0] = particle3.position[0]
particle_position3[0,1] = particle3.position[1]
particle_position3[0,2] = particle3.position[2]

# analytical position parameters
true_position = np.zeros((nt+1,2))
time = np.arange(0,time_horizon + dt , dt)
true_position[:,0] = actual_x(4,time,45)
true_position[:,1] = actual_y(4,time,45)
r_true = (true_position[:,0]**2+true_position[:,1]**2)**0.5

for t in range(nt):
    particle1.update('euler')
    particle2.update('rk2')
    particle3.boris_pusher()
    particle_position1[t+1,0] = particle1.x
    particle_position1[t+1,1] = particle1.y
    particle_position2[t+1,0] = particle2.x
    particle_position2[t+1,1] = particle2.y
    particle_position3[t+1,0] = particle3.position[0]
    particle_position3[t+1,1] = particle3.position[1]
    particle_position3[t+1,2] = particle3.position[2]
    # print particle3.velocity, particle3.position

# rk2_pos = (particle_position2[:,0]**2+particle_position2[:,1]**2)**0.5
# euler_pos = (particle_position1[:,0]**2+particle_position1[:,1]**2)**0.5

# error
# error(r_true[-1], rk2_pos[-1])
# error(r_true[-1], euler_pos[-1])
# error(true_position[-1,1], particle_position2[-1,1])
# error(true_position[-1,1], particle_position1[-1,1])


# plt.plot(particle_position1[:,0],particle_position1[:,1],linestyle='-')
# plt.plot(particle_position1[-1,0],particle_position1[-1,1],'ro')

# plt.plot(particle_position2[:,0],particle_position2[:,1],linestyle='-')
# plt.plot(particle_position2[-1,0],particle_position2[-1,1],'bo')

plt.plot(particle_position3[:,0],particle_position3[:,1],linestyle='-')
plt.plot(particle_position3[-1,0], particle_position3[-1,0],'go')

# plt.plot(true_position[:,0], true_position[:,1])


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot(particle3.position[:,0], particle3.position[:,1], particle3.position[:,2], label='parametric curve')
# ax.legend()

# plt.show()

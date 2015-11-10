import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import Interpolation as interp1d

class Distribution:
	def __init__(self, name, m):
		# n - particle density
		self.n = 1;
		# m - particle mass
		self.m = m;
		# k - boltzmann constant
		self.k = 1.38e-23;
		# set name
		self.name = name
		
		# speeds, velocities, energies
		self.speed = 0
		self.velocity = 0
		self.energy = 0
		
		# distribution functions
		self.f_speed = 0
		self.f_velocity = 0
		self.f_ener = 0
		
		# temperature
		self.temperature = 0

	def set_temperature(self, temperature):
		self.temperature = temperature
	
	def get_scale_param(self):
		return np.sqrt(1/(self.k*self.temperature))

	def velocity_dist(self, velocity):
		a1 = np.exp(-self.m*np.dot(velocity,velocity)/(2*self.k*self.temperature))
		a2 = self.n*(self.m/(2*np.pi*self.k*self.temperature))**1.5
		f = a1*a2;
		return f

	def speed_dist(self, speed):
		a1 = np.exp(-self.m*speed**2/(2*self.k*self.temperature))
		a2 = self.n*(self.m/(2*np.pi*self.k*self.temperature))**1.5
		a3 = 4*np.pi*speed**2
		f = a1*a2*a3;
		return f

	def energy_dist(self, energy):
		a1 = self.n*(np.pi/(np.pi*self.k*self.temperature)**1.5)*np.sqrt(energy)
		a2 = np.exp(-energy/(self.k*self.temperature))
		f = a1*a2;
		return f

	def get_distribution(self, variable):
		if variable=='speed':
			self.speed = np.linspace(0,3000,3000)
			self.f_speed = np.zeros_like(self.speed)
			for i in range(len(self.speed)):
				self.f_speed[i] = self.speed_dist(self.speed[i])
		
		elif variable=='velocity':
			self.velocity = np.linspace(-2000,2000,4000)
			self.f_velocity = np.zeros_like(self.velocity)
			for i in range(len(self.velocity)):
				self.f_velocity[i] = self.velocity_dist(self.velocity[i])
		
		elif variable=='energy':
			self.energy = np.linspace(0,2000,3000)
			self.energy = self.m*self.energy**2/2
			self.f_energy = np.zeros_like(self.energy)
			for i in range(len(self.energy)):
				self.f_energy[i] = self.energy_dist(self.energy[i])

		else:
			print 'error'

	def get_mean(self, variable):
		if variable=='speed':
			return np.sum(np.multiply(self.f_speed, self.speed))
		elif variable=='velocity':
			return np.sum(np.multiply(self.f_velocity, self.velocity))
		elif variable=='energy':
			return np.sum(np.multiply(self.f_energy, self.energy))
		else:
			print "error"

	def get_rms(self, variable):
		if variable=='speed':
			return np.sqrt(np.sum(np.multiply(self.f_speed, self.speed**2)))
		elif variable=='velocity':
			return np.sqrt(np.sum(np.multiply(self.f_velocity, self.velocity**2)))
		elif variable=='energy':
			return np.sqrt(np.sum(np.multiply(self.f_energy, self.energy**2)))
		else:
			print 'error'
		
	def get_most_probable(self, variable):
		if variable=='speed':
			return self.speed[np.argmax(self.f_speed)]
		elif variable=='velocity':
			return self.velocity[np.argmax(self.f_velocity)]
		elif variable=='energy':
			return self.energy[np.argmax(self.f_energy)]
		else:
			print 'error'

	def get_characteristic_speeds(self, variable):
		mean = self.get_mean(variable)
		rms = self.get_rms(variable)
		most_probable = self.get_most_probable(variable)
		print "--------- T = ", self.temperature, "K ---------------"
		print 'mean for '+variable+' is: ', mean, 'm/s'
		print 'rms for '+variable+' is: ', rms, 'm/s'
		print 'most_probable for '+variable+' is: ', most_probable, 'm/s'



	def plot_distribution(self, variable):
		plt.rc('text', usetex=True)
		plt.rc('font', family='helvetica')
		scaling = self.get_scale_param()
		if variable=='speed':
			scaling *= np.sqrt(self.m)
			plt.plot(self.speed*scaling, self.f_speed, label='T='+str(self.temperature)+' K')
			plt.ylabel(r'$f(v)$', fontsize=16)
			plt.xlabel(r'$\displaystyle v\sqrt{\frac{m}{kT}}$',fontsize=16)
			plt.subplots_adjust(bottom=0.15)
			plt.title(variable)
		elif variable=='velocity':
			scaling *= np.sqrt(self.m)
			plt.plot(self.velocity*scaling, self.f_velocity, label='T='+str(self.temperature)+' K')
			plt.ylabel(r'$f(v_x)$', fontsize=16)
			plt.xlabel(r'$\displaystyle v_x\sqrt{\frac{m}{kT}}$',fontsize=16)
			plt.subplots_adjust(bottom=0.15)
			plt.title(variable)
		elif variable=='energy':
			scaling *= scaling
			scaling *= self.energy
			plt.plot(scaling, self.f_energy, label='T='+str(self.temperature)+' K')
			plt.ylabel(r'$f(\epsilon)$', fontsize=16)
			plt.xlabel(r'$\displaystyle \frac{\epsilon}{kT}$',fontsize=16)
			plt.subplots_adjust(bottom=0.15)
			plt.title(variable)
		else:
			print "plotting error"
	

	def compute_moments(self, a, b, variable, moment):
		[weights, roots] = interp1d.compute_lgl_weights(interp1d.legendre_interp, interp1d.lobatto_interp)
		if variable=='speed':
			temp = (b-a)/2
			integral = 0
			new_x = roots*(b-a)/2 + (a+b)/2
			for i in range(len(new_x)):
				integral += weights[i]*self.speed_dist(new_x[i])*(new_x[i]**moment)
			integral *= temp
			if moment==2:
				integral = integral*self.m/(3*self.n*self.k)
				self.computed_temperature = integral
				return integral
			else:
				self.computed_density = integral
				return integral
	
	def compute_error_norms(self):
		self.temperature_l2_err = abs(self.temperature - self.computed_temperature)/self.temperature
		self.density_l2_err = abs(self.n - self.computed_density)/self.n



def test_distribution_function(plotter):
	nitrogen = Distribution('nitrogen', m=4.64e-26)
	temperature = [300];
	variable = 'speed'
	for T in temperature:
		nitrogen.set_temperature(T)
		nitrogen.get_distribution(variable)
		nitrogen.get_characteristic_speeds(variable)
		if variable=='speed':
			print 0,'moment of speed = ', nitrogen.compute_moments(0, 3000, variable, 0)
			print 2,'moment of speed = ', nitrogen.compute_moments(0, 3000, variable, 2), 'K'
			nitrogen.compute_error_norms()
			print 'error in numerical calculation of temperature: ', nitrogen.temperature_l2_err
			print 'error in numerical calculation of density: ', nitrogen.density_l2_err
		if plotter:
			nitrogen.plot_distribution(variable)
	if plotter:
		plt.legend()
		plt.show()

test_distribution_function(False)


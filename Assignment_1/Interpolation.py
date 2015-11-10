import numpy as np
import scipy
import math
from scipy import interpolate
import scipy.special as sp
from scipy.special import legendre
import matplotlib.pylab as plt

runge = lambda x:1/(1+50*x**2)
cosine = lambda x:np.cos(np.pi*x/2)
gamma = lambda n:math.sqrt(2.0/(2*n+1))



def plot_grid(grid):
    plt.scatter(grid[:,0],grid[:,1])
    plt.grid()



def plot_error(x1, y1, method):
    y = function(x1)
    error = np.linalg.norm(y-y1)/np.linalg.norm(y)
    name = method.__class__.__name__
    plt.subplot(2,1,2)

    plt.plot(method.order, error, 'o', label=method.__class__.__name__)
    plt.legend()

def plot_interpolated_function(x1, y1, name):
    '''
    return plot object
    '''
    plt.subplot(2,1,1)
    plt.plot(x1, y1, label=name)
    plt.legend()

def compute_lg_weights(legendre_method):
    roots = legendre_method.compute_roots()
    weights = np.zeros(len(roots))
    poly = legendre_method.compute_poly()
    poly_derivative = np.polyder(poly,1)
    for i in range(len(weights)):
        numerator = 2
        denominator = np.polyval(poly_derivative,roots[i])**2
        denominator *= (1-roots[i]**2)
        weights[i] = 2/denominator
    return weights, roots

def compute_lgl_weights(legendre_method, lobatto_method):
    order = lobatto_method.order
    roots = lobatto_method.compute_roots()
    weights = np.zeros(len(roots))
    poly = legendre_method.compute_poly()
    for i in range(len(weights)):
        numerator = 2.0/(order*(order+1))
        denominator = np.polyval(poly,roots[i])**2
        weights[i] = numerator/denominator

    return weights, roots

def evaluate_interpolation_function(lagrange_polynomial, grid_param):
    nx = grid_param[2]
    xmin = grid_param[0]
    xmax = grid_param[1]
    x1 = np.linspace(xmin, xmax, nx)
    y1 = np.zeros(nx)
    for i in xrange(nx):
        for j in xrange(len(lagrange_polynomial)):
            y1[i] += np.polyval(lagrange_polynomial[j],x1[i])
    return x1, y1


def vandermonde1D(order, x):
    '''Returns a vandermonde matrix for 1D setting
       using given order and points x
    '''
    v1D = np.zeros((len(x),order+1))
    for i in range(order+1):
        v1D[:,i] = np.polyval(sp.jacobi(i,0,0),x)/gamma(i)
    return v1D



def interp(method):
    '''Interpolates the given function in the domain [-1,1]
        using a given method (which uses a set of specialised
        set of points)
    '''
    xprime = method.compute_roots()
    yprime = method.function(xprime)
    lagrange_polynomial = lagrange.naive(xprime, yprime)
    [x1, y1] = evaluate_interpolation_function(lagrange_polynomial, grid_param)
    plot_interpolated_function(x1, y1, method.__class__.__name__)
    plot_error(x1, y1, method)

def integrate(method, flag):
    '''Integrate the given function in the domain [-1,1]
       using the Legendre Gauss Lobatto/ Legendre Gauss
       quadrature points
    '''
    if flag == 1: #lgl weights
        xprime = method[1].compute_roots()
        yprime = method[1].function(xprime)
        [weights, roots] = compute_lgl_weights(method[0], method[1])
	print vandermonde1D(1, roots)
        
    else:         #lg weights
        xprime = method[0].compute_roots()
        yprime = method[0].function(xprime)
        [weights, roots] = compute_lg_weights(method[0])
    
    lagrange_polynomial = lagrange.naive(xprime, yprime)
    
    integration = 0
    
    for i in range(len(weights)):
        temp_integrate = 0
        for j in range(len(lagrange_polynomial)):
            temp_integrate += np.polyval(lagrange_polynomial[j],roots[i])
        temp_integrate *= weights[i]
        integration += temp_integrate
    return integration

class Interpolation:
    def __init__(self, x, y, function):
        self.x = x
        self.y = y
        self.n = len(x)
        self.order = self.n
        self.function = function

    def compute_roots(self):
        return self.x
    
class InterpLegendre(Interpolation):
    '''Interpolate using Lagrange polynomials
        using the roots of (N+1)th order Chebyshev polynomial
        where N is the number of points on the x axis
    '''
        
    def __init__(self, x, y, function):
        Interpolation.__init__(self, x, y, function)

    def compute_roots(self):
        poly = np.zeros(self.order+1)
        poly[-1] = 1
        return np.polynomial.legendre.Legendre(poly).roots()

    def compute_poly(self):
        return legendre(self.order)
    
class InterpLobatto(Interpolation):
    '''Interpolate using Lagrange polynomials
        using the roots of (N+1)th order Lobatto polynomial
        where N is the number of points on the x axis
    '''
    def __init__(self, x, y, function):
        Interpolation.__init__(self, x, y, function)
    
    def compute_roots(self):
        poly = np.zeros(self.order+1)
        poly[-1] = 1
        Legendre_derivative = np.polynomial.legendre.Legendre(poly).deriv(1)
        Legendre_roots = np.polynomial.legendre.Legendre.roots(Legendre_derivative)

        lobatto_roots = list(Legendre_roots)
        lobatto_roots.append(1)
        lobatto_roots = [-1] + lobatto_roots
        lobatto_roots = np.asarray(lobatto_roots)
        return lobatto_roots

class InterpChebyshev(Interpolation):
    '''Interpolate using Lagrange polynomials
        using the roots of (N+1)th order Chebyshev polynomial
        where N is the number of points on the x axis
    '''
    def __init__(self, x, y, function):
        Interpolation.__init__(self, x, y, function)

    def compute_roots(self):
        poly = np.zeros(self.order)
        poly[-1] = 1
        chebyshev_roots = np.polynomial.chebyshev.Chebyshev(poly).roots()
        return chebyshev_roots   

class InterpEquiSpaced(Interpolation):
    '''Interpolate using Lagrange polynomials
        using equally spaced points
    '''
    def __init__(self, x, y, function):
        Interpolation.__init__(self, x, y, function)


class LagrangePolynomial:
    def __init__(self, method):
        self.method = method

    def naive(self, x, y):
        '''A very naive approach to compute the lagrange
        polynomial coefficients using a direct algorithm
        '''
        N = len(x) 
        lagrange_polynomial = []
        
        for i in range(N):
            poly_numerator = np.array([0, 1.0])
            poly_denominator = 1.0
            for j in range(N):
                if j!=i:
                    temp_numerator = np.array([1.0, -x[j]])
                    poly_numerator = np.polymul(poly_numerator, temp_numerator)
                    poly_denominator *= x[i]-x[j]

            lagrange_polynomial.append(poly_numerator*y[i]/poly_denominator)
        return lagrange_polynomial


lagrange = LagrangePolynomial('naive')
function = cosine
nx = 20
N = nx
grid_param = [-1, 1, 51]


x = np.linspace(-1.0, 1.0, nx)
y = function(x)

legendre_interp = InterpLegendre(x, y, function)
lobatto_interp = InterpLobatto(x, y, function)


# interp(lobatto_interp)
# interp(legendre_interp)
# interp_2d(2)

# compute_legendre_gauss_weights(legendre_interp)
# print integrate([legendre_interp,lobatto_interp],1)
# plt.show()













import numpy as np
import pylab as py
import scipy.special as sp

def airy0(x):
    f = 1.14
    a = np.sin(x*f)/((x*f)**2)-np.cos(x*f)/(x*f)
    return 2.25*((2*a/(x*f))**2)
def airy(x):
    return (2*sp.j1(x)/x)**2
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x = np.linspace(-10, 10, 2000)
py.plot(x, airy(x),'r--',
        x, airy0(x),'b--',
        x, gaussian(x,0,3.83/3),'g--')

py.xlim((-10, 10))
py.ylim((-0.5, 1.1))
py.legend(('$(2\mathcal{J}_1(x)/x)^2$',
           r'$\mathrm{2.25}(\frac{sin(x)/x^2-cos(x)/x}{x})^2$',
           '$gaussian$'))
py.xlabel('$x$')
py.ylabel('Intensity')
                              
py.grid(True)
                                     
py.show()



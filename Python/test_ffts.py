import sys
if '/Users/NicoG/Dropbox/python_useful' not in sys.path:
    sys.path.append('/Users/NicoG/Dropbox/python_useful')
import numpy as np
import matplotlib.pyplot as plt
import ngobfft as tf

plt.close("all")

# ============================================================================|

nx, ny = 64, 32
Lx, Ly = 10., 2.

dx, dy = Lx/nx, Ly/ny
x = np.linspace(0, Lx-dx, nx)
y = np.linspace(0, Ly-dy, ny)

Y, X = np.meshgrid(y, x)

Fct = (np.cos(8*np.pi*X/Lx)*np.sin(4*np.pi*Y/Ly)
       + 0.5*np.cos(16*np.pi*Y/Ly + 0.3465) + 0.25)

kk = tf.k_of_x(x)
ll = tf.k_of_x(y)
L, K = np.meshgrid(ll, kk)

# ============================================================================|

FFF = tf.obfft2(x, y, Fct)

plt.figure(3)

plt.subplot(121)
plt.pcolor(X, Y, Fct)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis([x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
# plt.grid()

plt.subplot(122)
plt.pcolor(K, L, abs(FFF))
plt.xlabel('$k$')
plt.ylabel('$l$')
plt.axis([k.min(), k.max(), l.min(), l.max()])
plt.colorbar()

plt.tight_layout()
plt.show()

# ============================================================================|

Fct2 = np.real(tf.obifft2(kk, ll, FFF))

plt.figure(4)

plt.subplot(122)
plt.pcolor(X, Y, Fct2)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.axis([x.min(), x.max(), y.min(), y.max()])
plt.colorbar()
# plt.grid()

plt.subplot(121)
plt.pcolor(K, L, abs(FFF))
plt.xlabel('$k$')
plt.ylabel('$l$')
plt.axis([kk.min(), kk.max(), ll.min(), ll.max()])
plt.colorbar()

plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import math

# ---- Real Y_lm ----
def Y_lm(l, m, theta, phi):
    P = sp.lpmv(abs(m), l, np.cos(theta))
    N = np.sqrt((2*l+1)/(4*np.pi) * math.factorial(l-abs(m))/math.factorial(l+abs(m)))
    if m > 0:
        return N * P * np.cos(m*phi)
    elif m < 0:
        return N * P * np.sin(abs(m)*phi)
    else:
        return N * P

# ---- Hydrogen energy ----
def energy(n):
    return -13.6/(n*n)

# ---- PLOT ORBITAL ----
def plot_orbital(n, l, m):
    theta = np.linspace(0, np.pi, 200)
    phi   = np.linspace(0, 2*np.pi, 200)
    phi, theta = np.meshgrid(phi, theta)

    Y = Y_lm(l, m, theta, phi)
    r = np.abs(Y)

    # convert to cartesian
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(x, y, z, cmap="viridis", linewidth=0)

    ax.set_title(f"n={n}, l={l}, m={m}\nEnergy = {energy(n):.2f} eV")
    ax.set_axis_off()
    try: ax.set_box_aspect([1,1,1])
    except: pass

    plt.show()

# ---- USER INPUT ----
n = int(input("Enter n: "))
l = int(input("Enter l: "))
m = int(input("Enter m: "))

plot_orbital(n, l, m)


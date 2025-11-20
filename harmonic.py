import numpy as np
import matplotlib.pyplot as plt
from math import factorial
from numpy.polynomial.hermite import hermval

L = 100          
N = 5000          
dx = L / (N + 1)
x = np.linspace(-L/2, L/2, N + 2)  
x_int = x[1:-1]                    

# CONSTRUCT HAMILTONIAN MATRIX
H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2.0 / dx**2
    if i > 0:
        H[i, i - 1] = -1.0 / dx**2
    if i < N - 1:
        H[i, i + 1] = -1.0 / dx**2

V = 0.5 * x_int**2
H += np.diag(V)
eigvals, eigvecs = np.linalg.eigh(H)

# ANALYTICAL WAVEFUNCTION (Harmonic Oscillator)
def psi_analytical(n, x):
    coeff = [0]*n + [1]  # Hermite polynomial coefficients
    norm = 1.0 / np.sqrt((2**n) * factorial(n) * np.sqrt(np.pi))
    return norm * hermval(x, coeff) * np.exp(-x**2 / 2)

n_levels = [0, 1,2,10]
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

print("-------------------------------------------------")
print("n    E_numerical     E_analytical     ")
print("-------------------------------------------------")

for i, n in enumerate(n_levels):
    ax = axes.flatten()[i]

    # Numerical wavefunction
    psi_num = eigvecs[:, n].copy()
    psi_num /= np.sqrt(np.trapezoid(psi_num**2, x_int))
    psi_num = np.concatenate(([0.0], psi_num, [0.0]))

    # Analytical wavefunction
    psi_ana = psi_analytical(n, x)
   

    # Align sign 
    if np.trapezoid(psi_num * psi_ana, x) < 0:
        psi_num *= -1

         
    print(f"{n:<4d} {eigvals[n]:<15.6f} {n + 0.5:<15.6f} ")

    # Plot
    ax.plot(x, psi_ana, 'r--', label='Analytical')
    ax.plot(x, psi_num, 'b-', label='Numerical')
    ax.set_title(f"n={n}   E_num={eigvals[n]:.3f}   E_ana={n+0.5:.3f}")
    ax.set_xlabel("x")
    ax.set_ylabel(r'$\psi(x)$')
    ax.legend()
    ax.grid(True)

print("-------------------------------------------------")
plt.tight_layout()
plt.show()


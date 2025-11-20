import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib.pyplot as plt

def morse_potential(r):
    return D * (np.exp(-2 * a * r) - 2 * np.exp(-a * r))

a = 1.94         
D = 4.75         
re = 0.74        
hcut = 1973.27   
m = (940 / 2) * 1e6  


N = 2000
r = np.linspace(0, 2.5, N) - re   
h = r[1] - r[0]
V = morse_potential(r)

main_diag = (hcut**2 / (m * h**2)) + V
off_diag = np.full(N - 1, -hcut**2 / (2 * m * h**2))

E, psi = eigh_tridiagonal(main_diag, off_diag)

n_states = [0, 1, 2, 8]
for n in n_states:
    print(f"State n={n}: Energy = {E[n]:.4f} eV")
    y = psi[:, n]
    y /= np.sqrt(np.trapezoid(y**2, r))  # normalize
    plt.figure(figsize=(6, 5))
    plt.plot(r + re, y**2 + E[n], label=f'n={n}')
    plt.plot(r + re, V, 'k-', label='V(r)')
    plt.axvline(re, color='black', lw=1.0, ls='--', label='rₑ')
    plt.xlabel("r (Å)")
    plt.ylabel("Energy (eV)")
    plt.title("Morse Potential - H₂ Molecule")
    plt.legend()
    plt.grid(True)
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e, epsilon_0, pi

# -------------------------
# Grid
# -------------------------
Z = 1
r_max = 5e-10
N = 1200

r = np.linspace(1e-12, r_max, N)
dr = r[1] - r[0]

# -------------------------
# Coulomb Potential
# -------------------------
V = -Z * e**2 / (4 * pi * epsilon_0 * r)

# -------------------------
# Hamiltonian (explicit matrix form)
# -------------------------
H = np.zeros((N, N))

coef = -hbar**2 / (2 * m_e * dr**2)

for i in range(N):

    # Diagonal term
    H[i, i] = coef * (-2) + V[i]

    # Upper diagonal
    if i < N - 1:
        H[i, i+1] = coef * 1

    # Lower diagonal
    if i > 0:
        H[i, i-1] = coef * 1

# -------------------------
# Solve eigenvalue problem
# -------------------------
E, psi = np.linalg.eigh(H)

for idx, n in enumerate([0,1]):   # n=1,2

    u = psi[:, n]
    R = u / r

    P = r**2 * R**2
    P /= np.trapezoid(P, r)

    fig, ax1 = plt.subplots(figsize=(7,5))

    ax1.plot(r*1e10, P, label=f"n={n+1}")
    ax1.set_xlabel("r (Å)")
    ax1.set_ylabel("Radial Probability")
    ax1.legend()

    ax2 = ax1.twinx()
    ax2.plot(r*1e10, V/e, 'r--')
    ax2.set_ylabel("Potential (eV)", color='r')

    plt.title(f"Hydrogen Atom: n = {n+1}")
    plt.tight_layout()
    plt.show()

    print(f"Energy for n={n+1}: {E[n]/e:.4f} eV")


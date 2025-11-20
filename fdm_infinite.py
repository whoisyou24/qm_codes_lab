import numpy as np
import matplotlib.pyplot as plt

# --- Grid ---
L = 1.0
N = 200
dx = L / (N + 1)
x = np.linspace(0, L, N + 2)
x_int = x[1:-1]

# --- Hamiltonian (your loop version) ---
H = np.zeros((N, N))
for i in range(N):
    H[i, i] = 2.0 / dx**2
    if i > 0:
        H[i, i - 1] = -1.0 / dx**2
    if i < N - 1:
        H[i, i + 1] = -1.0 / dx**2

# --- Solve ---
eigvals, eigvecs = np.linalg.eigh(H)

def psi_analytical(n, x):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)

def extend(v):
    return np.concatenate(([0], v, [0]))

plt.figure(figsize=(10,4))

for idx, n in enumerate([1, 2], 1):
    v = eigvecs[:, n-1]

    # Sign fix
    if np.trapezoid(v * psi_analytical(n, x_int), x_int) < 0:
        v = -v

    # Normalize
    v /= np.sqrt(np.trapezoid(v*v, x_int))

    # --- Energy ---
    E = eigvals[n-1]

    # --- Plot ---
    plt.subplot(1, 2, idx)
    plt.plot(x, psi_analytical(n, x), label=f"Analytical")
    plt.plot(x, extend(v), '--', label="Numerical")

    plt.title(f"n={n}\nEnergy = {E:.4f}")
    plt.grid()
    plt.legend()

plt.tight_layout()
plt.show()


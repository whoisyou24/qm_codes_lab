import numpy as np
import matplotlib.pyplot as plt

# constants
hbarc = 197.3
m = 940.0
k = 100.0
b_values = [0.0, 10.0, 30.0]
pref = hbarc**2 / (2*m)

# grid
L = 6.0
N_total = 4000
r_full = np.linspace(-15, L, N_total)
dr = r_full[1] - r_full[0]
r = r_full[1:-1]
N = len(r)

# FULL MATRIX HAMILTONIAN (FDM 2nd derivative + V)
def build_H(V):
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = 2*pref/dr**2 + V[i]
        if i > 0:
            H[i, i-1] = -pref/dr**2
        if i < N-1:
            H[i, i+1] = -pref/dr**2
    return H

for b in b_values:
    V = 0.5*k*r**2 + (b/3)*r**3
    H = build_H(V)

    # solve
    E, psi = np.linalg.eigh(H)
    u0, u1 = psi[:,0], psi[:,1]

    # normalize + sign fix
    u0 /= np.sqrt(np.sum(u0**2)*dr)
    u1 /= np.sqrt(np.sum(u1**2)*dr)
    if u0[0] < 0: u0 = -u0
    if u1[0] < 0: u1 = -u1

    print(f"b={b:.1f}  E0={E[0]:.3f} MeV  E1={E[1]:.3f} MeV")

    # plot
    fig, axs = plt.subplots(1,2,figsize=(10,4))
    fig.suptitle(f"b={b:.1f}  |  E0={E[0]:.3f}  E1={E[1]:.3f}")

    # ground state
    axs[0].plot(r, u0, label='u0')
    axs[0].plot(r, V/np.max(V)*np.max(u0), 'k--', label='V (scaled)')
    axs[0].set_title("n=0")
    axs[0].grid()
    axs[0].legend()

    # excited state
    axs[1].plot(r, u1, label='u1')
    axs[1].plot(r, V/np.max(V)*np.max(u1), 'k--', label='V (scaled)')
    axs[1].set_title("n=1")
    axs[1].grid()
    axs[1].legend()

    plt.tight_layout()
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def vel_to_omega(u, v, dx, dy):  
    omega = (v[1:-1, 2:] - v[1:-1, :-2]) / 2 / dx - (u[2:, 1:-1] - u[:-2, 1:-1]) / 2 / dy  
    return omega
def compute_stream(omega_int, dx, dy, tol=1e-5, psi_init=None):
    """
    Jacobi iteration
    """
    if psi_init is None:
        psi_init = np.zeros((Ny, Nx))
    change = 1.
    iter = 0
    while change > tol:
        psi_new = psi_init.copy()
        psi_new[1:-1, 1:-1] = (omega_int \
                               + (psi_init[2:, 1:-1] + psi_init[:-2, 1:-1]) / dy / dy \
                               + (psi_init[1:-1, 2:] + psi_init[1:-1, :-2]) / dx / dx) / (2 / dx / dx + 2 / dy / dy)
        change = np.max(np.abs(psi_new - psi_init))
        psi_init = psi_new
        iter += 1
    end_signal = True if iter == 1 else False
    return psi_init, end_signal
def apply_bnd_omega_(psi, u, v, omega):
    omega[:, 0] = 2 / dx / dx * (psi[:, 0] - psi[:, 1]) - 2 / dx * v[:, 0]
    omega[:, -1] = 2 / dx / dx * (psi[:, -1] - psi[:, -2]) + 2 / dx * v[:, -1]
    omega[0, :] = 2 / dy / dy * (psi[0, :] - psi[1, :]) + 2 / dy * u[0, :]
    omega[-1, :] = 2 / dy / dy * (psi[-1, :] - psi[-2, :]) - 2 / dy * u[-1, :]
def forward_omega_FTCS(psi, omega):
    uw_x = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / 4 / dx / dy * (omega[1:-1, 2:] - omega[1:-1, :-2])
    vw_y = +(psi[1:-1, 2:] - psi[1:-1, :-2]) / 4 / dx / dy * (omega[2:, 1:-1] - omega[:-2, 1:-1])
    diff_y = nu * (omega[2:, 1:-1] + omega[:-2, 1:-1] - 2 * omega[1:-1, 1:-1]) / dy / dy
    diff_x = nu * (omega[1:-1, 2:] + omega[1:-1, :-2] - 2 * omega[1:-1, 1:-1]) / dx / dx
    return omega[1:-1, 1:-1] + dt * (uw_x + vw_y + diff_y + diff_x)

Nx = 101
Ny = 101
Lx = Ly = 1
x = np.linspace(0, Lx, Nx)
dx = x[-1] - x[-2]
y = np.linspace(0, Ly, Ny)
dy = y[-1] - y[-2]
nu = 0.001
dt = 0.002

xx, yy = np.meshgrid(x, y)

u = np.zeros_like(xx)
v = np.zeros_like(yy)
u[-1, 1:-1] = 1.

omega = np.zeros_like(xx)
omega[1:-1, 1:-1] = vel_to_omega(u, v, dx, dy)
psi, _ = compute_stream(omega[1:-1, 1:-1], dx, dy)

end_signal = False
out_iter = 0
while not end_signal:
    apply_bnd_omega_(psi, u, v, omega)
    omega[1:-1, 1:-1] = forward_omega_FTCS(psi, omega)
    psi, end_signal = compute_stream(omega[1:-1, 1:-1], dx, dy, psi_init=psi, tol=1e-5)

plt.figure(figsize=(10, 8))
plt.contourf(x, y, psi, levels=20, cmap='jet')
plt.colorbar(label='Vorticity')
plt.title('Vorticity Distribution')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
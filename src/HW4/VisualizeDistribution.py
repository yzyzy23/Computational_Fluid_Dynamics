import numpy as np
import matplotlib.pyplot as plt

def plot_temperature_contour(filename, grid_size=(15, 12)):
    origin_data = np.loadtxt(filename)
    data = origin_data.T
    x = np.linspace(0, grid_size[0], data.shape[1])
    y = np.linspace(0, grid_size[1], data.shape[0])
    X, Y = np.meshgrid(x, y)

    plt.figure(figsize=(10, 8))
    levels = np.arange(293.15, 373.15, 5)
    cs = plt.contourf(X, Y, data, levels=levels, cmap='jet', extend='both')
    contour = plt.contour(X, Y, data, levels=levels, colors='k', linewidths=0.5)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%dK')
    plt.xlabel('x (cm)', fontsize=12)
    plt.ylabel('y (cm)', fontsize=12)
    cbar = plt.colorbar(cs)
    cbar.set_label('Temperature (K)', rotation=270, labelpad=15)
    plt.title('Steady-State Temperature Distribution', pad=15)
    plt.grid(linestyle=':', color='gray', alpha=0.5)
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW4/TemperatureDistribution.png')
    plt.show()

if __name__ == "__main__":
    plot_temperature_contour("./Computational_Fluid_Dynamics/src/HW4/TemperatureDistribution.txt", grid_size=(15, 12))
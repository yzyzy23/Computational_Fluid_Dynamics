import numpy as np
import matplotlib.pyplot as plt

def plot_convergence_rate(filename):
    omega = []
    convergence_iter = []
    plt.figure(figsize=(10, 6))
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=0, vmax=1)

    with open(filename, 'r') as file:
        lines = file.readlines()
        h = 0.0
        color = cmap(norm(h))
        for line in lines:
            data = line.split()
            if len(data) == 1:
                if omega:
                    plt.plot(omega, convergence_iter, marker='o', linestyle='-', color=color, label=f'h={h:.2f}', linewidth=0.5, markersize=1)

                h = float(data[0])
                color = cmap(norm(h))
                omega = []
                convergence_iter = []
            else:
                omega.append(float(data[0]))
                convergence_iter.append(int(data[1]))
    plt.plot(omega, convergence_iter, marker='o', linestyle='-', color=color, label=f'h={h:.2f}', linewidth=0.5,markersize=1)
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Grid Space h (cm)', fontsize=12)
    plt.xlabel('Relaxation Factor ω', fontsize=12)
    plt.ylabel('Number of Iterations', fontsize=12)
    plt.title('Convergence Rate of Successive Over-Relaxation Method', fontsize=14)
    plt.legend(loc='best', title='h values', fontsize=10, framealpha=0.8)
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW4/ConvergenceRate.png')
    plt.show()

if __name__ == "__main__":
    plot_convergence_rate("./Computational_Fluid_Dynamics/src/HW4/ConvergenceRate.txt")
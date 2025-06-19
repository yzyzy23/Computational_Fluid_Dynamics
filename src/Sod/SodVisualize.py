import numpy
import matplotlib.pyplot as plt

def read_data(data_name):
    p_array = []
    rho_array = []
    u_array = []
    x_array = []
    with open(data_name, "r") as f:
        lines = f.readlines()
        for line in lines:
            x = float(line.split()[0])
            rho = float(line.split()[1])
            u = float(line.split()[2])
            p = float(line.split()[3])
            p_array.append(p)
            rho_array.append(rho)
            u_array.append(u)
            x_array.append(x)
    return numpy.array(x_array), numpy.array(rho_array), numpy.array(u_array), numpy.array(p_array)

def plot_data(x, rho, u, p):
    plt.figure(figsize=(12, 8))

    plt.subplot(3, 1, 1)
    plt.plot(x, rho, label='Density (rho)', color='blue')
    plt.xlabel('Position (x)')
    plt.ylabel('Density (rho)')
    plt.title('Density Distribution')
    plt.grid()
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(x, u, label='Velocity (u)', color='green')
    plt.xlabel('Position (x)')
    plt.ylabel('Velocity (u)')
    plt.title('Velocity Distribution')
    plt.grid()
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(x, p, label='Pressure (p)', color='red')
    plt.xlabel('Position (x)')
    plt.ylabel('Pressure (p)')
    plt.title('Pressure Distribution')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.savefig('./Computational_Fluid_Dynamics/photo/Sod/SodSolve.png')
    plt.show()

if __name__ == "__main__":
    x, rho, u, p = read_data("./Computational_Fluid_Dynamics/src/Sod/SodResult.txt")
    plot_data(x, rho, u, p)
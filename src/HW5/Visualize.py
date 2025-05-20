import numpy
import matplotlib.pyplot as plt

def plot_stream_function(data_name, h):
    data = numpy.zeros((100, 100))
    x = numpy.linspace(0, 1, int(1/h))
    y = numpy.linspace(0, 1, int(1/h))
    with open(data_name, "r") as f:
        lines = f.readlines()
        for line in lines:
            i = int(line.split()[0])
            j = int(line.split()[1])
            data[i][j] = float(line.split()[2])

    plt.figure(figsize=(10, 8))
    levels = numpy.arange(-1, 0, 0.01)
    contour = plt.contour(x, y, data, levels=levels, cmap="jet", linewidths=0.5)
    plt.clabel(contour, inline=True, fontsize=8, fmt='%d')
    plt.xlabel('x (cm)', fontsize=12)
    plt.ylabel('y (cm)', fontsize=12)
    plt.title('Stream Function Distribution', pad=15)
    plt.show()

if __name__ == "__main__":
    h = 0.01
    data_name = "./Computational_Fluid_Dynamics/src/HW5/LidDrivenCavity.txt"
    plot_stream_function(data_name, h)
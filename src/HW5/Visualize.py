import numpy
import matplotlib.pyplot as plt

def read_data(data_name, h):
    data = numpy.zeros((int(1/h), int(1/h)))
    x = numpy.linspace(0, 1, int(1/h))
    y = numpy.linspace(0, 1, int(1/h))
    with open(data_name, "r") as f:
        lines = f.readlines()
        for line in lines:
            i = int(line.split()[0])
            j = int(line.split()[1])
            data[j][i] = float(line.split()[2])
    u = numpy.gradient(data, h, axis=0)
    v = -numpy.gradient(data, h, axis=1)
    return data, x, y, u, v

def plot_contour(x, y, data, name=None):
    plt.figure(figsize=(10, 8))
    plt.contourf(x, y, data, levels=50, cmap="jet")
    plt.colorbar(label=name, orientation='vertical')
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title(name, pad=15)
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW5/' + name + '.png')
    plt.close()

def plot_stream_line(x, y, u, v):
    speed = numpy.sqrt(u**2 + v**2)

    plt.figure(figsize=(10, 8))
    strm = plt.streamplot(x, y, u, v, color=speed, cmap="jet", density=1.5)
    plt.colorbar(strm.lines, label='Velocity', orientation='vertical')
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.title('Streamlines', pad=15)
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW5/StreamlinePlot.png')
    plt.close()

def plot_midline_velocity(x, y, u, v):
    plt.figure(figsize=(10, 8))
    midline_velocity = v[int(len(v)/2), :]
    plt.plot(x, midline_velocity, label='Midline Velocity', color='blue')
    plt.xlabel('x', fontsize=12)
    plt.ylabel('v', fontsize=12)
    plt.title('MidlineVelocity', pad=15)
    plt.legend()
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW5/MidlineVelocityV.png')
    plt.close()

    plt.figure(figsize=(10, 8))
    midline_velocity = u[:, int(len(u)/2)]
    plt.plot(y, midline_velocity, label='MidlineVelocity', color='blue')
    plt.xlabel('y', fontsize=12)
    plt.ylabel('u', fontsize=12)
    plt.title('MidlineVelocity', pad=15)
    plt.legend()
    plt.savefig('./Computational_Fluid_Dynamics/photo/HW5/MidlineVelocityU.png')
    plt.close()

def locate_vortices(data, h):
    for i in range(1, len(data) - 1):
        for j in range(1, len(data[0]) - 1):
            if (data[i][j] < data[i-1][j] and
                data[i][j] < data[i+1][j] and
                data[i][j] < data[i][j-1] and
                data[i][j] < data[i][j+1]):
                print(f"Vortex located at ({j * h}, {i * h}) with stream function value {data[i][j]}")

if __name__ == "__main__":
    h = 0.01
    data_name = "./Computational_Fluid_Dynamics/src/HW5/StreamFunction.txt"
    data, x, y, u, v = read_data(data_name, h)
    plot_contour(x, y, data, name="StreamFunction")
    plot_stream_line(x, y, u, v)
    plot_midline_velocity(x, y, u, v)
    locate_vortices(data, h)
    data_name = "./Computational_Fluid_Dynamics/src/HW5/Vorticity.txt"
    data, x, y, _, _ = read_data(data_name, h)
    plot_contour(x, y, data, name="Vorticity")
    print("All figures saved in ./Computational_Fluid_Dynamics/photo/HW5/")
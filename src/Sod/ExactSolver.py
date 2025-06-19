import numpy as np
import matplotlib.pyplot as plt


def pressure_function(p, state, gamma):
    rho, v, p_k = state
    if p > p_k:  # 激波情况
        A = 2 / ((gamma + 1) * rho)
        B = (gamma - 1) / (gamma + 1) * p_k
        f = (p - p_k) * np.sqrt(A / (p + B))
        df = np.sqrt(A / (p + B)) * (1 - (p - p_k) / (2 * (p + B)))
    else:  # 稀疏波情况
        a = np.sqrt(gamma * p_k / rho)
        f = 2 * a / (gamma - 1) * ((p / p_k) ** ((gamma - 1) / (2 * gamma)) - 1)
        df = 1 / (rho * a) * (p / p_k) ** (- (gamma + 1) / (2 * gamma))
    return f, df


def initial_pressure_guess(state_left, state_right, gamma):
    _, _, p_left = state_left
    _, _, p_right = state_right
    return 0.5 * (p_left + p_right)


def solve_contact_discontinuity(state_left, state_right, gamma):
    v_left, v_right = state_left[1], state_right[1]
    delta_v = v_right - v_left

    max_iterations = 100
    tolerance = 1e-12
    pressure_old = initial_pressure_guess(state_left, state_right, gamma)

    for _ in range(max_iterations):
        f_left, df_left = pressure_function(pressure_old, state_left, gamma)
        f_right, df_right = pressure_function(pressure_old, state_right, gamma)

        pressure_new = pressure_old - (f_left + f_right + delta_v) / (df_left + df_right)
        pressure_new = max(pressure_new, tolerance)  # 防止负压力

        if 2 * abs(pressure_new - pressure_old) / (pressure_new + pressure_old) < tolerance:
            velocity = 0.5 * (v_left + v_right + f_right - f_left)
            return pressure_new, velocity

        pressure_old = pressure_new

    raise RuntimeError("Not converged after maximum iterations.")


def sample_solution(p_m, v_m, state_left, state_right, gamma, s):
    rho_left, v_left, p_left = state_left
    a_left = np.sqrt(gamma * p_left / rho_left)

    rho_right, v_right, p_right = state_right
    a_right = np.sqrt(gamma * p_right / rho_right)

    if s < v_m:  # 左侧区域
        if p_m < p_left:  # 左稀疏波
            s_left = v_left - a_left
            a_m_left = a_left * (p_m / p_left) ** ((gamma - 1) / (2 * gamma))
            s_m_left = v_m - a_m_left
            if s < s_left:
                return rho_left, v_left, p_left
            elif s < s_m_left:
                rho = rho_left * (2 / (gamma + 1) + (gamma - 1) / ((gamma + 1) * a_left) * (v_left - s)) ** (2 / (gamma - 1))
                v = 2 / (gamma + 1) * (a_left + (gamma - 1) / 2 * v_left + s)
                p = p_left * (2 / (gamma + 1) + (gamma - 1) / ((gamma + 1) * a_left) * (v_left - s)) ** (2 * gamma / (gamma - 1))
                return rho, v, p
            else:
                return rho_left * (p_m / p_left) ** (1 / gamma), v_m, p_m
        else:  # 左激波
            s_shock = v_left - a_left * np.sqrt((gamma + 1) * p_m / (2 * gamma * p_left) + (gamma - 1) / (2 * gamma))
            if s < s_shock:
                return rho_left, v_left, p_left
            else:
                rho = rho_left * (p_m / p_left + (gamma - 1) / (gamma + 1)) / ((gamma - 1) * p_m / ((gamma + 1) * p_left) + 1)
                return rho, v_m, p_m
    else:  # 右侧区域
        if p_m < p_right:  # 右稀疏波
            s_right = v_right + a_right
            a_m_right = a_right * (p_m / p_right) ** ((gamma - 1) / (2 * gamma))
            s_m_right = v_m + a_m_right
            if s > s_right:
                return rho_right, v_right, p_right
            elif s > s_m_right:
                rho = rho_right * (2 / (gamma + 1) - (gamma - 1) / ((gamma + 1) * a_right) * (v_right - s)) ** (2 / (gamma - 1))
                v = 2 / (gamma + 1) * (-a_right + (gamma - 1) / 2 * v_right + s)
                p = p_right * (2 / (gamma + 1) - (gamma - 1) / ((gamma + 1) * a_right) * (v_right - s)) ** (2 * gamma / (gamma - 1))
                return rho, v, p
            else:
                return rho_right * (p_m / p_right) ** (1 / gamma), v_m, p_m
        else:  # 右激波
            s_shock = v_right + a_right * np.sqrt((gamma + 1) * p_m / (2 * gamma * p_right) + (gamma - 1) / (2 * gamma))
            if s > s_shock:
                return rho_right, v_right, p_right
            else:
                rho = rho_right * (p_m / p_right + (gamma - 1) / (gamma + 1)) / ((gamma - 1) * p_m / ((gamma + 1) * p_right) + 1)
                return rho, v_m, p_m


def compute_domain(state_left, state_right, gamma, t):
    v_left, p_left = state_left[1], state_left[2]
    v_right, p_right = state_right[1], state_right[2]
    a_left = np.sqrt(gamma * p_left / state_left[0])
    a_right = np.sqrt(gamma * p_right / state_right[0])
    return 2.0 * max(abs(v_left) + a_left, abs(v_right) + a_right) * t


def draw_solution(state_left, state_right, gamma, t=0.2, num_points=100):
    L = compute_domain(state_left, state_right, gamma, t)
    x = np.linspace(-L, L, num=num_points)
    solution = np.empty((num_points, 3), dtype=float)

    p_m, v_m = solve_contact_discontinuity(state_left, state_right, gamma)

    for i in range(num_points):
        solution[i, :] = sample_solution(p_m, v_m, state_left, state_right, gamma, x[i] / t)

    np.savetxt('solution.csv', solution, delimiter=',')
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1)
    plt.plot(x, solution[:, 0], label="Density")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(x, solution[:, 1], label="Velocity")
    plt.xlabel("x")
    plt.ylabel("Velocity")
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(x, solution[:, 2], label="Pressure")
    plt.xlabel("x")
    plt.ylabel("Pressure")
    plt.legend()

    plt.tight_layout()
    plt.savefig("./Computational_Fluid_Dynamics/photo/Sod/ExactSolve.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    gamma = 1.4
    input_left = [1.0, 0.0, 1.0]
    input_right = [0.125, 0.0, 0.1]
    draw_solution(input_left, input_right, gamma)
import matplotlib.pyplot as plt
import numpy as np

def plot_rod(q, ctime):
    plt.figure()  # This line ensures a new figure is created each time the function is called

    x = q[::2]
    y = q[1::2]
    plt.clf()
    plt.plot(x, y, 'ko-')
    plt.title(f"t = {ctime:.4f}s")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.pause(0.01)

def plot_results(totalTime, pos, vel, angle, saveImage=False):
    t = np.linspace(0, totalTime, len(pos))

    plt.figure()
    plt.plot(t, pos, 'ko-')
    plt.xlabel("Time (s)")
    plt.ylabel("Middle Node Position (m)")
    plt.title("Middle Node Position vs. Time")
    if saveImage: plt.savefig("data/results/position.png")

    plt.figure()
    plt.plot(t, vel, 'ko-')
    plt.xlabel("Time (s)")
    plt.ylabel("Middle Node Velocity (m/s)")
    plt.title("Middle Node Velocity vs. Time")
    if saveImage: plt.savefig("data/results/velocity.png")

    plt.figure()
    plt.plot(t, angle, 'ko-')
    plt.xlabel("Time (s)")
    plt.ylabel("Turning Angle (°)")
    plt.title("Turning Angle vs. Time")
    if saveImage: plt.savefig("data/results/angle.png")

    plt.show()




dt = 0.1 # Time step size
maxTime = 10   # total time of simulation
t = np.arange(0, maxTime + dt, dt)

# free indices
free_DOF = np.concatenate((np.arange(2, 4), np.arange(6, 8)))

# A list of the specific time points where you want to plot
plot_times = np.array([0, 0.1, 10, 100]) # This tells my plotting function which times stamps to plot

# Container to store y-coordinate of middle node 1
y_middle1 = np.zeros(len(t))
y_middle1[0] = x_old[3] # y-coordinate of middle node 1

# Container to store y-coordinate of middle node 3
y_middle3 = np.zeros(len(t))
y_middle3[0] = x_old[3] # y-coordinate of middle node 3

plot(x_old, index_matrix, t[0])
# print(x_old)
# print(index_matrix)
# print(t)

for k in range(len(t)-1):
  t_new = t[k+1]
  x_new, u_new = myInt(t_new, x_old, u_old, free_DOF, stiffness_matrix, index_matrix, m, dt)
  # if k % 10 == 0:
  if np.any(np.isclose(t_new, plot_times)):
    plot(x_new, index_matrix, t_new)

  y_middle1[k+1] = x_new[3]
  y_middle3[k+1] = x_new[7]

  x_old = x_new
  u_old = u_new

# Plot Node 1
plt.figure()
plt.plot(t, y_middle1, 'ro-')
plt.xlabel('Time (s)')
plt.ylabel('Middle Node 1 y-coordinate')
plt.title('Middle Node 1 y-coordinate vs. Time')
plt.grid(True)
plt.show()

# Plot Node 3
plt.figure()
plt.plot(t, y_middle3, 'ro-')
plt.xlabel('Time (s)')
plt.ylabel('Middle Node 3 y-coordinate')
plt.title('Middle Node 3 y-coordinate vs. Time')
plt.grid(True)
plt.show()
import numpy as np
import matplotlib.pyplot as plt



def generate_random_point_cloud():
    return np.random.rand(1000, 3)



d = generate_random_point_cloud()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(d[:,0], d[:,1], d[:,2])

plt.show()

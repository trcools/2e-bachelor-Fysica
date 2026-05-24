import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

def normalize(v):
    return v / np.linalg.norm(v)

# Tetrahedron vertices (symmetric choice)
verts = np.array([
    [1, 1, 1],
    [1, -1, -1],
    [-1, 1, -1],
    [-1, -1, 1],
], dtype=float)
verts = np.array([normalize(v) for v in verts])

# Twofold axes (midpoints of opposite edges) are along coordinate axes
axes2 = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
], dtype=float)
axes2 = np.array([normalize(v) for v in axes2])

# Prepare sphere
phi, theta = np.mgrid[0.0:np.pi:60j, 0.0:2.0*np.pi:120j]
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, rstride=4, cstride=4, color='lightblue', alpha=0.15, linewidth=0)

# Plot 3-fold vertices (4 points)
ax.scatter(verts[:,0], verts[:,1], verts[:,2], color='red', s=80, label='3-fold vertices (4)')  # type: ignore[arg-type]
for i,v in enumerate(verts):
    ax.text(v[0]*1.05, v[1]*1.05, v[2]*1.05, f'v{i+1}', color='red')

# Plot their negatives (face centers)
neg = -verts
ax.scatter(neg[:,0], neg[:,1], neg[:,2], color='orange', s=50, label='face centers (negatives)')  # type: ignore[arg-type]

# Plot 2-fold axes points (± for each axis)
pts2 = np.vstack([axes2, -axes2])
ax.scatter(pts2[:,0], pts2[:,1], pts2[:,2], color='green', s=60, label='2-fold axes points (6)')  # type: ignore[arg-type]
for i,v in enumerate(pts2):
    ax.text(v[0]*1.05, v[1]*1.05, v[2]*1.05, f'a{i+1}', color='green')

ax.set_box_aspect([1,1,1])
ax.view_init(elev=20, azim=30)
ax.set_axis_off()
ax.legend(loc='upper left')

plt.tight_layout()
output_path = Path(__file__).resolve().parents[1] / 'Figuren' / 'orbits_tetrahedron.png'
plt.savefig(output_path, dpi=200)
print(f'Saved {output_path}')

if __name__ == '__main__':
    pass

import vpython as vp
import numpy as np

n = 4

class atom(vp.sphere):
    def __init__(self, k, m, pos_, color_):
        self.k = k
        self.m = m
        self.sphere = vp.sphere(radius=0.2, pos=pos_, color=color_)

class lattice(atom):
    def __init__(self, k_1, k_2, m_1, m_2):
        self.k1 = k_1
        self.k2 = k_2
        self.m1 = m_1
        self.m2 = m_2
        self.grid = [
            [
                atom(self.k1, self.m1, vp.vector(0, 0, 1)*i + vp.vector(0, 1, 0)*j, vp.color.red)
                if (i+j) % 2 == 0 else
                atom(self.k2, self.m2, vp.vector(0, 0, 1)*i + vp.vector(0, 1, 0)*j, vp.color.blue)
                for i in range(0, n+2)
            ] for j in range(0, n+2)
        ]

    def GetPos(self):
        return np.array([[(self.grid[i][j].sphere.pos) for i in range(n+2)] for j in range(n+2)])

l = lattice(1, 2, 3, 2)
k = [10, 10]
r_eq = l.GetPos()
r_new = r_eq.copy()
a_new = np.array([vp.vector(0, 0, 0)] * (n+2)**2).reshape((n+2, n+2))
v_new = a_new.copy()
r_new[4, 1] = r_eq[4, 1] + vp.vector(0, 0, 0.3)
m = [1, 5]

def UpdateValues(a, r, v, dt):
    for i in range(1, n+1):
        for j in range(1, n+1):
            a[i, j] = (
                k[(i+j) % 2] / m[(i+j) % 2] * ((1, 1) <= (i+1, j) <= (4, 4)) *
                ((vp.mag(r[i+1, j] - r[i, j]) - 1) * (vp.norm(r[i+1, j] - r[i, j]))) +
                k[(i+j) % 2] / m[(i+j) % 2] * ((1, 1) <= (i-1, j) <= (4, 4)) *
                (vp.mag(r[i-1, j] - r[i, j]) - 1) * (vp.norm(r[i-1, j] - r[i, j])) +
                k[(i+j) % 2] / m[(i+j) % 2] * ((1, 1) <= (i, j+1) <= (4, 4)) *
                (vp.mag(r[i, j+1] - r[i, j]) - 1) * (vp.norm(r[i, j+1] - r[i, j])) +
                k[(i+j) % 2] / m[(i+j) % 2] * ((1, 1) <= (i, j-1) <= (4, 4)) *
                (vp.mag(r[i, j-1] - r[i, j]) - 1) * (vp.norm(r[i, j-1] - r[i, j]))
            )
    v_ = v + a * dt
    r_ = r + v * dt
    return v_, r_

while True:
    vp.rate(1000)
    v_new, r_new = UpdateValues(a_new, r_new, v_new, 0.001)
    for i in range(1, n+1):
        for j in range(1, n+1):
            l.grid[i][j].sphere.pos = r_new[i, j]

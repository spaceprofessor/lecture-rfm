"""
This program was created exclusively for the course "Raumflugmechanik"
(Spaceflight Dynamics) at FH Aachen University of Applied Sciences
by Prof. Dr. Bernd Dachwald on 03 November 2021.
"""

import numpy as np
import matplotlib.pyplot as plt


class Earth:
    mu = 398600.4  # Gravitational parameter in [km^3/s^2]
    R = 6378.14  # Radius in [km]
    J2 = 0.00108263  # J2


e = Earth()  # Create Earth object


class Orbit:
    def __init__(self):
        self.a = 0.0  # Semi-major axis in [km]
        self.ecc = 0.0  # Eccentricity
        self.rp = 0.0  # Pericenter radius in [km]
        self.ra = 0.0  # Apocenter radius in [km]
        self.n = 0.0  # Mean motion in [rad/s]

    def set_circ_orbit(self, r):
        self.a = r
        self.rp = r
        self.ra = r
        self.n = np.sqrt(e.mu / np.power(r, 3))

    def set_elliptic_orbit(self, hp, ha):
        self.rp = e.R + hp
        self.ra = e.R + ha
        self.a = (self.rp + self.ra) / 2
        self.ecc = 1 - self.rp / self.a
        self.n = np.sqrt(e.mu / np.power(self.a, 3))


def to_deg_per_day(angle):  # Converts from [rad/s] to [°/d]
    return angle * (180 / np.pi) * 86400


def calculate_prefactor(orb):  # Calculates prefactor in [°/d]
    pref = -3 * orb.n * e.J2 * np.power(e.R, 2) / (2 * np.power(orb.a, 2) * np.power(1 - np.power(orb.ecc, 2), 2))
    return to_deg_per_day(pref)


def rotation_nodal_line(i, orb):  # Calculates dOmega / dt in [°/d]
    return calculate_prefactor(orb) * np.cos(np.deg2rad(i))


def rotation_apsidal_line(i, orb):  # Calculates domega / dt in [°/d]
    return calculate_prefactor(orb) * (5/2 * np.power(np.sin(np.deg2rad(i)), 2) - 2)


def plot_rot_nodal_line(inc, orb):
    fig, ax = plt.subplots()
    for j in range(0, len(orb)):
        h = int(np.ceil(orb[j].a - e.R))
        ax.plot(inc, rotation_nodal_line(inc, orb[j]), '-', label='h = {} km'.format(h))
    ax.legend()


def plot_rot_apsidal_line(inc, orb):
    fig, ax = plt.subplots()
    for j in range(len(orb)):
        ha = int(np.ceil(orb[j].ra - e.R))
        ax.plot(inc, rotation_apsidal_line(inc, orb[j]), '-', label='ha = {} km'.format(ha))
    ax.legend()


h = [200.0, 500.0, 1000.0, 2500.0, 5000.0]  # List of altitudes in [km]
orbit_nodal_rotation = [Orbit() for j in range(0, len(h))]
for j in range(0, len(h)):
    orbit_nodal_rotation[j].set_circ_orbit(e.R + h[j])

hp = 200.0  # Pericenter altitude for apsidal rotation in [km]
ha = [200.0, 1000.0, 2500.0, 5000.0]  # List of apocenters for apsidal rotation in [km]
orbit_apsidal_rotation = [Orbit() for j in range(0, len(ha))]
for j in range(0, len(ha)):
    orbit_apsidal_rotation[j].set_elliptic_orbit(hp, ha[j])

i = np.arange(0.0, 180.0, 1.0)  # Range of inclinations in [°]

plot_rot_nodal_line(i, orbit_nodal_rotation)
plt.grid()
plt.xlabel("Inclination [°]")
plt.ylabel("Rotation of nodal line [°/d]")
plt.show()

plot_rot_apsidal_line(i, orbit_apsidal_rotation)
plt.grid()
plt.xlabel("Inclination [°]")
plt.ylabel("Rotation of apsidal line [°/d]")
plt.show()

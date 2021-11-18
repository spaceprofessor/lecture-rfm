"""
This program was created exclusively for the course "Raumflugmechanik"
(Spaceflight Dynamics) at FH Aachen University of Applied Sciences
by Prof. Dr. Bernd Dachwald on 18 November 2021.
"""

import numpy as np
import matplotlib.pyplot as plt


mu = 398600.4  # Gravitational parameter of Earth in [km^3/s^2]
 

def sin(x): return np.sin(x)  # sin and cos are part of the NumPy library. This just makes the reading of the code easier
def cos(x): return np.cos(x)  # because one can write sin() instead of np.sin() and cos() instead of np.cos()
def norm(x): return np.linalg.norm(x)  # ... and the norm / magnitude of a vector
pi = np.pi  # also not part of the Python core language but the NumPy library


# Input spacecraft state to be transformed
# === Begin user input ===========================================
rvec = np.array([10000.0, 40000.0, -5000.0])  # in [km]
vvec = np.array([-1.5, 1.0, -0.1])            # in [km/s]
# === End user input =============================================

r = norm(rvec)
v = norm(vvec)


# Step 1 - semi-major axis:
a = r / (2 - r * v**2 / mu)
print("Semi-major axis                : %f km" %a)

# Step 2 - eccentricity:
evec = (v**2 / mu - 1 / r) * rvec - 1 / mu * np.dot(rvec, vvec) * vvec
e = norm(evec)
print("Eccentricity                   : %f" %e)

# Step 3 - inclination:
hvec = np.cross(rvec, vvec)
h = norm(hvec)
i = np.arccos(hvec[2]/h)
print("Inclination                    : %f°" %np.rad2deg(i))

# Step 4 - longitude of the ascending node:
nvec = np.cross(([0.0 ,0.0, 1.0]), hvec)
n = norm(nvec)
Omega = np.arccos(nvec[0]/n)
if nvec[1] < 0.0: Omega = 2*pi - Omega
print("Longitude of the ascending node: %f°" %np.rad2deg(Omega))

# Step 5 - argument of pericenter:
omega = np.arccos(np.dot(nvec, evec)/(n * e))
if evec[2] < 0.0: omega = 2*pi - omega
print("Argument of pericenter         : %f°" %np.rad2deg(omega))

# Step 6 - true anomaly
f = np.arccos(np.dot(evec, rvec)/(e * r))
if np.dot(rvec, vvec) < 0.0: f = 2*pi - f
print("True anomaly                   : %f°" %np.rad2deg(f))

# Reverse conversion to check for correctness of first conversion

# Step 1 - radius, speed, and orbital angular momentum:
r = a * (1 - e**2) / (1 + e * np.cos(f))
v = (2*mu / r - mu / a)**0.5
h = (mu * a * (1 - e**2))**0.5


# Step 2 - position and velocity vector:
theta = omega + f
rvec = r * np.array([
    cos(Omega) * cos(theta) - sin(Omega) * sin(theta) * cos(i),
    sin(Omega) * cos(theta) + cos(Omega) * sin(theta) * cos(i),
                                           sin(theta) * sin(i)])
print("r =", rvec)
f1 = sin(theta) + e * sin(omega)
f2 = cos(theta) + e * cos(omega)
vvec = - mu / h * np.array([
    cos(Omega) * f1 + sin(Omega) * f2 * cos(i),
    sin(Omega) * f1 - cos(Omega) * f2 * cos(i),
                    -              f2 * sin(i)])
print("v =", vvec)
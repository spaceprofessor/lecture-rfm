"""
This program was created exclusively for the course "Raumflugmechanik"
(Spaceflight Dynamics) at FH Aachen University of Applied Sciences
by Prof. Dr. Bernd Dachwald on 13 October 2021.
"""

import numpy as np
import matplotlib.pyplot as plt


class Earth:
    mu = 398600.4  # Gravitational parameter in [km^3/s^2]
    R = 6378.14  # Radius in [km]


h = 400  # Altitude in [km]

e = Earth()

r = e.R + h
V = np.sqrt(e.mu/r)

print("The velocity in a circular Earth orbit with an altitude of", h, "km is", np.around(V, decimals=4), "km/s.")

# Plot velocity over altitude

h = np.linspace(0, 1200)
r = e.R + h
V = np.sqrt(e.mu/r)
plt.plot(h, V)
plt.grid()
plt.xlabel("Altitude [km]")
plt.ylabel("Velocity [km/s]")
plt.show()

# Plot orbital period over radius (up to 45000 km)

r = np.arange(0, 45000)
P = 2 * np.pi * np.sqrt(r**3/e.mu) / 3600  # 3600 is for conversion from seconds to hours
plt.plot(r, P)
plt.grid()
plt.xlabel("Orbit Radius [km]")
plt.ylabel("Orbit Period [h]")
plt.show()

# Plot orbital period over altitude (up to 1200 km)

h = np.arange(0, 1200)
r = e.R + h
P = 2 * np.pi * np.sqrt(r**3/e.mu) / 60  # 60 is for conversion from seconds to minutes
plt.plot(h, P)
plt.grid()
plt.xlabel("Altitude [km]")
plt.ylabel("Orbit Period [min]")
plt.show()

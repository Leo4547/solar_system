from datetime import datetime
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import math

# Current UTC time
UTC_time = datetime.utcnow()
t = Time(UTC_time)

plt.style.use("_mpl-gallery")

# Orbital periods for ephemeris sampling
# a, e, i, w, Ω, T
orb_elements = {
    "mercury": (0.3871, 0.20564, 7.006, 77.46, 48.34, 0.241),
    "venus": (0.7233, 0.00676, 3.398, 131.77, 76.67, 0.615),
    "earth": (1.0, 0.01673, 0.0, 102.93, 0.0, 1.0),
    "mars": (1.5237, 0.09337, 1.852, 336.08, 49.71, 1.881),
    "jupiter": (5.2025, 0.04854, 1.299, 14.27, 100.29, 11.87),
    "saturn": (9.5415, 0.05551, 2.492, 92.86, 113.64, 29.47),
    "uranus": (19.188, 0.04689, 0.773, 172.43, 73.96, 84.05),
    "neptune": (30.070, 0.00895, 1.770, 46.68, 131.79, 164.9)
}

planets = list(orb_elements.keys())

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 10))

# Plot Sun with text
ax.scatter(0, 0, color="yellow", s=200)
ax.text(0, 0, "Sun", fontsize=10, ha="center", va="bottom")

def orbit_xyz(elements, num_points=500):
    # Convert degrees into radians
    # arbitrary change
    i = np.radians(elements[2])
    w = np.radians(elements[3])
    Ω = np.radians(elements[4])
    a = elements[0]
    e = elements[1]

    # Create a list of evenly spaced values from 0 to 2 * pi
    v = np.linspace(0, 2 * np.pi, num_points)

    # Radius in orbital plane
    r = (a * (1 - e ** 2)) / (1 + e * np.cos(v))

    # 3D coordinates
    x = r * (np.cos(Ω) * np.cos(w + v) - np.sin(Ω) * np.sin(w + v) * np.cos(i))
    y = r * (np.sin(Ω) * np.cos(w + v) + np.cos(Ω) * np.sin(w + v) * np.cos(i))
    z = r * (np.sin(w + v) * np.sin(i))

    return x, y, z

planets_coordinates = {}

for planet_name, planet_elements in orb_elements.items():
    print(planet_name)
    print(planet_elements)
    x_orbit, y_orbit, z_orbit = orbit_xyz(planet_elements)
    ax.plot(x_orbit, y_orbit, linestyle="--", alpha=0.6, label=f"{planet_name} orbit")

for planet in planets:

    a, e, i, w, Ω, T = orb_elements[planet]

    pos, vel = get_body_barycentric_posvel(planet, t)
    x_now = pos.x.to(u.AU).value
    y_now = pos.y.to(u.AU).value
    planets_coordinates[planet] = (x_now, y_now)

planets_x = [coord[0] for coord in planets_coordinates.values()]
planets_y = [coord[1] for coord in planets_coordinates.values()]


ax.scatter(planets_x, planets_y)

for planet, (x, y) in planets_coordinates.items():
    ax.text(x, y, planet.capitalize(), fontsize=10, ha="center")

ax.set_title("2D Solar System")
plt.legend(loc="upper right")
plt.show()
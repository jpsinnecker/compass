import numpy as np

density = 7850.0
thickness = 0.26e-3
needle_len = 10e-3
needle_width = 3e-3
pivot_radius = 1e-3
pivot_thickness = 2e-3
pivot_density = 8500.0
Ms = 1.59e6

# Blade volume & mass
area_solid = 0.5 * needle_len * needle_width
area_hole = np.pi * (pivot_radius**2)
area_net = area_solid - area_hole
volume_blade = area_net * thickness
mass_blade = density * volume_blade

# Blade moment of inertia
mass_solid = density * area_solid * thickness
I_solid = (1.0 / 24.0) * mass_solid * (needle_len**2 + needle_width**2)

mass_hole = density * area_hole * thickness
I_hole = 0.5 * mass_hole * (pivot_radius**2)

I_blade = I_solid - I_hole

# Pivot volume & mass
mass_pivot = pivot_density * np.pi * (pivot_radius**2) * pivot_thickness
I_pivot = 0.5 * mass_pivot * (pivot_radius**2)

I_total = I_blade + I_pivot
m_total = Ms * volume_blade

print(f"Blade net area = {area_net:.6e} m2")
print(f"Blade net volume = {volume_blade:.6e} m3")
print(f"Blade mass = {mass_blade:.6e} kg")
print(f"I_solid = {I_solid:.6e} kg.m2")
print(f"I_hole = {I_hole:.6e} kg.m2")
print(f"I_blade = {I_blade:.6e} kg.m2")
print(f"Pivot mass = {mass_pivot:.6e} kg")
print(f"I_pivot = {I_pivot:.6e} kg.m2")
print(f"Total Inertia I_total = {I_total:.6e} kg.m2")
print(f"Total Magnetic Moment m = {m_total:.6e} A.m2")

import math

from matplotlib import pyplot as plt
from pyXSteam.XSteam import XSteam


def average_velocity(mass_flow, density, cross_section) -> float:
    return mass_flow / (cross_section * density)

def reynolds_number(density, average_velocity, hydraulic_diameter, dynamic_viscosity) -> float:
    return density * average_velocity * hydraulic_diameter / dynamic_viscosity

def prantll_number(dynamic_viscosity, specific_heat, thermal_conductivity) -> float:
    return dynamic_viscosity * specific_heat / thermal_conductivity

def nusselt_number(reynolds_number, prantll_number) -> float:
    return 0.023 * reynolds_number ** 0.8 * prantll_number** 0.4

def heat_transfer_coefficient(nusselt_number, hydraulic_diameter, thermal_conductivity) -> float:
    return nusselt_number * thermal_conductivity / hydraulic_diameter

def heat_flux(heat_transfer_coefficient, surface_temp, bulk_fluid_temp, hydraulic_diameter, pipe_length) -> float:
    area = pipe_length * (math.pi * hydraulic_diameter/2 )**2
    return heat_transfer_coefficient * (surface_temp - bulk_fluid_temp) * area

def outlet_temp(heat_flux, mass_flow, specific_heat, inlet_temp) -> float:
    return heat_flux / (mass_flow * specific_heat) + inlet_temp

def calculate_node(
    mass_flow,
    hydraulic_diameter,
    pipe_diameter,
    surface_temp,
    inlet_temp,
    tube_pressure,
    condenser_pressure,
    node_length
) -> tuple:
    
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    density=steam_table.rho_pt(tube_pressure, inlet_temp)
    dynamic_viscosity=steam_table.my_pt(tube_pressure, inlet_temp)
    specific_heat=steam_table.Cp_pt(tube_pressure, inlet_temp)
    thermal_conductivity=steam_table.tc_pt(tube_pressure, inlet_temp)
    surface_temp=steam_table.tsat_p(condenser_pressure)

    cross_section = math.pi * pipe_diameter ** 2 / 4

    avg_vel = average_velocity(
        mass_flow=mass_flow,
        density=density,
        cross_section=cross_section
    )

    reynolds = reynolds_number(
        density=density,
        average_velocity=avg_vel,
        hydraulic_diameter=hydraulic_diameter,
        dynamic_viscosity=dynamic_viscosity
    )

    prantll = prantll_number(
        dynamic_viscosity=dynamic_viscosity,
        specific_heat=specific_heat,
        thermal_conductivity=thermal_conductivity
    )

    nusselt = nusselt_number(reynolds_number=reynolds, prantll_number=prantll)

    heat_transfer_coeff = heat_transfer_coefficient(
        nusselt_number=nusselt, 
        hydraulic_diameter=hydraulic_diameter,
        thermal_conductivity=thermal_conductivity
    )

    flux = heat_flux(
        heat_transfer_coefficient=heat_transfer_coeff,
        surface_temp=surface_temp,
        bulk_fluid_temp=inlet_temp,
        hydraulic_diameter=hydraulic_diameter,
        pipe_length=node_length
    )

    outlet= outlet_temp(
        heat_flux=flux,
        specific_heat=specific_heat,
        mass_flow=mass_flow,
        inlet_temp=inlet_temp
    )

    return outlet

def calculate_pipe_outlet_temp(
    n_nodes,
    pipe_length,
    mass_flow,
    pipe_diameter,
    hydraulic_diameter,
    surface_temp,
    inlet_temp,
    tube_pressure,
    condenser_pressure
) -> float:
    """Calculates pipe outlet temperature by 
    splitting up the pipe into n_nodes and calculating each
    node's outlet temperature based on previous
    node's outlet temperature. 
    """
    node_length = pipe_length / n_nodes
    for _ in range(n_nodes):
        inlet_temp = calculate_node(
            mass_flow=mass_flow,
            hydraulic_diameter=hydraulic_diameter,
            pipe_diameter=pipe_diameter,
            surface_temp=surface_temp,
            inlet_temp=inlet_temp,
            tube_pressure=tube_pressure,
            condenser_pressure=condenser_pressure,
            node_length=node_length
        )
    return inlet_temp

# PROBLEM 1

pipe_diameter = 0.01587 #meters
init_fluid_temp = 15 #Celsius
init_pressure = 4 #bar
mass_flow = 0.5 #kg/s per tube
surface_temp = 344.79 #Celsius
hydraulic_diameter = 0.01587 #meters
condenser_pressure=0.474 #bar
length=8

#Calculate pipe outlet temperature for increasing number of nodes
n_nodes = list(range(1, 100))
outlet_temps = [
    calculate_pipe_outlet_temp(
        n_nodes = n,
        pipe_length=length,
        mass_flow=mass_flow,
        pipe_diameter=pipe_diameter,
        hydraulic_diameter=hydraulic_diameter,
        surface_temp=surface_temp,
        inlet_temp=init_fluid_temp,
        tube_pressure=init_pressure,
        condenser_pressure=condenser_pressure
    )
    for n in n_nodes
]

print(f"Pipe outlet temperature is {outlet_temps[-1]:.2f} Celsius with {n_nodes[-1]} nodes with length {length/n_nodes[-1]:.4f}m")

#Plotting to see convergence

fig, ax = plt.subplots(figsize=(8, 8))
ax.plot(n_nodes, outlet_temps)
ax.set_title("Pipe Outlet Temperature Simulation")
ax.set_xlabel("Number of Nodes")
ax.set_ylabel("Temperature (Celsius)")
ax.axhline(outlet_temps[-1] - 0.2, xmin=0, xmax=n_nodes[-1], color="red", linestyle="dashed")
ax.spines['top'].set_visible(False)
ax.spines["right"].set_visible(False)
ax.text(n_nodes[-5], outlet_temps[-1] + 2, f"Pipe outlet temp -> {outlet_temps[-1]:.2f}C")
fig.savefig("problem_1.png")
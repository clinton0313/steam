#%%
import math

from pyXSteam.XSteam import XSteam

steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)
#%%

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
    bulk_fluid_temp,
    inlet_temp,
    tube_pressure
) -> tuple:
    
    steam_table = XSteam(XSteam.UNIT_SYSTEM_MKS)

    density=steam_table.rho_pt(tube_pressure, inlet_temp)
    dynamic_viscosity=steam_table.my_pt(tube_pressure, inlet_temp)
    specific_heat=steam_table.Cp_pt(tube_pressure, inlet_temp)
    thermal_conductivity=steam_table.tc_pt(tube_pressure, inlet_temp)
    enthalpy = steam_table.h_pt(tube_pressure, inlet_temp)

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
        bulk_fluid_temp=bulk_fluid_temp,
        hydraulic_diameter=hydraulic_diameter,
        pipe_length=pipe_length
    )

    outlet= outlet_temp(
        heat_flux=flux,
        specific_heat=specific_heat,
        mass_flow=mass_flow,
        inlet_temp=inlet_temp
    )

    return outlet
# %%
n_pipes = 1000
pipe_length = 1 #meters
pipe_diameter = 0.01587 #meters
init_fluid_temp = 15 #Celsius
init_pressure = 4 #bar
mass_flow = 0.5 #kg/s per tube
surface_temp = 344.79 #Celsius
hydraulic_diameter = 0.01587 #meters
bulk_fluid_temp=15


node_1 = calculate_node(
    mass_flow=mass_flow,
    hydraulic_diameter=hydraulic_diameter,
    pipe_diameter=pipe_diameter,
    surface_temp=surface_temp,
    bulk_fluid_temp=bulk_fluid_temp,
    inlet_temp=init_fluid_temp,
    tube_pressure=init_pressure

)

print(f"Node 1 outlet temp: {node_1:.2f}")
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

def fermi_dirac_velocity_distribution(v_norm, T_norm):
    """
    Berekent de genormaliseerde snelheidsverdeling.
    v_norm: v / v_F (snelheid gedeeld door Fermi-snelheid)
    T_norm: k_B * T / E_F (gereduceerde temperatuur)
    """
    # Voorkom deling door nul bij het absolute nulpunt
    if T_norm == 0:
        T_norm = 1e-10 
        
    # De toestandsdichtheid schaalt met v^2
    density_of_states = v_norm**2
    
    # De Fermi-Dirac kansfunctie, vertaald naar snelheid
    # (v^2 - 1) komt van (E - E_F) / E_F
    x = np.clip((v_norm**2 - 1) / T_norm, -700, 700)
    probability = 1 / (np.exp(x) + 1)
    
    return density_of_states * probability

# Stel de X-as in (genormaliseerde snelheid van 0 tot 1.5)
v = np.linspace(0, 1.5, 500)

# Definieer een paar verschillende temperaturen om te vergelijken
temperatures = [0.001, 0.05, 0.15]
colors = ['blue', 'orange', 'red']

plt.figure(figsize=(10, 6))

for T, color in zip(temperatures, colors):
    n_v = fermi_dirac_velocity_distribution(v, T)
    plt.plot(v, n_v, label=f'T/T_F = {T}', color=color, linewidth=2)
    plt.fill_between(v, n_v, alpha=0.1, color=color)

# Grafiek opmaken
plt.title('Snelheidsverdeling in een Elektronengas (Fermi-Dirac)')
plt.xlabel('Genormaliseerde snelheid ($v / v_F$)')
plt.ylabel('Waarschijnlijkheidsdichtheid $n(v)$')
plt.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, label='Fermi-snelheid ($v_F$)')
plt.xlim(0, 1.5)
plt.ylim(0, 1.2)
plt.legend()
plt.grid(True, alpha=0.3)

output_path = Path(__file__).with_name("fermi_dirac_plot.png")
plt.savefig(output_path, dpi=200, bbox_inches="tight")

if os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"):
    plt.show()
else:
    print(f"Geen grafische display gevonden. Plot opgeslagen in: {output_path}")

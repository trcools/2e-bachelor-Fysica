import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# --- 1. KOSMOLOGISCHE SETUP (De schaal van het heelal bepalen) ---
# We gebruiken H0 = 70 km/s/Mpc en een standaard dichtheid van materie (Omega_m0 = 0.3)
# Geef een expliciete type-annotatie en onderdruk hier statische type-waarschuwingen van Pylance/pyright.
cosmo: FlatLambdaCDM = FlatLambdaCDM(H0=70*u.km/u.s/u.Mpc, Om0=0.3)  # type: ignore[call-arg]

z_L = 0.444  # Redshift van de lens
z_S = 2.379  # Redshift van de bron

# Bereken de hoekdiameter-afstanden in Megaparsec (Mpc)
D_L = cosmo.angular_diameter_distance(z_L)  # type: ignore[attr-defined]
D_S = cosmo.angular_diameter_distance(z_S)  # type: ignore[attr-defined]
D_LS = cosmo.angular_diameter_distance_z1z2(z_L, z_S)  # type: ignore[attr-defined]

# Bereken de 'Time Delay Distance' (D_dt)
D_dt = (1 + z_L) * (D_L * D_S / D_LS)

# Omrekenningsfactor: D_dt is in Mpc, we willen de tijdvertraging in dagen
# 1 Mpc = 3.086e22 meter, c = 299792458 m/s, 1 dag = 86400 seconden
# Hoek moet van boogseconden (arcsec) naar radialen! (1 arcsec = 4.848e-6 rad)
arcsec_to_rad = np.pi / (180 * 3600)
c = 299792458 * u.m / u.s

time_delay_factor = (D_dt / c).to(u.day).value * (arcsec_to_rad ** 2)

print(f"Time Delay Distance (D_dt): {D_dt:.2f}")
print(f"1 arcsec^2 verschil in potentiaal komt overeen met {time_delay_factor:.2f} dagen.")

# --- 2. LENS MODEL PARAMETERS (Evans & Witt 2003) ---
# Fourier coëfficiënten voor het linker paneel (4 beelden: A, B, C, D)
a0 = 9.89
a2 = 0.090
b2 = -0.11
a3 = 0.02
b3 = -0.04

# Positie van de bron (beta). 
# We zetten deze iets uit het midden om de waargenomen asymmetrie (A, B, C, D) te genereren.
#beta_x = 0.2
#beta_y = 0.5 
beta_x = 0.1
beta_y = 0.1 

# --- 3. HET RASTER (GRID) MAKEN ---
# We bekijken een gebied van -6 tot +6 boogseconden (zoals de assen in Figuur 5)
grid_size = 5000
x = np.linspace(-6, 6, grid_size)
y = np.linspace(-6, 6, grid_size)
X, Y = np.meshgrid(x, y)

# Zet om naar poolcoördinaten (radius r en hoek theta)
R = np.sqrt(X**2 + Y**2)
Theta = np.arctan2(Y, X)

# --- 4. DE FERMAT POTENTIAAL BEREKENEN ---
# Stap A: De vorm van de lens F(theta) berekenen via de Fourierreeks
# (Bij dit specifieke type model stelt a0/2 ruwweg de Einstein-straal voor)
#F_theta = (a0 / 2) + (1**2-2**2)*a2*np.cos(2*Theta) + (1**2-2**2)*b2*np.sin(2*Theta) + (1**2-3**2)*a3*np.cos(3*Theta) + (1**2-3**2)*b3*np.sin(3*Theta)
F_theta = (a0 / 2) + a2*np.cos(2*Theta) + b2*np.sin(2*Theta) + a3*np.cos(3*Theta) + b3*np.sin(3*Theta)

# Stap B: De zwaartekrachtspotentiaal (psi) van de lens
Psi = R * F_theta

# Stap C: De geometrische vertraging
# (Hoeveel langer is het pad vergeleken met een rechte lijn?)
Geom_delay = 0.5 * ((X - beta_x)**2 + (Y - beta_y)**2)

# Stap D: De totale Fermat Potentiaal (Geometrie - Zwaartekracht)
# Let op: de eenheid hier is in (boogseconden)^2
Tau = Geom_delay - Psi

# --- Waargenomen Positie maxima van luminositeit ---
# Er zijn 4 maxima in Tau die overeenkomen met de 4 beelden (A, B, C, D).
# A: (3".0,4".6), B: (-1".1,5".2), C: (-4".7,2".2) en D: (2."0,-4".0)
# We kunnen deze punten markeren op de plot van Tau om te laten zien waar de beelden zich bevinden in het 'landschap' van de Fermat potentiaal.
image_positions = np.array([[3.0, 4.6], [-1.1, 5.2], [-4.7, 2.2], [2.0, -4.0]])  # in arcsec
image_indices = []
 
# Als je `image_positions` in RA-conventie hebt opgegeven (positief naar links),
# zet `positions_in_ra_convention = True`. De plot gebruikt standaard positieve x naar rechts
# en de x-coördinaten van de markers worden automatisch gespiegeld.
positions_in_ra_convention = True

# De fouten op de posities zijn =< 0".4, dus we kunnen een cirkel van 0".4 rond elk punt tekenen om de onzekerheid aan te geven.
# We kunnen ook de deze foutenvlaggen op de punten markeren op de plot van tau, bijvoorbeeld met een cirkel of een kruis.
error_radius = 0.4  # in arcsec
theta_circle = np.linspace(0, 2*np.pi, 100)





# --- 5. TIJDVERTRAGINGEN BEREKENEN (OPTIONEEL) ---
# We kunnen ook de tijdvertragingen tussen verschillende punten in het raster berekenen
# Bijvoorbeeld, de tijdvertraging tussen het centrum (0,0) en een punt (x,y) is:
# time_delay = time_delay_factor * (Tau - Tau_center)
# Waar Tau_center = Tau op (0,0) is. Dit geeft ons een idee van de tijdvertragingen in dagen tussen verschillende delen van het 'landschap' van de Fermat potentiaal.
Tau_center = Tau[grid_size//2, grid_size//2]  # Tau op het midden van het raster
Time_Delay_Map = time_delay_factor * (Tau - Tau_center)  # In dagen 





# --- 6. PLOTTEN EN VISUALISEREN ---
plt.figure(figsize=(8, 8))
plt.title("Gereconstrueerde Fermat Potentiaal (Cosmic Horseshoe)", fontsize=14)

# We tekenen de contouren, vergelijkbaar met Figuur 5 in de paper
# We gebruiken 50 niveaus om het 'landschap' goed in kaart te brengen
levels = np.linspace(np.min(Tau), np.max(Tau), 30)
contour = plt.contour(X, Y, Tau, levels=levels, colors='black', linewidths=0.8)

# Assen netjes maken zoals in de publicatie
plt.xlabel("$\\Delta\\alpha$ [arcsec]", fontsize=12)
plt.ylabel("$\\Delta\\delta$ [arcsec]", fontsize=12)
# Laat de standaard-orientatie (positieve x naar rechts) zien.
# Als je de RA-conventie wilt (Oost links), zet je de volgende regel terug:
# plt.gca().invert_xaxis()
plt.xlim(-6, 6)
plt.ylim(-6, 6)

for pos in image_positions:
    plot_x = -pos[0] if positions_in_ra_convention else pos[0]
    idx_x = np.argmin(np.abs(x - plot_x))
    idx_y = np.argmin(np.abs(y - pos[1]))
    image_indices.append((idx_y, idx_x))  # (row, col) in de Tau array 

    circle_x = plot_x + error_radius * np.cos(theta_circle)
    circle_y = pos[1] + error_radius * np.sin(theta_circle)
    plt.plot(circle_x, circle_y, color='blue', linestyle='--', linewidth=0.5)  # Foutcirkel
    plt.plot(plot_x, pos[1], 'x', color='blue')  # Positie van het beeld
    

# De Kritische Kromme (De rode lijn) tekenen
# Bij dit specifieke type lens ligt de kromme ongeveer op de plekken waar r = F(theta)
x_crit = F_theta * np.cos(Theta)
y_crit = F_theta * np.sin(Theta)
# (Om hem netjes te tekenen gebruiken we een simpele lijnplot voor r = F(theta) ipv het grid)
theta_line = np.linspace(0, 2*np.pi, 1000)
f_line = (a0 / 2) + a2*np.cos(2*theta_line) + b2*np.sin(2*theta_line) + a3*np.cos(3*theta_line) + b3*np.sin(3*theta_line)
plt.plot(f_line*np.cos(theta_line), f_line*np.sin(theta_line), color='red', label="Kritische Kromme")

plt.legend()
plt.grid(False)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
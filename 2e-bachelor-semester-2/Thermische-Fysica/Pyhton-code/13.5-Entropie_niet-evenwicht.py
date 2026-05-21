# imports
import numpy as np
import matplotlib.pyplot as plt

# We hebben een staaf met lengte L.
# We veronderstellen de begintoestand stationair.

# beginconties
T_i0 = 100 # Temperatuur aan het begin van de staaf
T_iL = 0.1 # Temperatuur aan het einde van de staaf
length = 1 # Lengte van de staaf
position = np.linspace(0, length, 100) # Positie langs de staaf

#De begintemperatuur op elke plaats x van de staaf is:
def initial_temperature(x=position, L=length, T_i0=T_i0, T_iL=T_iL):
    return T_i0 - (T_i0-T_iL)/L * x

# Voorstelling
plt.plot(position, initial_temperature(position), label='Begintoestand')
plt.xlabel('Positie langs de staaf')
plt.ylabel('Temperatuur')
plt.title('Begintoestand van de temperatuur langs de staaf')
plt.legend()
plt.grid()
plt.show()

# De plot is lineair


# ==============================================
# Entropiebepaling
# ==============================================

warmte_capaciteit = 1 # Warmtecapaciteit van het materiaal
massdichtheid = 1 # Massadichtheid van het materiaal
oppervlakte = 1 # Oppervlakte van de staaf

def entropy_niet_evenwicht(T_0=T_i0, T_L=T_iL, T_f=T_iL, L=length, c=warmte_capaciteit, rho=massdichtheid, A=oppervlakte):
    # Integraal: int_0^L ln(T_0/T_f - (T_0-T_L)/(L*T_f)*x) dx
    def integrand(x):
        return np.log(T_0/T_f - (T_0 - T_L)/(L * T_f) * x)
    
    # Numerieke integratie
    from scipy.integrate import quad
    integral, error = quad(integrand, 0, L)
    
    # Bereken de verandering in entropie
    delta_S = -c * rho * A * integral
    return delta_S

# Standalone interactive plot using matplotlib sliders (works outside Jupyter)
from scipy.integrate import quad
from matplotlib.widgets import Slider, Button


def plot_with_sliders():
    fig, ax = plt.subplots(figsize=(8, 5))
    plt.subplots_adjust(left=0.15, bottom=0.30)

    x = np.linspace(0, length, 400)

    # initial values
    T0_val = T_i0
    TL_val = T_iL
    Tf_val = T_iL

    def integrand_vals(T0, TL, Tf):
        return np.log(T0 / Tf - (T0 - TL) / (length * Tf) * x)

    y = integrand_vals(T0_val, TL_val, Tf_val)
    l, = ax.plot(x, y, lw=2)
    fill = ax.fill_between(x, y, alpha=0.3)
    ax.set_xlabel('Positie langs de staaf')
    ax.set_ylabel(r'$\ln\left(\frac{T_0}{T_f} - \frac{T_0-T_L}{L T_f} x\right)$')
    ax.set_title('Integrand voor de verandering in entropie')
    ax.grid(True)

    # text box for delta S
    txt = ax.text(0.02, 0.95, '', transform=ax.transAxes, va='top')

    # Slider axes
    axcolor = 'lightgoldenrodyellow'
    ax_T0 = plt.axes([0.15, 0.20, 0.70, 0.03], facecolor=axcolor)
    ax_TL = plt.axes([0.15, 0.14, 0.70, 0.03], facecolor=axcolor)
    ax_Tf = plt.axes([0.15, 0.08, 0.70, 0.03], facecolor=axcolor)

    sT0 = Slider(ax_T0, 'T0', 50.0, 150.0, valinit=T0_val, valstep=1.0)
    sTL = Slider(ax_TL, 'TL', 0.0, 150.0, valinit=TL_val, valstep=1.0)
    sTf = Slider(ax_Tf, 'Tf', 0.01, 200.0, valinit=max(Tf_val, 0.01), valstep=1.0)

    def update(val):
        T0 = sT0.val
        TL = sTL.val
        Tf = sTf.val
        # avoid invalid values that would make log arg <= 0
        # guard against Tf <= 0 which causes division by zero
        if Tf <= 0:
            Tf = 1e-8
        arr = integrand_vals(T0, TL, Tf)
        valid = np.isfinite(arr) & (np.isreal(arr))
        l.set_ydata(arr)
        # remove previous collections (robust across matplotlib versions)
        for coll in list(ax.collections):
            try:
                coll.remove()
            except Exception:
                pass
        ax.fill_between(x, arr, alpha=0.3)

        # compute delta S safely using quad on continuous integrand
        try:
            def integrand_scalar(xi):
                # guard against division by zero
                Tf_local = Tf if Tf > 0 else 1e-8
                return np.log(T0 / Tf_local - (T0 - TL) / (length * Tf_local) * xi)

            integral, err = quad(integrand_scalar, 0, length)
            delta_S = -warmte_capaciteit * massdichtheid * oppervlakte * integral
            txt.set_text(f'ΔS = {delta_S:.3f} J/K')
        except Exception as e:
            txt.set_text('ΔS: invalid params')

        fig.canvas.draw_idle()

    sT0.on_changed(update)
    sTL.on_changed(update)
    sTf.on_changed(update)

    # Reset button
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    def reset(event):
        sT0.reset()
        sTL.reset()
        sTf.reset()

    button.on_clicked(reset)

    # initial update to compute and show ΔS
    update(None)
    plt.show()


if __name__ == '__main__':
    plot_with_sliders()

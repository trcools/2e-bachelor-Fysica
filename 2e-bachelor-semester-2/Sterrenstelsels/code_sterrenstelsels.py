import numpy as np
import matplotlib.pyplot as plt


def show_parsec():
    # Parsec: afstand waarbij 1 AU een parallaxhoek van 1 boogseconde heeft.
    # Dit visualiseert het geometrische verband: p (in boogseconden) en d (in parsec).

    fig, ax = plt.subplots(figsize=(9, 4))

    # Basislijn 1 AU, afstand 1 pc
    d_pc = 1
    p_arcsec = 1

    # Gebruik kleine-hoek benadering: p (rad) = 1 AU / d
    # Voor visualisatie nemen we een grote schaal in de tekening.
    baseline = 1.0  # 1 AU (genormaliseerd)
    distance = 12.0  # schaal voor visueel effect
    p_rad = np.deg2rad(p_arcsec / 3600)

    # Observatiepunten en sterpositie (horizontale layout)
    sun = np.array([0, 0])
    earth1 = np.array([0, -baseline])
    earth2 = np.array([0, baseline])
    star = np.array([distance, 0])

    # Lijnen van Aarde naar ster
    ax.plot([earth1[0], star[0]], [earth1[1], star[1]], color="tab:blue", lw=2.5)
    ax.plot([earth2[0], star[0]], [earth2[1], star[1]], color="tab:blue", lw=2.5)

    # Middenlijn (referentie)
    ax.plot([0, star[0]], [0, star[1]], color="tab:gray", lw=1.2, ls="--", alpha=0.7)

    # Basislijn (1 AU)
    ax.plot([earth1[0], earth2[0]], [earth1[1], earth2[1]], color="black", lw=2.2)
    ax.annotate(
        "1 AU",
        xy=(sun[0], sun[1]),
        xytext=(-1.2, baseline / 2),
        textcoords="data",
        ha="right",
        va="center",
        fontsize=11,
        arrowprops=dict(arrowstyle="-[", lw=1.4, color="black"),
    )

    # Parallaxhoek p
    angle_radius = 1.8
    theta = np.linspace(-p_rad, p_rad, 80)
    ax.plot(angle_radius * np.cos(theta), angle_radius * np.sin(theta), color="tab:red", lw=2.2)
    ax.text(angle_radius + 0.25, 0.25, "p = 1\"", color="tab:red", fontsize=11)
    ax.annotate(
        "parallaxhoek",
        xy=(angle_radius, 0),
        xytext=(angle_radius + 1.1, 1.0),
        color="tab:red",
        fontsize=10,
        arrowprops=dict(arrowstyle="->", lw=1.2, color="tab:red"),
    )

    # Labels
    ax.scatter([sun[0]], [sun[1]], color="gold", s=140, label="Zon")
    ax.scatter([earth1[0], earth2[0]], [earth1[1], earth2[1]], color="tab:green", s=80, label="Aarde (6 maanden apart)")
    ax.scatter([star[0]], [star[1]], color="tab:orange", s=120, label="Ster")
    ax.text(sun[0] - 0.2, sun[1] - 0.4, "Zon", fontsize=10, ha="right")
    ax.text(star[0] + 0.4, star[1], "Ster", va="center", fontsize=11)
    ax.text(0.2, earth2[1] + 0.3, "Aarde (juni)", fontsize=10)
    ax.text(0.2, earth1[1] - 0.6, "Aarde (dec)", fontsize=10)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-2, 13)
    ax.set_ylim(-4, 4)
    ax.axis("off")
    ax.set_title("Definitie van 1 parsec (parallax 1\")", fontsize=12)
    ax.legend(loc="upper right", frameon=True, fontsize=10)

    plt.show()
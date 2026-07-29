# -*- coding: utf-8 -*-
#
#     ╔═══════════════════════════════════════════════════════════╗
#     ║ TIPE 2025-2026, MPI Lycée Louis Thuillier: Cardot Clément ║
#     ╚═══════════════════════════════════════════════════════════╝
#
# ══════════════════════ Objectif de ce code: ══════════════════════
# Permettre le tracé des composantes x(t) et y(t) du contour correspondant au fichier echantillonage.csv genéré par Calcul_2.c
# et la comparaison avec la reconstruction par les coefficients de Fourier calculés par Coeffs_Fourrier.c
#
# ════════════════════ Bibliothèques utilisées: ════════════════════
#    - matplotlib: pour le tracé des graphiques.
#    - numpy: pour la manipulation efficace des tableaux de données.
#  
# ═══════════════════════ Structure du code: ═══════════════════════
#       1. Lecture et stockage des données du fichier CSV.
#       2. Création de deux graphiques:
#           a. Un graphique pour le contour (x,y) avec les points Moore Neighbour et la reconstruction par Fourier.
#           b. Un graphique pour l'évolution de x(t) et y(t) en fonction de t, pour les deux méthodes.
#       3. Affichage des courbes.

print("     ╔═══════════════════════════════════════════════════════════╗\n     ║ TIPE 2025-2026, MPI Lycée Louis Thuillier: Cardot Clément ║\n     ╚═══════════════════════════════════════════════════════════╝")
print("\n ⟶ Objectif: Plot des composantes x(t) et y(t) à partir des données du fichier \n    echantillonage.csv et comparaison avec la reconstruction par Fourier.\n")


import matplotlib.pyplot as plt
import numpy as np

# ══════════════════════ Création d'un créneau ═════════════════════
# N = 1000
# periode = 0.2
# phase = 0.05
#
# with open("echantillonage.csv", "w") as f:
#     f.write("t,x,y\n")
#
#     for k in range(N):
#         t = k / N
#
#         # Créneau x(t)
#         if (t % periode) < (periode / 2):
#             x = 200
#         else:
#             x = 0
#
#         # Créneau y(t) déphasé
#         if ((t + phase) % periode) < (periode / 2):
#             y = 300
#         else:
#             y = 0
#
#         f.write(f"{t:.3f},{x},{y}\n")


# ══════ Lecture echantillonage.csv (points Moore Neighbour) ═══════
t_points, x_points, y_points = [], [], []

with open("echantillonage.csv") as f:
    next(f)  # sauter l'en-tête
    for line in f:
        line = line.strip()
        if not line:
            continue  # ignore ligne vide
        
        parts = line.split(",")
        if len(parts) != 3:
            print("Ligne invalide :", line)
            continue

        t_val, x_val, y_val = parts
        t_points.append(float(t_val))
        x_points.append(float(x_val))
        y_points.append(float(y_val))

t_points = np.array(t_points)
x_points = np.array(x_points)
y_points = np.array(y_points)
print(" ",len(t_points), "points lus dans echantillonage.csv.")


# ═══════════ Lecture coeffs_fourier.csv (coefficients) ════════════
n_vals, re_vals, im_vals = [], [], []

with open("coeffs_fourier.csv") as f:
    next(f)  # sauter l'en-tête
    for line in f:
        n, re, im = line.strip().split(",")
        n_vals.append(int(n))
        re_vals.append(float(re))
        im_vals.append(float(im))

n_vals = np.array(n_vals)
c_vals = np.array(re_vals) + 1j*np.array(im_vals)

nb_coeffs = len(n_vals)
print(" ",nb_coeffs, "coefficients lus dans coeffs_fourier.csv.")


# ═══════════════════ Reconstruction par Fourier ═══════════════════
def reconstruire(t, n_vals, c_vals):
    """Calcule f(t) = somme c_n * exp(i*2pi*n*t)."""
    somme = np.zeros_like(t, dtype=complex)
    for n, c in zip(n_vals, c_vals):
        somme += c * np.exp(1j * 2*np.pi * n * t)
    return somme

t_reconstruit = np.linspace(0, 1, 2000)
f_t = reconstruire(t_reconstruit, n_vals, c_vals)
x_reconstruit = f_t.real
y_reconstruit = f_t.imag
print("  Reconstruction par Fourier effectuée.")


# ════════════════ Tracés avec proportions adaptées ════════════════
fig, axes = plt.subplots(
    3, 1, figsize=(10, 12),
    gridspec_kw={"height_ratios": [2, 1, 1]}  # plot (x,y) plus grand
)

# (1) Courbe des composantes (x,y) : Fourier vs Moore Neighbour
axes[0].scatter(x_points, y_points, s=0.7, color="gray", label="Points (Moore Neighbour)")
axes[0].plot(x_reconstruit, y_reconstruit, color="blue", linewidth=1.2, alpha=0.7, label="Reconstruction Fourier")
axes[0].set_aspect("equal")
axes[0].invert_yaxis()
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
axes[0].set_title("Contour (x,y) : Moore Neighbour vs Fourier")
axes[0].legend(loc="center left", bbox_to_anchor=(1, 0.5))
axes[0].grid(True, alpha=0.3)

# (2) Courbes x(t), y(t) Moore Neighbour
axes[1].scatter(t_points, x_points, s=0.5, c='red', label='x(t)')
axes[1].scatter(t_points, y_points, s=0.5, c='green', label='y(t)')
axes[1].set_title("Composantes x(t) et y(t) (Moore Neighbour Tracing)")
axes[1].set_xlabel("t")
axes[1].set_ylabel("Valeurs")
axes[1].legend(fontsize=9, markerscale=5, loc="upper right")
axes[1].grid(True, alpha=0.5)

# (3) Courbes x(t), y(t) par Fourier
axes[2].plot(t_reconstruit, x_reconstruit, c='red', linewidth=1.4, label='x(t)')
axes[2].plot(t_reconstruit, y_reconstruit, c='green', linewidth=1.4, label='y(t)')
axes[2].set_title(rf"Reconstruction par Fourier avec {nb_coeffs} coefficients : "
                  r"$f(t)=\sum_{n=-M}^M c_n e^{i 2\pi n t}$, "
                  r"$x(t)=Re(f(t)),\; y(t)=Im(f(t))$")
axes[2].set_xlabel("t")
axes[2].set_ylabel("Valeurs")
axes[2].legend(fontsize=9, loc="upper right")
axes[2].grid(True, alpha=0.4)

plt.tight_layout()
print("  Affichage en cours.")
plt.show()

print("\n  Fermeture du programme.\n")
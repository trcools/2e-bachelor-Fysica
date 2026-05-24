# 1. Eindige Groepen (Hoofdstuk 2)
from sympy.combinatorics import DihedralGroup

G = DihedralGroup(3) # Symmetriegroep van een gelijkzijdige driehoek
print("Is de groep Abels?:", G.is_abelian)
print("Aantal elementen (orde):", G.order())

# 2. Clebsch-Gordan Coëfficiënten (Hoofdstuk 4.6)
from sympy.physics.quantum.cg import CG

# Bereken de CG coëfficiënt voor j1=1, m1=1 en j2=1, m2=0 -> J=2, M=1
cg = CG(j1=1, m1=1, j2=1, m2=0, j3=2, m3=1)
print("Clebsch-Gordan waarde:", cg.doit())


from scipy.spatial.transform import Rotation as R
import numpy as np

# Maak een rotatie aan (bijv. 90 graden rond de z-as)
rot = R.from_euler('z', 90, degrees=True)

# Bekijk de matrix representatie (Hoofdstuk 4)
print("Rotatiematrix SO(3):\n", rot.as_matrix())
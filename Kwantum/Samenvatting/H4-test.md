# kwantum-samenvatting-H4

> Bronbestand: `kwantum-samenvatting-H4.pdf`

---

## 4.1 General formulae

**TDSE**  
\[
i\hbar\,\frac{\partial}{\partial t}\Psi(x,t)=
\left[
-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}
+V(x)
\right]\Psi(x,t)
\]

**TISE**  
\[
\left[
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
+V(x)
\right]\psi(x)=
E\,\psi(x)
\]

**Probability current (1D)**  
\[
j(x,t)=\frac{\hbar}{2mi}
\left[
\Psi^*(x,t)\frac{\partial\Psi(x,t)}{\partial x}
-\Psi(x,t)\frac{\partial\Psi^*(x,t)}{\partial x}
\right]
\]

---

## 4.2 The free particle

Neem \(V(x)=0\).

\[
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}=E\psi
\quad\Longrightarrow\quad
\frac{d^2\psi}{dx^2}+k^2\psi=0
\]
met
\[
E=\frac{\hbar^2k^2}{2m},
\qquad
k=\sqrt{\frac{2mE}{\hbar^2}}.
\]

**General solution**
\[
\psi(x)=A e^{ikx}+B e^{-ikx}.
\]

**Interpretation**
- \(A e^{ikx}\): vrij deeltje met impuls in de \(+\)x-richting  
- \(B e^{-ikx}\): vrij deeltje met impuls in de \(-\)x-richting

**Stationary state**
\[
\Psi(x,t)=\psi(x)\,e^{-iEt/\hbar},
\qquad
\omega=\frac{E}{\hbar}.
\]

---

## Special cases (plane waves)

### Geval 1: \(B=0\) (vrije deeltje in de positieve x-richting)
\[
\psi(x)=A e^{ikx}
\]
\[
\rho(x,t)=|\Psi|^2=|A|^2=\rho_0
\]
\[
j=\frac{\hbar k}{m}|A|^2>0
\]
\[
j=v\rho
\quad\text{met}\quad
v=\frac{\hbar k}{m}.
\]

### Geval 2: \(A=0\) (analoge situatie, negatieve x-richting)
\[
\psi(x)=B e^{-ikx}
\]
\[
\rho=|B|^2
\]
\[
j=-\frac{\hbar k}{m}|B|^2<0.
\]

### Geval 3: \(A=B\) (twee golven met gelijke amplitude in tegengestelde richting)
\[
\psi(x)=A(e^{ikx}+e^{-ikx})=2A\cos(kx)
\]
**nodes**
\[
x_m=\frac{(2m+1)\pi}{2k},\qquad m=0,1,2,\dots
\]
\[
\rho(x)=|\Psi|^2=4|A|^2\cos^2(kx)
\]
\[
j=0
\quad(\text{want } \psi \text{ reëel}).
\]

### Geval 4: \(A=-B\)
\[
\psi(x)=A(e^{ikx}-e^{-ikx})=2iA\sin(kx)
\]
\[
\rho(x)=4|A|^2\sin^2(kx),
\qquad
j=0.
\]

---

## Algemene superpositie \(A,B\neq 0\)

\[
\psi(x)=A e^{ikx}+B e^{-ikx}
\]

\[
\rho(x)=|\Psi|^2
=|A|^2+|B|^2 + A B^*e^{2ikx}+A^*B e^{-2ikx}
\]

\[
j=\frac{\hbar k}{m}\left(|A|^2-|B|^2\right)
\]

*(Interpretatie in scattering: amplitude \(A\) = “incoming” naar rechts, amplitude \(B\) = “reflected” naar links.)*

---

## Vrij deeltje in een 1D doos (periodieke randvoorwaarden)

In de praktijk wordt de beweging van het vrije deeltje vaak beperkt tot 1D door muren.

Neem een doos met lengte \(L\) en leg **periodieke randvoorwaarden** op:
\[
\psi(0)=\psi(L).
\]

Dan:
\[
k_n=\frac{2\pi n}{L},\qquad n=0,\pm1,\pm2,\dots
\]
\[
p_n=\hbar k_n=\frac{2\pi\hbar n}{L}
\]
\[
E_n=\frac{\hbar^2k_n^2}{2m}
=\frac{2\pi^2\hbar^2}{mL^2}\,n^2.
\]

Genormaliseerde eigenfuncties:
\[
\psi_n(x)=\frac{1}{\sqrt{L}}e^{ik_n x},
\qquad
\int_0^L|\psi_n(x)|^2dx=1.
\]

---

## 4.4 The potential barrier

**Potential**
\[
V(x)=
\begin{cases}
0, & x<0\\
V_0, & 0<x<a\\
0, & x>a
\end{cases}
\]

**Region I** \((x<0,\ V=0)\)
\[
\psi_I(x)=A e^{ikx}+B e^{-ikx},
\qquad
k=\sqrt{\frac{2mE}{\hbar^2}}.
\]

**Region II** \((0<x<a,\ V=V_0)\)

- situatie 1: \(E<V_0\)
\[
\psi_{II}(x)=C e^{\kappa x}+D e^{-\kappa x},
\qquad
\kappa=\sqrt{\frac{2m(V_0-E)}{\hbar^2}}.
\]

- situatie 2: \(E>V_0\)
\[
\psi_{II}(x)=C e^{ik'x}+D e^{-ik'x},
\qquad
k'=\sqrt{\frac{2m(E-V_0)}{\hbar^2}}.
\]

**Region III** \((x>a,\ V=0)\) (scattering van links, enkel transmitted naar rechts)
\[
\psi_{III}(x)=F e^{ikx}.
\]

---

### Reflection & transmission

**Definities (via current/flux)**
\[
R=\frac{j_{\text{reflected}}}{j_{\text{incident}}},
\qquad
T=\frac{j_{\text{transmitted}}}{j_{\text{incident}}}.
\]

Voor bovenstaande golfvormen:
\[
R=\frac{|B|^2}{|A|^2},
\qquad
T=\frac{|F|^2}{|A|^2},
\qquad
R+T=1.
\]

---

### Resultaat case 1: \(E<V_0\) (tunneling)

\[
T=\left[1+\frac{(k^2+\kappa^2)^2}{4k^2\kappa^2}\sinh^2(\kappa a)\right]^{-1}
\]
Equivalent geschreven in \(E,V_0\):
\[
T=\left[1+\frac{V_0^2}{4E(V_0-E)}\sinh^2(\kappa a)\right]^{-1}
\]
\[
R=1-T.
\]

**Thick barrier approximation** (als \(\kappa a\gg 1\))
\[
T \approx \frac{16E(V_0-E)}{V_0^2}\,e^{-2\kappa a}.
\]

---

### Resultaat case 2: \(E>V_0\)

\[
T=\left[1+\frac{(k^2-k'^2)^2}{4k^2k'^2}\sin^2(k'a)\right]^{-1}
\]
Equivalent in \(E,V_0\):
\[
T=\left[1+\frac{V_0^2}{4E(E-V_0)}\sin^2(k'a)\right]^{-1}
\]
\[
R=1-T.
\]

**Resonantie**
\[
k'a=n\pi \ \Rightarrow\ \sin(k'a)=0 \Rightarrow T=1,\ R=0.
\]

**Classical vs quantum**
- klassiek (voor \(E>V_0\)): \(R=0,\ T=1\)
- kwantum: in het algemeen \(R\neq 0,\ T<1\) (behalve bij resonantie)

**Dimensionless parameters (zoals genoteerd)**
\[
\kappa^2 a^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(1-\frac{E}{V_0}\right)
\]
\[
k'^2 a^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(\frac{E}{V_0}-1\right)
\]

---

## 4.5 The infinite square well

\[
V(x)=
\begin{cases}
0, & 0<x<a\\
\infty, & x\le 0 \text{ of } x\ge a
\end{cases}
\]

Binnen de put (\(0<x<a\)):
\[
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}=E\psi
\quad\Rightarrow\quad
\psi(x)=A\cos(kx)+B\sin(kx),
\quad
k=\sqrt{\frac{2mE}{\hbar^2}}.
\]

Randvoorwaarden:
\[
\psi(0)=0 \Rightarrow A=0
\]
\[
\psi(a)=0 \Rightarrow \sin(ka)=0 \Rightarrow ka=n\pi,\ n=1,2,3,\dots
\]

Dus:
\[
k_n=\frac{n\pi}{a}
\]
\[
E_n=\frac{\hbar^2k_n^2}{2m}
=\frac{n^2\pi^2\hbar^2}{2ma^2}
\]

Genormaliseerde eigenfuncties:
\[
\psi_n(x)=\sqrt{\frac{2}{a}}\sin\left(\frac{n\pi x}{a}\right)
\quad (0<x<a).
\]

---

## Notes (energy levels, ground state, uncertainty estimate)

- Spectrum: oneindig veel **discrete** energieniveaus.
- quantum number \(n\) labelt de eigenstate.
- ground state: \(n=1\).

Uncertainty-principle estimate (orde van grootte):
\[
\Delta x \sim a
\Rightarrow
\Delta p \sim \frac{\hbar}{2a}
\Rightarrow
E\sim \frac{(\Delta p)^2}{2m}\sim \frac{\hbar^2}{8ma^2}.
\]

---

## Parity / symmetry van potential

Als
\[
V(x)=V(-x)
\]
dan is de Hamiltoniaan invariant onder \(x\mapsto -x\).

### Case 1: eigenvalue \(E\) is non-degenerate
Als \(\psi(x)\) een eigenfunctie is, dan is \(\psi(-x)\) ook een eigenfunctie met dezelfde \(E\).
Non-degenerate \(\Rightarrow\) ze moeten proportioneel zijn:
\[
\psi(-x)=c\,\psi(x)
\]
Nog eens \(x\mapsto -x\):
\[
\psi(x)=c\,\psi(-x)=c^2\psi(x)\Rightarrow c^2=1 \Rightarrow c=\pm 1.
\]
Dus:
\[
\psi(-x)=+\psi(x)\ (\text{even})\quad\text{of}\quad \psi(-x)=-\psi(x)\ (\text{odd}).
\]

### Case 2: eigenvalue \(E\) is degenerate
Dan bestaan er meerdere lineair onafhankelijke eigenfuncties met dezelfde \(E\).
Je kan lineaire combinaties bouwen die wél definite parity hebben, bv.
\[
\psi_\pm(x)=\psi(x)\pm\psi(-x).
\]

---

## Wave function normalization / time evolution in eigenbasis

Algemene oplossing TDSE (voor discrete basis):
\[
\Psi(x,t)=\sum_n c_n\,\psi_n(x)\,e^{-iE_n t/\hbar}
\]
met coëfficiënten uit beginconditie:
\[
c_n=\int \psi_n^*(x)\,\Psi(x,0)\,dx.
\]

Voor de infinite square well: \(E_n=n^2E_1\).  
Dan:
\[
\Psi(x,t+T)=\Psi(x,t)\quad\text{voor}\quad T=\frac{2\pi\hbar}{E_1}
\]
want
\[
e^{-iE_n T/\hbar}=e^{-in^2E_1(2\pi\hbar/E_1)/\hbar}=e^{-i2\pi n^2}=1.
\]

---

## 4.7 The quantum harmonic oscillator

Klassiek:
\[
E=\frac{p^2}{2m}+\frac{1}{2}m\omega^2x^2.
\]

Potentiaal near minimum (Taylor):
\[
V(x)\approx V(x_0)+\frac{1}{2}k(x-x_0)^2,
\quad
k=\left.\frac{d^2V}{dx^2}\right|_{x_0}.
\]

### Classical mechanics (LHO)
\[
m\ddot x + kx=0,\qquad \omega=\sqrt{\frac{k}{m}}
\]
\[
x(t)=x_0\cos(\omega t+\phi)
\]
\[
T=\frac{2\pi}{\omega}.
\]

### Quantum mechanics (LHO)
\[
\left[-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+\frac{1}{2}m\omega^2x^2\right]\psi(x)=E\psi(x).
\]

Introduce dimensionless variable:
\[
\xi=\sqrt{\frac{m\omega}{\hbar}}\,x,
\qquad
\varepsilon=\frac{2E}{\hbar\omega}.
\]

Dan wordt de vergelijking:
\[
\frac{d^2\psi}{d\xi^2}+(\varepsilon-\xi^2)\psi=0.
\]

---

## Asymptotics → Hermite equation

Voor grote \(|\xi|\):
- oplossingen gedragen zich als \(e^{-\xi^2/2}\) of \(e^{+\xi^2/2}\)
- normaliseerbaarheid \(\Rightarrow\) hou enkel \(e^{-\xi^2/2}\).

Ansatz:
\[
\psi(\xi)=e^{-\xi^2/2}H(\xi).
\]

Dan volgt:
\[
\frac{d^2H}{d\xi^2}-2\xi\frac{dH}{d\xi}+(\varepsilon-1)H=0
\]
(= Hermite equation).

Normaliseerbaarheid \(\Rightarrow H\) moet een **polynoom** zijn:
\[
\varepsilon=2n+1,\qquad n=0,1,2,\dots
\]

### Energy levels
\[
E_n=\left(n+\frac{1}{2}\right)\hbar\omega.
\]

### Hermite polynomials and wavefunctions

\[
\psi_n(\xi)=\mathcal N_n\,e^{-\xi^2/2}H_n(\xi).
\]

In \(x\)-variabele:

\[
\psi_n(x)=
\frac{1}{\sqrt{2^n n!}}
\left(\frac{m\omega}{\pi\hbar}\right)^{1/4}
\exp\left(-\frac{m\omega x^2}{2\hbar}\right)
H_n\!\left(\sqrt{\frac{m\omega}{\hbar}}\,x\right).
\]

---

## Hermite polynomials (eerste paar)

\[
H_0(\xi)=1
\]
\[
H_1(\xi)=2\xi
\]
\[
H_2(\xi)=4\xi^2-2
\]
\[
H_3(\xi)=8\xi^3-12\xi
\]
\[
H_4(\xi)=16\xi^4-48\xi^2+12
\]
\[
H_5(\xi)=32\xi^5-160\xi^3+120\xi
\]

**Orthogonality**
\[
\int_{-\infty}^{+\infty} e^{-\xi^2}H_m(\xi)H_n(\xi)\,d\xi
=\sqrt{\pi}\,2^n n!\,\delta_{mn}.
\]

---

## Bewijs van Gauss-integraal (zoals op het blad)

\[
I=\int_{-\infty}^{+\infty} e^{-x^2}\,dx
\]

\[
I^2=\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}
e^{-(x^2+y^2)}\,dx\,dy
\]

Ga naar poolcoördinaten:
\[
I^2=\int_0^{2\pi}\int_0^{+\infty} e^{-r^2}\,r\,dr\,d\theta=
2\pi\cdot\frac{1}{2}=\pi
\]
Dus:
\[
I=\sqrt{\pi}.
\]

---


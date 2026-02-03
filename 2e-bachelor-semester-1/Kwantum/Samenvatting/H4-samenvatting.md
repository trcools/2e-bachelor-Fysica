4.1 General formulae

TDSE:
\[
i\hbar \frac{\partial}{\partial t}\psi(x,t)=
\left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2}+V(x)\right]\psi(x,t)
\]

TISE:
\[
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}\psi(x)+V(x)\psi(x)=E\psi(x)
\]

\[
\vec{j}(x,t)=\frac{\hbar}{2mi}\left[\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x}-\psi(x,t)\frac{\partial\psi^*(x,t)}{\partial x}\right]
\]


4.2 The free particle

\[
E=\frac{\hbar^2k^2}{2m}
\]
\[
(V=0)\Rightarrow -\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}=E\psi(x)
\Rightarrow \frac{d^2\psi(x)}{dx^2}+\frac{2mE}{\hbar^2}\psi(x)=0
\Rightarrow \frac{d^2\psi(x)}{dx^2}+k^2\psi(x)=0
\]

dimensions: \[onleesbaar\] (iets met \(\left(\frac{2mE}{\hbar^2}\right)^{1/2}\))

general solution:
\[
\psi(x)=Ae^{ikx}+Be^{-ikx}
\]

interpretatie:

1) \(Ae^{ikx}\): vrij deeltje met impuls \(+\hbar k\) in de positieve \(x\)-richting

2) \(Be^{-ikx}\): vrij deeltje met impuls \(-\hbar k\) in de negatieve \(x\)-richting

\[
\psi(x,t)=Ae^{i(kx-\omega t)}+Be^{-i(kx+\omega t)}
\]

stationaire toestand: meting van de energie van het vrij deeltje zal steeds de waarde \(E=\hbar\omega\) opleveren.


Geval 1: \(B=0\) (vrij deeltje in de positieve \(x\)-richting)

- \(P(x,t)=P(x)=|A|^2\)

\[
\Delta x\to\infty \quad \text{want}\quad \Delta p_x=0
\]

- waarschijnlijkheidsdichtheid:
\[
j(x,t)=j(x)=\frac{\hbar k}{m}|A|^2=v|A|^2=vP(x)
\]

\[
j=v\cdot P \qquad \text{(flux = snelheid}\cdot\text{dichtheid)}
\]


Geval 2: \(A=0\) analoge resultaten

\[
P(x,t)=P(x)=|B|^2
\]

\[
j(x,t)=j(x)=-v|B|^2=-vP(x)
\]


Geval 3: \(A=B\) (twee vlakke golven die met gelijke amplitude in tegengestelde richting bewegen)

- \[
\psi(x,t)=2(A\cos(kx))\,e^{-i\omega t}
\]

\[
x_m\ \text{nodes}=\pm\frac{\left(\frac{\pi}{2}+m\pi\right)}{k}\ ,\quad m=0,1,2,\ldots
\]

- \[
P(x,t)=P(x)=|2A|^2\cos^2(kx)
\]

- \[
j(x,t)=j(x)=0
\]
(geen netto stroom, want \(\mathrm{im}(\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x})=0\))


Geval 4: \(A=-B\)

\[
P(x,t)=P(x)=|2A|^2\sin^2(kx)
\]

\[
j(x,t)=j_0=0
\]


Geval 5: \(A\) en \(B\) arbitrair (superpositie van 2 vlakke golven die in tegengestelde richting bewegen resp. met amplitude A en B)

\[
\psi(x,t)=Ae^{i(kx-\omega t)}+Be^{-i(kx+\omega t)}
\]

\[
P(x,t)=P(x)=|A|^2+|B|^2+\left(AB^*e^{2ikx}+A^*Be^{-2ikx}\right)
\]

met interferentie effecten tussen de 2 golven

\[
j(x,t)=j(x)=\frac{\hbar}{m}\ \mathrm{Im}\left(\psi^*(x,t)\frac{\partial\psi(x,t)}{\partial x}\right)
=v\left(|A|^2-|B|^2\right)
\]

netto flux van 2 golven: 1 in de positieve \(x\)-richting (amplitude A) en 1 in de negatieve \(x\)-richting (amplitude B)


* De golffunctie voor een vrij deeltje is niet kwadraat integreerbaar!!

In de praktijk wordt de beweging van het vrij deeltje vaak beperkt tot 1D met behulp van box normalisatie.

vrij deeltje in eindige doos:
met periodieke randvoorwaarden opgelegd: \(\psi_n(x=0)=\psi_n(x=L)\)

mogelijke impuls eigenwaarden en energie eigenwaarden:
\[
k_n=\frac{2\pi n}{L}\ ,\qquad
E_n=\frac{\hbar^2k_n^2}{2m}=\frac{2\pi^2\hbar^2}{mL^2}n^2
\quad (n=0,\pm1,\pm2,\ldots)
\]

de golffunctie kan genormeerd worden in de basis-doos \((0\le x\le L)\):
\[
\psi_n(x)=\frac{1}{\sqrt{L}}e^{ik_nx}
\qquad
\text{met}\quad
\int_0^L \psi_n^*(x)\psi_n(x)\,dx=1
\]


4.4 The potential barrier

(I)            (II)                 (III)
─────| 0        | a ───────────────────────> x
     |----------|
        V0

\[
V(x)=
\begin{cases}
0, & x<0 \\
V_0, & 0<x<a \\
0, & x>a
\end{cases}
\]

**Region I:** \(-\infty<x\le 0\)

\[
\psi(x)=Ae^{i\ell x}+Be^{-i\ell x}
\]
\[
E=\frac{\hbar^2\ell^2}{2m}
\]
\[
j_0=v(|A|^2-|B|^2)=j(\text{incoming})-j(\text{reflected})
\]

**Region II:** \(0\le x\le a\)

2 gevallen:
- **geval 1:** \(E<V_0\)
\[
\psi(x)=Fe^{Kx}+Ge^{-Kx}
\qquad \text{met}\quad
K=\left[\frac{2m(V_0-E)}{\hbar^2}\right]^{1/2}
\]

- **geval 2:** \(E>V_0\)
\[
\psi(x)=Fe^{i\ell' x}+Ge^{-i\ell' x}
\qquad \text{met}\quad
\ell'=\left[\frac{2m(E-V_0)}{\hbar^2}\right]^{1/2}
\]

**Region III:** \(a<x<+\infty\)

\[
\psi(x)=Ce^{i\ell x}
\quad
(\text{want deeltje alleen in de positieve }x\text{-richting beweegt})
\]
\[
E=\frac{\hbar^2\ell^2}{2m}
\]
\[
j_0=v|C|^2=j(\text{doorgelaten})
\]

* Reflectie-coëfficiënt \(R\):
\[
R=\frac{j(\text{reflected})}{j(\text{incoming})}
\]

* Transmissie-coëfficiënt \(T\):
\[
T=\frac{j(\text{doorgelaten})}{j(\text{incoming})}
\]

R en T worden berekend via de waarschijnlijkheidsstromen \(j(x,t)\)
(niet via de positie waarschijnlijkheid \(P(x,t)\)).


\[
R=\frac{|B|^2}{|A|^2}
\qquad ; \qquad
T=\frac{|C|^2}{|A|^2}
\]
\[
R+T=\frac{|B|^2+|C|^2}{|A|^2}=1
\]

Terug naar geval 1: \(E<V_0\)

De golffunctie op \(x=0\) en \(x=a\) is continu, dus:

met \(Ae^{i\ell x}+Be^{-i\ell x}=Fe^{Kx}+Ge^{-Kx}\) op \(x=0\) en \(x=a\) (lijnen), vinden we:
\[
\begin{cases}
A+B=F+G \\
Fe^{Ka}+Ge^{-Ka}=Ce^{i\ell a}
\end{cases}
\]

De \[onleesbaar\] afgeleiden op \(x=0,\ x=a\) moeten ook gelijk zijn om in dezelfde fase te blijven:
\[
\begin{cases}
i\ell(A-B)=K(F-G) \\
K(Fe^{Ka}-Ge^{-Ka})=i\ell Ce^{i\ell a}
\end{cases}
\]

We elimineren \(F\) en \(G\), en bekomen \(\frac{B}{A}\) en \(\frac{C}{A}\):

\[
\frac{B}{A}=
\frac{(\ell^2+K^2)(e^{2Ka}-1)}
{e^{2Ka}(\ell+iK)^2-(\ell-iK)^2}
\]

\[
\frac{C}{A}=
\frac{4i\ell K\,e^{-i\ell a}\,e^{Ka}}
{e^{2Ka}(\ell+iK)^2-(\ell-iK)^2}
\]

\[
T=\frac{|C|^2}{|A|^2}=
\left[1+\frac{(\ell^2+K^2)^2\sinh^2(Ka)}{4\ell^2K^2}\right]^{-1}=
\left[1+\frac{V_0^2\sinh^2(Ka)}{4E(V_0-E)}\right]^{-1}
\]

\[
R=\frac{|B|^2}{|A|^2}=
\left[1+\frac{4\ell^2K^2}{(\ell^2+K^2)^2\sinh^2(Ka)}\right]^{-1}=
\left[1+\frac{4E(V_0-E)}{V_0^2\sinh^2(Ka)}\right]^{-1}
\]

Nu kunnen we nog eens verifiëren dat \(R+T=1\) wat intuitief ook logisch is.


\[
\frac{mV_0a^2}{\hbar^2}
:\ \text{opaciteit van de potentiaal barriere}
\quad
(\text{maat van kwantum mechanische ondoordringbaarheid})
\]

\[
T \approx \frac{16E(V_0-E)}{V_0^2}\,e^{-2Ka}
\quad
(\text{toepassing in scanning tunneling microscope})
\]

Terug naar Case 2: \(E>V_0\)

\[
\psi(x)=Fe^{i\ell' x}+Ge^{-i\ell' x}
\qquad \text{met}\quad
\ell'=\left[\frac{2m(E-V_0)}{\hbar^2}\right]^{1/2}
\]

Analoge oplossing probleem vinden we:

\[
R=\frac{|B|^2}{|A|^2}=
\left[1+\frac{4\ell^2\ell'^2}{(\ell^2-\ell'^2)^2\sin^2(\ell'a)}\right]^{-1}=
\left[1+\frac{4E(E-V_0)}{V_0^2\sin^2(\ell'a)}\right]^{-1}
\]

\[
T=\frac{|C|^2}{|A|^2}=
\left[1+\frac{V_0^2\sin^2(\ell'a)}{4E(E-V_0)}\right]^{-1}
\]

gemakkelijk zien we dat \(R+T=1\)

klassiek: \(R=0\) en \(T=1\)

kwantum: \(R\ne 0\) en \(T<1\)

\[
(\ell a)^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(\frac{E}{V_0}\right)
\]
\[
(Ka)^2=\left(\frac{2mV_0a^2}{\hbar^2}\right)\left(1-\frac{E}{V_0}\right)
\]

T hangt enkel af van 2 dimensieloze grootheden:
de opaciteit \(\frac{mV_0a^2}{\hbar^2}\) en de ratio \(\frac{E}{V_0}\).


## 4.5 The infinite square well

\[
V(x)=
\begin{cases}
0, & -a<x<a \\
\infty, & |x|>a
\end{cases}
\]

\[
\psi(x=\pm a)=0
\]

Schrödinger for \(|x|<a\):
\[
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}=E\psi(x)
\Rightarrow \frac{d^2\psi(x)}{dx^2}+k^2\psi(x)=0
\qquad
\left(k=\left(\frac{2mE}{\hbar^2}\right)^{1/2}\right)
\]

general solution:
\[
\psi(x)=A\cos(kx)+B\sin(kx)
\]

if we apply the boundary conditions:
\[
A\cos(ka)=0 \quad \text{en} \quad B\sin(ka)=0
\]

1) \(B=0\Rightarrow \cos(ka)=0\)

The only allowed values of \(k\) are:
\[
k_m=\frac{m\pi}{2a}=\frac{m\pi}{L}
\qquad (m=1,3,5,\ldots)
\]

The corresponding eigenfunction
\[
\psi_m(x)=A_m\cos(k_m x)
\]
can be normalized so that:
\[
\int_{-a}^{a}\psi_m^*(x)\psi_m(x)\,dx=1
\]

The normalized eigenfunctions are:
\[
\psi_m(x)=\frac{1}{\sqrt{a}}\cos\!\left(\frac{m\pi}{2a}x\right)
\qquad (m=1,3,5,\ldots)
\]

2) \(A=0\Rightarrow \sin(ka)=0\)

The normalized eigenfunctions are:
\[
\psi_m(x)=\frac{1}{\sqrt{a}}\sin\!\left(\frac{m\pi}{2a}x\right)
\qquad (m=2,4,\ldots)
\]

\[
k_m=\frac{m\pi}{L}
\qquad
\lambda_m=\frac{2\pi}{k_m}=\frac{2L}{m}
\]

The energy is quantized, the energy eigenvalues being:
\[
E_m=\frac{\hbar^2 k_m^2}{2m}
=\frac{\hbar^2\pi^2 m^2}{8ma^2}
=\frac{\hbar^2\pi^2 m^2}{2mL^2}
\qquad (m=1,2,3,\ldots)
\]


## (p.8) vervolg: spectrum + orthogonaliteit + zero-point energy

The energy spectrum obtained consists of an infinite number of discrete energy levels corresponding to bound states.

eigenfunctions \(\psi_n(x)\) and \(\psi_m(x)\) corresponding to different eigenvalues \(E_n\) and \(E_m\) are orthogonal:
\[
\int_{-a}^{a}\psi_n^*(x)\psi_m(x)\,dx=0
\qquad n\ne m
\]

Note that the lowest energy (zero-point energy) is:
\[
E_1=\frac{\hbar^2\pi^2}{8ma^2}
\]
so there is no state of zero-energy!!

uncertainty principle:
\[
\Delta x\approx a \ \Rightarrow\ \Delta p_x \approx \frac{\hbar}{a}
\]
indicates an minimal kinetic energy:
\[
\frac{(\Delta p_x)^2}{2m}=\frac{\hbar^2}{2ma^2}=\frac{4E_1}{\pi^2}
\]

* Parity

The eigenfunctions of the cosine function are such that:
\[
\psi_n(-x)=\psi_n(x)
\]
and are therefore even function of \(x\),

while the eigenfunctions with a sine function in it are such that:
\[
\psi_n(-x)=-\psi_n(x)
\]
and hence are odd.

Reflection through the origin \(x\to -x\) is called the parity operation.

If the potential is symmetric the Hamiltonian
\[
H=-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V(x)
\]
does not change when \(x\) is replaced by \(-x\) (invariant under parity operation).


## (p.9) Case 1–2 (degeneracy)  [gedeeltelijk onleesbaar]

Case 1: The eigenvalue \(E\) is non-degenerate \([onleesbaar]\)

The eigenfunctions \(\psi_1(x)\) and \(\psi_2(x)\) corresponding to the same energy eigenvalue \(E\) can thus only differ by a multiplicative constant:
\[
\psi_1(x)=c\,\psi_2(x)
\]

Combining the above with the parity equation yields: \([onleesbaar]\)

and hence we may choose the eigenfunction \(\psi(x)\) have a definite parity being either even or odd.

Case 2: The eigenvalue \(E\) is degenerate \([onleesbaar]\)

In this case more than one linearly independent eigenfunction corresponds to the same energy eigenvalue. \([onleesbaar]\)

However, it is always possible to construct linear combinations of these eigenfunctions which do have definite parity. \([onleesbaar]\)

The Schrödinger equation becomes: \([onleesbaar]\)

\[
\left[-\frac{\hbar^2}{2m}\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}\right)+V(x,y,z)\right]\psi(x,y,z)=E\psi(x,y,z)
\]

Using parity operation \([onleesbaar]\)

We can have \(P\psi=\pm\psi\) \([onleesbaar]\)

The commutation equation becomes: \([onleesbaar]\)


* wave function representation

general solution TDSE:
\[
\psi(x,t)=\sum_n c_n\,\psi_n(x)\,e^{-iE_nt/\hbar}
\]

\[
c_n=\int \psi_n^*(x)\,\psi(x,t=0)\,dx
\]

\[
\psi(x,t=T)=\psi(x,t=0)
\]

Since \(E_n=nE_1\) we have
\[
\psi(x,t+T)=\sum_n c_n\psi_n(x)\exp\!\left(-i\frac{nE_1(t+T)}{\hbar}\right)
=\sum_n c_n\psi_n(x)e^{-inE_1t/\hbar}\exp\!\left(-i\frac{nE_1T}{\hbar}\right)
\]

Choose
\[
T=\frac{2\pi\hbar}{E_1}
\Rightarrow
\exp\!\left(-i\frac{nE_1T}{\hbar}\right)=\exp(-i2\pi n)=1
\]

and thus
\[
\psi(x,t+T)=\psi(x,t)
\]

----------------------------------------------------------------

4.7 The harmonic oscillator

Energy:
\[
E=\frac{1}{2}m\dot{x}^{\,2}+\frac{1}{2}kx^2
=\frac{1}{2}m\dot{x}^{\,2}+V(x)
\]

Potential:
\[
V(x)=\frac{1}{2}kx^2=\frac{1}{2}m\omega^2x^2
\qquad\text{met}\qquad
\omega=\sqrt{\frac{k}{m}}
\]

* LHO in de klassieke mechanica

\[
\frac{d^2x}{dt^2}=-\omega^2x
\Rightarrow x(t)=x_0\cos(\omega t)
\qquad
(\text{waar }x_0\text{ bepaald wordt door de startcondities})
\]

\[
v_{\max}=\omega x_0
\qquad
T=\frac{2\pi}{\omega}\ \ \text{(interval }(-x_0,x_0)\text{)}
\]

de “turning points” (waar \(v=0\)):
\[
\frac{1}{2}kx_0^2=E
\Rightarrow
x_0=\sqrt{\frac{2E}{k}}=\sqrt{\frac{2E}{m\omega^2}}
\]

* LHO in de quantum mechanica

TISE:
\[
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}+\frac{1}{2}m\omega^2x^2\,\psi(x)=E\psi(x)
\]

dimensionless eigenvalues:
\[
\lambda=\frac{2E}{\hbar\omega}
\qquad\text{where}\qquad
\omega=\left(\frac{k}{m}\right)^{1/2}
\]

we solicit the dimensionless variable \(\xi=\alpha x\) where
\[
\alpha=\left(\frac{m^2\omega^2}{\hbar^2}\right)^{1/4}
=\left(\frac{m\omega}{\hbar}\right)^{1/2}
\]

The Schrödinger equation then becomes
\[
\frac{d^2\psi(\xi)}{d\xi^2}+(\lambda-\xi^2)\psi(\xi)=0
\]

----------------------------------------------------------------

For large \(|\xi|\):
\[
\psi(\xi)=e^{-\xi^2/2}
\qquad
(\text{satisfy the equation as an asymptotic solution, } \xi\to\infty)
\]

For all values of \(\xi\):
\[
\psi(\xi)=e^{-\xi^2/2}H(\xi)
\]
(we obtain for \(H(\xi)\) the differential equation)

\[
\frac{d^2H}{d\xi^2}-2\xi\frac{dH}{d\xi}+(\lambda-1)H=0
\]
which is called the Hermite equation

* Energy levels

Energy spectrum of the LHO:
\[
E_n=\left(n+\frac{1}{2}\right)\hbar\omega
\qquad (n=0,1,2,\ldots)
\]

LHO with lowest state \((n=0)\):
\[
E_0=\frac{1}{2}\hbar\omega \qquad (\text{zero-point energy})
\]

* Hermite polynomials

\[
\psi_n(\xi)=\text{(const)}\ e^{-\xi^2/2}H_n(\xi)
\]

where the function \(H_n(\xi)\) is a polynomial of order \(n\).

Both \(\psi_n(\xi)\) and \(H_n(\xi)\) have parity \((-1)^n\).

The polynomials \(H_n(\xi)\) satisfy the Hermite equation with \(\lambda=2n+1\):
\[
\frac{d^2H_n}{d\xi^2}-2\xi\frac{dH_n}{d\xi}+2nH_n=0
\]

The Hermite polynomials:
\[
H_n(\xi)=(-1)^n e^{\xi^2}\frac{d^n}{d\xi^n}\left(e^{-\xi^2}\right)
\]

The first few Hermite polynomials \(H_n(\xi)\) are:
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

Orthogonaliteit:
\[
\int_{-\infty}^{\infty}e^{-\xi^2}H_m(\xi)H_n(\xi)\,d\xi=0
\qquad (m\ne n)
\]

Normering:
\[
\int_{-\infty}^{\infty}e^{-\xi^2}H_n^2(\xi)\,d\xi=\sqrt{\pi}\,2^n n!
\]

Energy-eigenfunctions in eigenwaarden van de lineaire harmonische oscillator:

Energy eigenfunctions:
\[
\psi_n(x)=\left(\frac{\alpha}{\sqrt{\pi}\,2^n n!}\right)^{1/2}
e^{-\alpha^2x^2/2}\,H_n(\alpha x)
\qquad
(\alpha=\sqrt{m\omega/\hbar})
\]
\[onleesbaar\ \text{alternatieve vorm met }x_0=\sqrt{\hbar/(m\omega)}\]

- omdat de Hamiltoniaan voor de LHO Hermitisch is geldt:
\[
\langle m|n\rangle=\int_{-\infty}^{\infty}\psi_m^*(x)\psi_n(x)\,dx=\delta_{mn}
\]

- verwachtingswaarde van \(x\):
\[
\langle x\rangle=\int_{-\infty}^{\infty}\psi_n^*(x)\,x\,\psi_n(x)\,dx=0
\]
(zoals \(\psi_n\) even of oneven is)

----------------------------------------------------------------

Bewijs van (4.178)

\[
I=\int_{-\infty}^{\infty}e^{-\xi^2}\,e^{(-s^2+2s\xi)}\,e^{(-t^2+2t\xi)}\,d\xi
\]

\[
=\int_{-\infty}^{\infty}e^{-\xi^2}\,e^{2\xi(s+t)}\,e^{-(s^2+t^2)}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\int_{-\infty}^{\infty}e^{-\left(\xi^2-2\xi(s+t)\right)}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\int_{-\infty}^{\infty}e^{-(\xi-(s+t))^2}\,e^{(s+t)^2}\,d\xi
\]

\[
=e^{-(s^2+t^2)}\,e^{(s+t)^2}\int_{-\infty}^{\infty}e^{-(\xi-(s+t))^2}\,d\xi
\]

zet \(x=\xi-(s+t)\):
\[
I=e^{2st}\int_{-\infty}^{\infty}e^{-x^2}\,dx
\]

Omdat
\[
\int_{-\infty}^{\infty}e^{-x^2}\,dx=\sqrt{\pi}
\]

zeggen we:
\[
I=\sqrt{\pi}\,e^{2st}
\]
\(\square\)

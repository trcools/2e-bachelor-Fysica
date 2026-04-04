# Problem 1: The Yukawa potential

We have the known potential:

$$
V(\mathbf r)=-g^2\frac{e^{-\alpha mr}}{r}
$$

## Electric field $\mathbf E(\mathbf r)$

We know that we can calculate the electric field from:

$$
\mathbf E=-\nabla V
$$

So:

$$
\begin{aligned}
\mathbf E
&= -\frac{\partial V(\mathbf r)}{\partial r}\,\hat{\mathbf r} \\
&= g^2\left[r^{-1}\frac{\partial (e^{-\alpha mr})}{\partial r} + e^{-\alpha mr}\frac{\partial (r^{-1})}{\partial r}\right]\hat{\mathbf r} \\
&= g^2\left[\frac{-\alpha m}{r}e^{-\alpha mr}-\frac{e^{-\alpha mr}}{r^2}\right]\hat{\mathbf r} \\
&= -g^2\frac{e^{-\alpha mr}}{r}\left[\alpha m + \frac{1}{r}\right]\hat{\mathbf r} \\
&= V(\mathbf r)\left[\alpha m + \frac{1}{r}\right]\hat{\mathbf r}
\end{aligned}
$$

In the last step we reused the given potential to keep the notation compact.

## Charge density $\rho(\mathbf r)$

To determine the charge density of our electric field, we use Gauss' law in differential form:

$$
\nabla\cdot\mathbf E(\mathbf r)=\frac{\rho(\mathbf r)}{\epsilon_0}
$$

where

$$
\mathbf E=-g^2\frac{e^{-\alpha mr}}{r}\left[\alpha m + \frac{1}{r}\right]\hat{\mathbf r}.
$$

In spherical coordinates:

$$
\nabla\cdot\mathbf E
= \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2E_r\right)
+ \frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\,E_\theta\right)
+ \frac{1}{r\sin\theta}\frac{\partial E_\phi}{\partial\phi}.
$$

Since our field only has an $r$-component:

$$
E_r=-g^2\frac{e^{-\alpha mr}}{r}\left[\alpha m+\frac{1}{r}\right],
$$

this reduces to

$$
\nabla\cdot\mathbf E = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2E_r\right)=\frac{\rho(\mathbf r)}{\epsilon_0}.
$$

So for $r>0$:

$$
\begin{aligned}
\rho(\mathbf r)
&= \frac{\epsilon_0}{r^2}\frac{\partial}{\partial r}\left(-g^2e^{-\alpha mr}(\alpha mr+1)\right) \\
&= \frac{-g^2\epsilon_0}{r^2}\left[\frac{\partial (e^{-\alpha mr})}{\partial r}(\alpha mr+1)+e^{-\alpha mr}\frac{\partial(\alpha mr+1)}{\partial r}\right] \\
&= \frac{-g^2\epsilon_0}{r^2}\left[-(\alpha m)e^{-\alpha mr}(\alpha mr+1)+(\alpha m)e^{-\alpha mr}\right] \\
&= \frac{g^2\epsilon_0}{r^2}\left[(\alpha m)^2e^{-\alpha mr}r\right] \\
&= \frac{(g\alpha m)^2\epsilon_0}{r}e^{-\alpha mr}
\end{aligned}
$$

This is the regular part ($r>0$). Because $V(\mathbf r)\sim -g^2/r$ for $r\to 0$, there is also a point-charge contribution at the origin. Therefore the full charge density is:

$$
\rho(\mathbf r)=\frac{(g\alpha m)^2\epsilon_0}{r}e^{-\alpha mr}-4\pi\epsilon_0 g^2\,\delta^{(3)}(\mathbf r).
$$

Short explanation of the delta term: $\delta^{(3)}(\mathbf r)$ is the 3D Dirac delta, i.e. a point-like distribution located at the origin. It is zero for $\mathbf r\neq 0$ and satisfies

$$
\int f(\mathbf r)\,\delta^{(3)}(\mathbf r)\,d^3r=f(\mathbf 0),
$$

so the term $-4\pi\epsilon_0 g^2\,\delta^{(3)}(\mathbf r)$ represents a point charge

$$
q_0=-4\pi\epsilon_0 g^2
$$

at $\mathbf r=0$.

## Total charge $Q$

$$
\int_0^\infty r e^{-\alpha mr}\,dr=\frac{1}{(\alpha m)^2}
\qquad (\alpha m>0).
$$

En voor de delta-term (in 3D, sferisch):

$$
\delta^{(3)}(\mathbf r)=\frac{\delta(r)}{4\pi r^2}
\Rightarrow
\int_0^\infty \delta^{(3)}(\mathbf r)\,r^2\,dr
=\frac{1}{4\pi}.
$$

Dus ook:

$$
4\pi\int_0^\infty \delta^{(3)}(\mathbf r)\,r^2\,dr=1.
$$

Daarom werkt je totale-ladingsterm netjes uit naar

$$
-16\pi^2\epsilon_0 g^2\left(\frac1{4\pi}\right)
=-4\pi\epsilon_0 g^2.
$$

To calculate the enclosed charge, we use Gauss' law in integral form on a spherical surface of radius $R$:

$$
\oint_S \mathbf E\cdot d\mathbf a = \frac{Q_{enclosed}(R)}{\epsilon_0}.
$$

With spherical symmetry:

$$
\oint_S\mathbf E\cdot d\mathbf a = 4\pi R^2E_r(R).
$$

Hence:

$$
Q_{enclosed}(R)=\epsilon_0\,4\pi R^2E_r(R)=-4\pi\epsilon_0 g^2 e^{-\alpha mR}(1+\alpha mR).
$$

So the total charge is:

$$
Q_{tot}=\lim_{R\to\infty}Q_{enclosed}(R)=0.
$$

Equivalent decomposition:

$$
Q_{regular}=\int\frac{(g\alpha m)^2\epsilon_0}{r}e^{-\alpha mr}\,d^3r=4\pi\epsilon_0 g^2,
$$

$$
Q_{delta}=\int-4\pi\epsilon_0 g^2\,\delta^{(3)}(\mathbf r)\,d^3r=-4\pi\epsilon_0 g^2,
$$

so indeed:

$$
Q_{tot}=Q_{regular}+Q_{delta}=0.
$$

# Problem 2

## Part a

We are given the vector potential

$$
\mathbf A(r,\theta,\phi)=
\begin{cases}
\dfrac{1}{3}\mu_0 R\omega\sigma\,r\sin\theta\,\hat{\mathbf\phi}, & r\le R,\\[4pt]
\dfrac{1}{3}\mu_0 R^4\omega\sigma\,\dfrac{1}{r^2}\sin\theta\,\hat{\mathbf\phi}, & r\ge R.
\end{cases}
$$

For the inside region ($r\le R$), only $A_\phi$ is nonzero:

$$
A_\phi=\frac{1}{3}\mu_0R\omega\sigma\,r\sin\theta.
$$

Using $\mathbf B=\nabla\times\mathbf A$ in spherical coordinates (with only $A_\phi$):

$$
B_r=\frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta A_\phi),
\qquad
B_\theta=-\frac{1}{r}\frac{\partial}{\partial r}(rA_\phi),
\qquad
B_\phi=0.
$$

Now compute:

$$
\begin{aligned}
B_r
&=\frac{1}{r\sin\theta}\frac{\partial}{\partial\theta}\left(\frac{1}{3}\mu_0R\omega\sigma\,r\sin^2\theta\right)
=\frac{2}{3}\mu_0R\omega\sigma\cos\theta,\\[4pt]
B_\theta
&=-\frac{1}{r}\frac{\partial}{\partial r}\left(\frac{1}{3}\mu_0R\omega\sigma\,r^2\sin\theta\right)
=-\frac{2}{3}\mu_0R\omega\sigma\sin\theta.
\end{aligned}
$$

$$
\hat{\mathbf r}
=\sin\theta\cos\phi,\hat{\mathbf x}
+\sin\theta\sin\phi,\hat{\mathbf y}
+\cos\theta,\hat{\mathbf z},
$$

$$
\hat{\boldsymbol\theta}
=\cos\theta\cos\phi,\hat{\mathbf x}
+\cos\theta\sin\phi,\hat{\mathbf y}
-\sin\theta,\hat{\mathbf z}.
$$

Als je deze twee combineert (of ($\hat{\mathbf z}$) isoleert), krijg je:

$$
\hat{\mathbf z}=\cos\theta,\hat{\mathbf r}-\sin\theta,\hat{\boldsymbol\theta}.
$$

Therefore

$$
\mathbf B_{\text{in}}
=\frac{2}{3}\mu_0R\omega\sigma\left(\cos\theta\,\hat{\mathbf r}-\sin\theta\,\hat{\boldsymbol\theta}\right).
$$

Using $\hat{\mathbf z}=\cos\theta\,\hat{\mathbf r}-\sin\theta\,\hat{\boldsymbol\theta}$, this is

$$
\boxed{\mathbf B_{\text{in}}=\frac{2}{3}\mu_0R\omega\sigma\,\hat{\mathbf z}}
$$

so the magnetic field inside the shell is uniform.

## Part b

Now consider a superconducting sphere in a uniform applied field

$$
\mathbf B_{\text{app}}=B\,\hat{\mathbf z}.
$$

For a superconductor (Meissner state), inside the sphere:

$$
\mathbf B_{\text{tot,in}}=0.
$$

By symmetry, the induced surface current has the form

$$
\mathbf K(\theta)=K_0\sin\theta\,\hat{\mathbf\phi}.
$$

From part a (identify $K_0=\sigma\omega R$), such a surface current produces inside:

$$
\mathbf B_{\text{ind,in}}=\frac{2}{3}\mu_0K_0\,\hat{\mathbf z}.
$$

Impose cancellation inside:

$$
\mathbf B_{\text{app}}+\mathbf B_{\text{ind,in}}=0
\quad\Rightarrow\quad
B+\frac{2}{3}\mu_0K_0=0,
$$

so

$$
K_0=-\frac{3B}{2\mu_0}.
$$

Hence the induced surface current distribution is

$$
\boxed{\mathbf K(\theta)=-\frac{3B}{2\mu_0}\sin\theta\,\hat{\mathbf\phi}}
$$

(same magnitude profile, sign depending on the chosen $\hat{\mathbf\phi}$ orientation).

Approximate field-line sketch description:

- No magnetic field lines inside the sphere ($\mathbf B=0$).
- Outside, field lines bend around the sphere.
- At the surface they are tangent (no normal component at the boundary).

# Problem 3

We have two infinite coaxial ideal conductors:

- inner conductor at radius $a$ with charge per unit length $+\lambda$ and current $\mathbf I_{\text{in}}=-I\,\hat{\mathbf z}$,
- outer conductor at radius $b>a$ with charge per unit length $-\lambda$ and return current $\mathbf I_{\text{out}}=+I\,\hat{\mathbf z}$.

## Part a: $\mathbf E$ and $\mathbf B$

By cylindrical symmetry, fields depend only on $r$.

### Electric field

Using Gauss:

$$
\oint \mathbf E\cdot d\mathbf a = \frac{Q_{\text{enc}}}{\epsilon_0}.
$$

So

$$
\mathbf E(r)=
\begin{cases}
0, & r<a,\\[4pt]
\dfrac{\lambda}{2\pi\epsilon_0 r}\,\hat{\mathbf r}, & a<r<b,\\[6pt]
0, & r>b.
\end{cases}
$$

(Outside, enclosed line charge is $\lambda+(-\lambda)=0$.)

### Magnetic field

Using Ampère:

$$
\oint \mathbf B\cdot d\mathbf l=\mu_0 I_{\text{enc}}.
$$

For $a<r<b$, enclosed current is only the inner one: $I_{\text{enc}}=-I$, hence

$$
B_\phi(2\pi r)=\mu_0(-I)
\quad\Rightarrow\quad
\mathbf B(r)=-\frac{\mu_0 I}{2\pi r}\,\hat{\boldsymbol\phi}.
$$

So

$$
\mathbf B(r)=
\begin{cases}
0, & r<a,\\[4pt]
-\dfrac{\mu_0 I}{2\pi r}\,\hat{\boldsymbol\phi}, & a<r<b,\\[8pt]
0, & r>b,
\end{cases}
$$

since for $r>b$ the enclosed current is $-I+I=0$.

## Part b: Poynting vector and its flux

$$
\mathbf S=\frac{1}{\mu_0}\mathbf E\times\mathbf B.
$$

In the region $a<r<b$:

$$
\mathbf S =\frac{1}{\mu_0}
\left(\frac{\lambda}{2\pi\epsilon_0 r}\hat{\mathbf r}\right) \times \left(-\frac{\mu_0 I}{2\pi r}\hat{\boldsymbol\phi}\right)
=-\frac{\lambda I}{4\pi^2\epsilon_0 r^2}\,\hat{\mathbf z}.
$$

So energy flows along $-\hat{\mathbf z}$ (for $\lambda>0$, $I>0$ with the given current directions).

Flux through a cross section (normal $+\hat{\mathbf z}$):

$$
\Phi_S=\iint \mathbf S\cdot d\mathbf a
=\int_a^b S_z\,(2\pi r\,dr)
=-\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

Hence the transmitted power magnitude is

$$
P=-\Phi_S=\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

$$
\mathbf S=\frac{1}{\mu_0}\mathbf E\times\mathbf B,\qquad
\mathbf E=\frac{\lambda}{2\pi\epsilon_0 r}\hat{\mathbf r},\qquad
\mathbf B=-\frac{\mu_0 I}{2\pi r}\hat{\boldsymbol\phi}
$$

$$
\Rightarrow\quad
\mathbf S=-\frac{\lambda I}{4\pi^2\epsilon_0 r^2}\hat{\mathbf z}.
$$

Neem nu een doorsnede (A) met normaal (+$\hat{\mathbf z}$):

$$
\Phi_S=\iint_A \mathbf S\cdot d\mathbf a
=\int_a^b S_z\,(2\pi r\,dr)
=-\frac{\lambda I}{2\pi\epsilon_0}\int_a^b \frac{dr}{r}
=-\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

Dus de **getransporteerde vermogensgrootte** is

$$
P=-\Phi_S=\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

Also, with

$$
V\equiv V(a)-V(b)=\int_a^b \frac{\lambda}{2\pi\epsilon_0 r}\,dr
=\frac{\lambda}{2\pi\epsilon_0}\ln\!\frac{b}{a},
$$

we get

$$
P=VI.
$$

## Part c: Resistor at $z=0$

After removing the $z<0$ half, the line is terminated at $z=0$ by a resistor $R$ between the two conductors.

### 1) Voltage across the resistor

At fixed $z=0$, the radial electric field in the dielectric region is still

$$
E_r(r)=\frac{\lambda}{2\pi\epsilon_0 r}, \qquad a<r<b.
$$

Hence the potential difference between inner and outer conductor is

$$
V\equiv V(a)-V(b)=\int_a^b E_r\,dr
=\frac{\lambda}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

This is exactly the voltage over the terminating resistor.

### 2) Power dissipated in the resistor

The same line current $I$ enters the termination, so

$$
P_R=VI=I^2R.
$$

### 3) Electromagnetic power arriving at $z=0$

From part b:

$$
S_z(r)=-\frac{\lambda I}{4\pi^2\epsilon_0 r^2}, \qquad a<r<b.
$$

Take the cross section $A$ at $z=0^+$ with normal $+\hat{\mathbf z}$:

$$
\Phi_S=\iint_A \mathbf S\cdot d\mathbf a
=\int_a^b S_z\,(2\pi r\,dr)
=-\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}.
$$

So the power flowing "+into" the resistor is

$$
P_{\text{in}}\equiv -\Phi_S
=\frac{\lambda I}{2\pi\epsilon_0}\ln\!\frac{b}{a}
=VI.
$$

Therefore

$$
\boxed{P_R=P_{\text{in}}=-\Phi_S}
$$

which proves that the resistor dissipates exactly the electromagnetic power transported by the fields.

# Problem 4: Geometrische schets

Onderstaande TikZ-figuur schetst de situatie van een lijnlading met ladingsdichtheid $\lambda$ die met constante snelheid $\mathbf v=v\hat{\mathbf z}$ langs zichzelf beweegt.

```latex



\begin{tikzpicture}[>=stealth,scale=1.1]
	% moving infinite line charge along z
	\draw[very thick] (0,-3) -- (0,3);
	\draw[->] (0,-3) -- (0,3.3) node[above] {$z$};
	\node[left] at (0,2.6) {$\lambda$};

	% velocity along line
	\draw[->,blue,thick] (0,1.4) -- (0,2.3) node[right] {$\mathbf v=v\hat{\mathbf z}$};

	% source element on line
	\coordinate (S) at (0,-1.2);
	\fill (S) circle (1.6pt);
	\node[left] at (S) {$dq=\lambda\,dz'$};

	% field point
	\coordinate (P) at (3.2,0.4);
	\fill (P) circle (1.6pt);
	\node[right] at (P) {$P(\mathbf r)$};

	% perpendicular distance rho from line to field point
	\coordinate (A) at (0,0.4);
	\draw[dashed] (A) -- (P);
	\node[above] at (1.6,0.4) {$\rho$};

	% separation vector R from source to field point
	\draw[->,red,thick] (S) -- (P) node[midway,below right] {$\mathbf R$};

	% angle theta between v (along +z) and R
	\draw (S) ++(90:0.55) arc[start angle=90,end angle=27,radius=0.55];
	\node at (0.55,-0.55) {$\theta$};
\end{tikzpicture}
```

Deze figuur maakt meteen duidelijk over welke variabele je integreert langs de lijn ($z'$), en hoe $\theta$ tussen $\mathbf v$ en $\mathbf R$ wordt gedefinieerd.

Met jouw hoekdefinitie geldt:

$$
\sin\theta=\frac{l}{R},\qquad \cos\theta=\frac{x}{R}
$$

dus

$$
x=R\cos\theta=\frac{l}{\sin\theta}\cos\theta=l\cot\theta.
$$

Daaruit volgt meteen:

$$
dx = -,l,\csc^2\theta,d\theta.
$$

Je zit heel dicht bij de juiste oplossing; de fout zit in de substitutiestap van \(I\).

$$
I=\int_0^\pi \frac{\sin\theta\,d\theta}{\left(1-\beta^2\sin^2\theta\right)^{3/2}},\qquad \beta\equiv \frac vc.
$$

Neem $(u=-\cos\theta$\), dan $du=\sin\theta\,d\theta$, en grenzen $θ\to\pi$ worden $-1\to 1$\.
De noemer wordt dan (hier zat je fout):

$$
1-\beta^2\sin^2\theta
=1-\beta^2(1-\cos^2\theta)
=1-\beta^2+\beta^2u^2.
$$

Dus:

$$
I=\int_{-1}^{1}\frac{du}{\left(1-\beta^2+\beta^2u^2\right)^{3/2}}.
$$

Factoriseer \(1-\beta^2\):

$$
I=\frac{1}{(1-\beta^2)^{3/2}}
\int_{-1}^{1}\frac{du}{\left(1+\frac{\beta^2}{1-\beta^2}u^2\right)^{3/2}}
=\frac{2}{1-\beta^2}.
$$

Daarmee:

$$
\mathbf E
=\frac{\lambda(1-\beta^2)}{4\pi\epsilon_0\,l}\,I\,\hat{\mathbf y}
=\frac{\lambda}{2\pi\epsilon_0\,l}\,\hat{\mathbf y}.
$$

En vervolgens:

$$
\mathbf B=\frac{1}{c^2}\mathbf v\times\mathbf E
\quad\Rightarrow\quad
|\mathbf B|=\frac{\mu_0 I}{2\pi l},\ \ I=\lambda v.
$$

Dus exact hetzelfde resultaat als elektrostatica/magnetostatica (voor gegeven \(\lambda\) en \(I\)).

# Problem 5: Moving magnetic monopole in a symmetric Maxwell theory

We use the generalized Maxwell equations (SI units):

$$
\nabla\cdot\mathbf E=\frac{\rho}{\epsilon_0},
\qquad
\nabla\cdot\mathbf B=\mu_0\rho_m,
$$

$$
\nabla\times\mathbf E+\frac{\partial\mathbf B}{\partial t}=-\mu_0\mathbf J_m,
\qquad
\nabla\times\mathbf B-\mu_0\epsilon_0\frac{\partial\mathbf E}{\partial t}=\mu_0\mathbf J.
$$

In het ruststelsel van de monopool:

$$
\mathbf B'(\mathbf r')=\frac{\mu_0 q_m}{4\pi}\frac{\mathbf r'}{r'^3},\qquad \mathbf E'=0.
$$

Voor constante snelheid $\mathbf v$ in het lab volgt (dual van Heaviside-veld):

$$
\mathbf B(\mathbf R,t)=
\frac{\mu_0 q_m}{4\pi}
\frac{1-\beta^2}{\left(1-\beta^2\sin^2\theta\right)^{3/2}}
\frac{\hat{\mathbf R}}{R^2}
\quad \beta=\frac{v}{c},
$$

met $\mathbf R=\mathbf r-\mathbf r_m(t)$ en $\theta=\angle(\mathbf R,\mathbf v)$, en

$$
\mathbf E(\mathbf R,t)= -\mathbf v\times \mathbf B(\mathbf R,t).
$$



Niet-relativistische limiet $(v\ll c)$:

$$
\mathbf B \approx \frac{\mu_0 q_m}{4\pi}\frac{\mathbf R}{R^3},\qquad
\mathbf E \approx -\frac{\mu_0 q_m}{4\pi}\frac{\mathbf v\times \mathbf R}{R^3}.
$$



For one monopole $q_m$ moving at constant velocity $\mathbf v$:

$$
\rho(\mathbf r,t)=0,
\qquad
\mathbf J(\mathbf r,t)=0,
$$

$$
\rho_m(\mathbf r,t)=q_m\,\delta^{(3)}\!\bigl(\mathbf r-\mathbf r_m(t)\bigr),
\qquad
\mathbf J_m(\mathbf r,t)=\rho_m\,\mathbf v,
$$

with

$$
\mathbf r_m(t)=\mathbf r_0+\mathbf v t.
$$

In the monopole rest frame:

$$
\mathbf E'(\mathbf r')=0,
\qquad
\mathbf B'(\mathbf r')=\frac{\mu_0 q_m}{4\pi}\frac{\hat{\mathbf r}'}{r'^2}.
$$

Transforming back to the lab frame (or by electric-magnetic duality with the moving electric charge field), define

$$
\mathbf R=\mathbf r-\mathbf r_m(t),
\qquad
R=|\mathbf R|,
\qquad
\beta=\frac{v}{c},
\qquad
\psi=\angle(\mathbf R,\mathbf v).
$$

Then

$$
\mathbf B(\mathbf r,t)=
\frac{\mu_0 q_m}{4\pi}
\frac{1-\beta^2}{\left(1-\beta^2\sin^2\psi\right)^{3/2}}
\frac{\hat{\mathbf R}}{R^2},
$$

$$
\mathbf E(\mathbf r,t)= -\,\mathbf v\times\mathbf B(\mathbf r,t).
$$

Equivalent explicit form:

$$
\mathbf E(\mathbf r,t)=
-\frac{\mu_0 q_m}{4\pi}
\frac{1-\beta^2}{\left(1-\beta^2\sin^2\psi\right)^{3/2}}
\frac{\mathbf v\times\hat{\mathbf R}}{R^2}.
$$

Low-velocity limit ($v\to 0$):

$$
\mathbf B\to\frac{\mu_0 q_m}{4\pi}\frac{\hat{\mathbf R}}{R^2},
\qquad
\mathbf E\to 0.
$$

Field-line sketch interpretation:

- $\mathbf B$-lijnen are radial from the monopole for $q_m>0$ (inward for $q_m<0$).
- $\mathbf E$-lijnen are closed loops around the velocity axis, with orientation from $-\mathbf v\times\mathbf B$.
- For large $v$, the field is strongest near $\psi\approx\pi/2$ and weaker near $\psi\approx 0,\pi$ (transverse compression).

# Extra opmerkingen cursus:

### P.466 (10.76):

Die (R)-term wordt niet “weggemoffeld”; hij wordt nul door het kruisproduct.

Uit

$$
\hat{\mathcal R}=\frac{\mathbf R}{\mathcal R}+\frac{\mathbf v}{c}
$$

volgt

$$
\mathbf B=\frac1c(\hat{\mathcal R}\times \mathbf E)
=\frac1c\left(\frac{\mathbf R}{\mathcal R}\times \mathbf E\right)+\frac1{c^2}(\mathbf v\times \mathbf E).
$$

Voor deze afleiding (lading met **constante** snelheid) geldt $\mathbf E \parallel \mathbf R) ((\mathbf E \propto \mathbf R)$. Daarom is

$$
\frac{\mathbf R}{\mathcal R}\times \mathbf E = \mathbf 0,
$$

en blijft alleen

$$
\mathbf B=\frac1{c^2}(\mathbf v\times \mathbf E)
$$

over.

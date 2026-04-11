# Stationary Analysis of a Battery-Buffered EV Charging Station
## A Piecewise-Deterministic Approach

*Andres Ferragut — Universidad ORT Uruguay*

---

## 1. A Battery-Buffered System Model

### 1.1 System Description

We consider a single EV charging station equipped with a battery buffer of capacity $B \in (0, \infty]$. The station has one charger connected to the electrical grid. Cars arrive and request service; a car that arrives when the charger is busy is turned away (loss system).

The key idea is that the battery mediates between the grid and the car charger. When the battery has energy, the car can be charged faster than the grid alone would allow. When the battery is depleted, the car receives energy directly from the grid at a lower rate.

### 1.2 Parameters

The system is characterized by the following rates:

- $\lambda > 0$: Poisson arrival rate of cars.
- $p_0$: grid-to-battery recharge rate when the charger is idle. Kept below $p_1$ to extend battery life.
- $p_1$: maximum grid supply rate.
- $p_2 > p_1$: fast charging rate when a car is present and the battery is positive.
- $c = p_2 - p_1 > 0$: net battery discharge rate during fast charging.
- $B$: battery capacity.
- $Y \sim \text{Exp}(\mu)$: energy demand of each car, exponentially distributed with mean $1/\mu$. We focus on this case throughout; general demands are discussed in Section 4.

The ordering is $p_0 \leq p_1 < p_2$, so $c > 0$.

### 1.3 Two-Phase Service

When a car is admitted, service proceeds in two phases:

- **Phase 1** (fast charging): while the battery is positive ($x > 0$), the car charges at rate $p_2$ and the battery depletes at rate $c$.
- **Phase 2** (slow charging): once the battery is exhausted ($x = 0$), the car charges at rate $p_1$ and the battery stays at zero.

A car may complete service entirely in Phase 1 if its demand is small relative to the available battery energy.

---

## 2. The PDMP Formulation

### 2.1 Piecewise-Deterministic Markov Processes

A piecewise-deterministic Markov process (PDMP) in the sense of Davis (1984) is a stochastic process $Z(t)$ evolving in $\mathbb{R}^d$ whose dynamics are characterized by three local objects:

- A **vector field** $r(z)$, which drives deterministic motion between jumps via $\dot{z} = r(z)$.
- A **jump rate** $\lambda(z) \geq 0$, giving the intensity of spontaneous jumps when the system is at state $z$.
- A **transition kernel** $Q(z;\,\cdot)$, a probability measure specifying the post-jump location given pre-jump state $z$.

The process follows the deterministic flow until a jump occurs, then relocates according to $Q$, and restarts.

### 2.2 State Space

The natural state of the system is $Z(t) = (n(t), x(t))$, where:

- $n(t) \in \{0, 1\}$: charger status ($0$ = idle, $1$ = busy).
- $x(t) \in [0, B]$: battery level.

The state space is $E = \{0, 1\} \times [0, B]$.

### 2.3 The State Space and the PDMP Framework

Our model is naturally formulated as a PDMP, but it sits outside the standard framework of Davis (1984) for two reasons. First, Davis requires the state space to be an **open** subset of $\mathbb{R}^d$. Second, Davis requires that the integrated jump rate satisfies $\Lambda(t^*(x), x) = \infty$ whenever $t^*(x) < \infty$ — meaning the process must jump spontaneously before the flow reaches the boundary. In our model, $\lambda$ is finite and constant, so $\Lambda(t^*(x), x) = \lambda \cdot t^*(x) < \infty$: the flow can reach $x=0$ or $x=B$ in finite time without a jump. Davis also explicitly prohibits self-jumps ($Q(\{x\};x) = 0$), which rules out forcing the process to jump back to the same boundary point.

Our model is closely related to **Markov-modulated fluid queues** (Rogers 1994, Kulkarni 1997), where a continuous level $x$ is driven by a finite-state Markov chain $n$, confined to $[0,B]$. The forward equations for such models are well-established. However, our model goes beyond the standard fluid queue setting: the jump rate of $n$ depends on $x$ — specifically, the departure rate switches between $\mu p_2$ and $\mu p_1$ at $x=0$. In a standard fluid queue the modulating chain is autonomous; here it is not.

We therefore proceed as follows. The process $Z(t) = (n(t), x(t))$ is constructed explicitly via its sample paths: starting from $(n_0, x_0)$, the trajectory follows the deterministic flow $\dot{x} = r(n,x)$ until the first event of a non-homogeneous Poisson process with rate $\lambda(n, x(t))$, at which point a jump occurs and the construction restarts. The flow is well-posed by Section 2.5. The inter-jump time has survival function:

$$P(\tau > t) = \exp\left(-\int_0^t \lambda(n, \phi(s, x_0))\,ds\right)$$

where $\phi(s, x_0)$ is the deterministic flow. Since the flow is deterministic given the current state, the distribution of $\tau$ depends only on $(n_0, x_0)$ and not on the history. This gives the **strong Markov property**: at any stopping time $\tau$, the future evolution depends only on $Z(\tau)$.

The resulting process is a **right process** in the sense of Sharpe (1988) — it has right-continuous paths, satisfies the strong Markov property, and admits a well-defined resolvent. This gives us access to the extended generator via the martingale problem (Davis 1993, Section 26): $\psi \in \mathcal{D}(\mathcal{A})$ if there exists a measurable $\mathcal{A}\psi$ such that

$$\psi(Z(t)) - \psi(Z(0)) - \int_0^t \mathcal{A}\psi(Z(s))\,ds$$

is a local martingale. The domain consists of functions $\psi$ that are absolutely continuous in $x$ for each fixed $n \in \{0,1\}$. No boundary condition is needed in the domain because $r$ vanishes at $x=0$ and $x=B$, killing the transport term there. The explicit form of $\mathcal{A}$ is given in Section 2.7, and we verify the martingale property directly for the test functions used in the analysis.

### 2.4 Vector Field

Between jumps, only $x$ evolves. The vector field is:

| $n$ | condition | $\dot{x} = r(n,x)$ |
|-----|-----------|----------------------|
| $0$ | $x < B$   | $+p_0$               |
| $0$ | $x = B$   | $0$                  |
| $1$ | $x > 0$   | $-c$                 |
| $1$ | $x = 0$   | $0$                  |

Setting $r = 0$ at the boundaries $x = 0$ (when $n=1$) and $x = B$ (when $n=0$) means the flow stops at the boundary rather than exiting $E$. The process dwells at the boundary until the next jump.

### 2.5 Well-Posedness of the Flow

The vector field $r(n, \cdot)$ is not Lipschitz at the boundary points: it drops from $p_0$ to $0$ at $x = B$ (when $n=0$) and from $-c$ to $0$ at $x = 0$ (when $n=1$). We nevertheless have a unique flow.

**Theorem (Uniqueness via Monotonicity).** *If $f : \mathbb{R} \to \mathbb{R}$ is non-increasing and bounded, then the ODE $\dot{x} = f(x)$ has at most one Carathéodory solution for each initial condition.*

*Proof.* Let $x(t)$ and $y(t)$ be two Carathéodory solutions with $x(0) = y(0) = x_0$. Since $f$ is bounded and measurable, both solutions are absolutely continuous, hence differentiable a.e. Set $d(t) = x(t) - y(t)$, so $d(0) = 0$ and $d$ is absolutely continuous. Writing the integral form:

$$d(t)^2 = d(0)^2 + \int_0^t \frac{d}{ds}d(s)^2\,ds = \int_0^t 2\,d(s)\,[f(x(s)) - f(y(s))]\,ds$$

Since $f$ is non-increasing, $d(s) > 0$ implies $x(s) > y(s)$ implies $f(x(s)) \leq f(y(s))$, so $d(s)[f(x(s))-f(y(s))] \leq 0$. Similarly for $d(s) < 0$. Hence the integrand is $\leq 0$ a.e., giving $d(t)^2 \leq 0$ for all $t \geq 0$, so $d(t) = 0$. $\square$

In our model, $r(0, \cdot)$ is non-increasing (drops from $p_0$ to $0$ at $x = B$) and $r(1, \cdot)$ is non-increasing (drops from $-c$ to $0$ at $x = 0$). Existence of solutions follows from the Carathéodory theorem (the vector field is bounded and measurable). Therefore the flow is uniquely defined.

### 2.6 Jumps

The only jumps are spontaneous:

- **Arrival** at $(0, x)$: rate $\lambda$. Kernel: $(0, x) \to (1, x)$. The battery level is unchanged; the charger becomes busy.
- **Departure** at $(1, x)$ with $x > 0$: rate $\mu p_2$. Kernel: $(1, x) \to (0, x)$. The charger becomes idle and the battery begins recharging.
- **Departure** at $(1, 0)$: rate $\mu p_1$. Kernel: $(1, 0) \to (0, 0)$.

**Remark.** An arriving car always finds $x > 0$. Since the battery charges at rate $p_0 > 0$ during idle periods, and the time until the next arrival is $\text{Exp}(\lambda)$ (hence positive a.s.), the battery is strictly positive at every arrival epoch. The state $(n=0, x=0)$ has zero stationary probability.

### 2.7 Generator

We derive the extended generator $\mathcal{A}$ via the martingale problem. Write $\psi_n(x) = \psi(n,x)$ for the restriction of $\psi$ to mode $n \in \{0,1\}$.

**Domain.** $\psi \in \mathcal{D}(\mathcal{A})$ if $\psi_n(\cdot)$ is absolutely continuous on $[0,B]$ for each $n \in \{0,1\}$, with derivative $\psi_n'$ defined Lebesgue-a.e. No boundary conditions are required since $r$ vanishes at $x=0$ and $x=B$.

**Derivation.** For $\psi \in \mathcal{D}(\mathcal{A})$, write the sample path decomposition:

$$\psi(Z(t)) - \psi(Z(0)) = \underbrace{\int_0^t r(n(s), x(s))\,\psi_{n(s)}'(x(s))\,ds}_{\text{flow contribution}} + \underbrace{\sum_{T_k \leq t} [\psi(Z(T_k)) - \psi(Z(T_k^-))]}_{\text{jump contribution}}$$

The flow term is well-defined a.s. since $r$ is bounded and $\psi_n' \in L^1[0,B]$ by absolute continuity.

For the jump term: in our model the jump kernel is a **point mass** — each jump changes only $n$ while leaving $x$ unchanged. So each jump contributes $\psi(1-n(T_k^-), x(T_k^-)) - \psi(n(T_k^-), x(T_k^-))$, which is a.s. bounded since $\psi_n$ is locally bounded (being absolutely continuous). The compensated jump sum:

$$\sum_{T_k \leq t}[\psi(Z(T_k)) - \psi(Z(T_k^-))] - \int_0^t \lambda(n(s),x(s))\,[\psi(1-n(s), x(s)) - \psi(n(s), x(s))]\,ds$$

is a local martingale by the standard theory of point processes (Davis 1984, Section 4), since the jumps are a.s. bounded.

Combining:

$$M_t^\psi := \psi(Z(t)) - \psi(Z(0)) - \int_0^t \mathcal{A}\psi(Z(s))\,ds$$

is a local martingale, where:

$$\mathcal{A}\psi(0, x) = r(0, x)\,\psi_0'(x) + \lambda\,[\psi(1, x) - \psi(0, x)]$$

$$\mathcal{A}\psi(1, x) = r(1, x)\,\psi_1'(x) + \mu p(x)\,[\psi(0, x) - \psi(1, x)]$$

with $p(x) = p_2\,\mathbf{1}_{x>0} + p_1\,\mathbf{1}_{x=0}$.

**Remark ($B = \infty$).** When $B = \infty$ the domain condition and derivation are unchanged: $\psi_n$ absolutely continuous on $[0,\infty)$ ensures $\psi_n' \in L^1_{\text{loc}}[0,\infty)$, making the flow term well-defined a.s. along any trajectory. The jump term remains a.s. bounded since the kernel is still a point mass leaving $x$ unchanged. Hence $M_t^\psi$ is a local martingale for all $\psi \in \mathcal{D}(\mathcal{A})$, with no growth conditions on $\psi$ required.

### 2.8 Kolmogorov Forward Equation

**General formulation.** Let $\pi_t$ denote the distribution of $Z(t)$, viewed as a general measure on $E = \{0,1\} \times [0,B]$. For any $\psi \in \mathcal{D}(\mathcal{A})$, the local martingale $M_t^\psi$ gives (assuming sufficient integrability to take expectations):

$$\frac{d}{dt}\int \psi\,d\pi_t = \int \mathcal{A}\psi\,d\pi_t$$

A stationary measure $\pi$ satisfies $\int \mathcal{A}\psi\,d\pi = 0$ for all $\psi \in \mathcal{D}(\mathcal{A})$. This is the weak form of $\mathcal{A}^*\pi = 0$.

**Structural ansatz.** We look for stationary measures of the following mixed form:

$$\pi = f_0(x)\,dx\big|_{(0,B)} + f_1(x)\,dx\big|_{(0,B)} + \pi_{0,B}\,\delta_{(0,B)} + \pi_{1,0}\,\delta_{(1,0)}$$

where $f_0, f_1 \geq 0$ are densities on $(0,B)$ and $\pi_{0,B}, \pi_{1,0} \geq 0$ are atoms at the two boundary points where the flow accumulates. Note that $(0,0)$ and $(1,B)$ do not support atoms since the flow points away from these points — any mass placed there is immediately transported into the interior.

In principle four atoms could appear at the four corners of $E$, but atoms at $(0,0)$ and $(1,B)$ are ruled out by the flow direction: $r(0,0^+) = p_0 > 0$ pushes mass away from $(0,0)$, and $r(1,B^-) = -c < 0$ pushes mass away from $(1,B)$. The existence of a stationary measure with this structure will be confirmed a posteriori by the explicit solution in Section 3.

**Deriving the balance equations.** Substituting $\pi$ into $\int \mathcal{A}\psi\,d\pi = 0$:

$$\int_0^B \left[r(0,x)\psi_0'(x) + \lambda(\psi_1-\psi_0)\right]f_0\,dx + \int_0^B\left[r(1,x)\psi_1'(x) + \mu p_2(\psi_0-\psi_1)\right]f_1\,dx$$
$$+ \lambda(\psi_1(B)-\psi_0(B))\pi_{0,B} + \mu p_1(\psi_0(0)-\psi_1(0))\pi_{1,0} = 0$$

Integrating the transport terms by parts. For the $n=0$ term, $r(0,0^+) = p_0 > 0$ (flow points away from $x=0$, no outward flux there) and $r(0,B^-) = p_0 > 0$ (flow arrives at $x=B$):

$$\int_0^B r(0,x)\psi_0'(x)f_0\,dx = p_0\psi_0(B^-)f_0(B^-) - \int_0^B\psi_0\,\partial_x[r(0,x)f_0]\,dx$$

For the $n=1$ term, $r(1,0^-) = -c < 0$ (flow arrives at $x=0$) and $r(1,B^-) = -c < 0$ (flow points away from $x=B$, no outward flux there):

$$\int_0^B r(1,x)\psi_1'(x)f_1\,dx = c\,\psi_1(0^+)f_1(0^+) - \int_0^B \psi_1\,\partial_x[r(1,x)f_1]\,dx$$

Requiring the result to hold for all $\psi \in \mathcal{D}(\mathcal{A})$ and identifying coefficients of $\psi_0(x)$, $\psi_1(x)$ on $(0,B)$, and $\psi_0(B)$, $\psi_1(0)$ separately:

**Densities** on $(0,B)$:

$$0 = -\frac{\partial}{\partial x}[r(0,x)f_0] - \lambda f_0 + \mu p_2 f_1 \tag{E0}$$

$$0 = -\frac{\partial}{\partial x}[r(1,x)f_1] + \lambda f_0 - \mu p_2 f_1 \tag{E1}$$

**Atom at $(0,B)$** (coefficient of $\psi_0(B)$):

$$0 = p_0\, f_0(B^-) - \lambda\, \pi_{0,B} \tag{EB}$$

**Atom at $(1,0)$** (coefficient of $\psi_1(0)$):

$$0 = c\, f_1(0^+) - \mu p_1\, \pi_{1,0} \tag{E0'}$$

together with the normalization $\int_0^B(f_0+f_1)\,dx + \pi_{0,B} + \pi_{1,0} = 1$.

---

## 3. Computing the Stationary Distribution

### 3.1 Irreducibility and Ergodicity

We argue that the process $Z(t)$ is **irreducible**: for every starting point $z \in E$ and every open set $U \subset E$, there exists $t > 0$ such that $P_z(Z(t) \in U) > 0$. Since any open set in $E$ contains an open rectangle $\{n_0\} \times (a,b)$, it suffices to show these are reachable from any starting point.

From any $(0, x_0)$: the battery drifts at rate $p_0 > 0$, so it reaches any level in $(x_0, B]$ in deterministic finite time. The probability of no arrival during this time is $e^{-\lambda t} > 0$, so $\{0\} \times (a,b)$ is reachable with positive probability for any $(a,b) \subset (0,B)$.

From any $(1, x_0)$: a departure occurs at rate at least $\mu p_1 > 0$, so mode $n=0$ is entered in finite time a.s. Then apply the argument above.

For $n_0 = 1$: an arrival occurs at rate $\lambda > 0$ from any idle state, and the battery drifts at rate $-c$ while busy, reaching any level in $(0, x_0)$ before the next departure with positive probability.

Since $E$ is compact and $Z(t)$ is irreducible, the process is **ergodic**: there exists a unique stationary distribution $\pi$, and for any initial condition $z \in E$ and any bounded measurable $f$:

$$\frac{1}{t}\int_0^t f(Z(s))\,ds \xrightarrow{t\to\infty} \int_E f\,d\pi \qquad \text{a.s.}$$

This follows from the general theory of Harris-recurrent Markov processes (Meyn and Tweedie, 1993, Chapter 13). The performance metrics in Section 3.3 are interpreted as expectations under $\pi$ via this ergodic theorem.

We characterize $\pi$ explicitly in the next section by solving the Kolmogorov forward equations $\mathcal{A}^*\pi = 0$. The uniqueness of the solution confirms the uniqueness of $\pi$.

### 3.2 Stationarity Conditions and Solution

Setting $\partial/\partial t = 0$ in the Kolmogorov forward equations of Section 2.8 gives the following system.

**Interior equations** for $x \in (0,B)$:
$$p_0\, f_0'(x) = -\lambda\, f_0(x) + \mu p_2\, f_1(x) \tag{E0}$$
$$-c\, f_1'(x) = \lambda\, f_0(x) - \mu p_2\, f_1(x) \tag{E1}$$

**Atom equations** from (EB) and (E0') of Section 2.8:
$$p_0\, f_0(B^-) = \lambda\, \pi_{0,B} \tag{B2}$$
$$c\, f_1(0^+) = \mu p_1\, \pi_{1,0} \tag{B1}$$

**Flux balance at $(0,0)$.** The point $(0,0)$ carries no atom (the flow $\dot{x}=p_0>0$ pushes mass away immediately). However, mass arrives at $(0,0)$ via departures from the atom $\pi_{1,0}$ at rate $\mu p_1\pi_{1,0}$, and is immediately carried into the interior at rate $p_0 f_0(0^+)$. Flux conservation gives:
$$p_0\, f_0(0^+) = \mu p_1\, \pi_{1,0} \tag{B0}$$

Combining (B0) and (B1): $p_0\, f_0(0^+) = c\, f_1(0^+)$.

Similarly, flux balance at $(1,B)$ (no atom, flow $\dot{x}=-c<0$ pushes mass away, mass arrives via arrivals from $\pi_{0,B}$ at rate $\lambda\pi_{0,B}$):
$$c\, f_1(B^-) = \lambda\, \pi_{0,B} \tag{B3}$$

Combining (B2) and (B3): $p_0\, f_0(B^-) = c\, f_1(B^-)$.

**Step 1: Proportionality.** Adding (E0) and (E1):
$$\frac{d}{dx}\bigl[p_0\, f_0(x) - c\, f_1(x)\bigr] = 0$$
so $p_0 f_0 - c f_1 = \text{const}$. Integrating (E0) and (E1) from $0$ to $x$ and adding gives the constant equals $p_0 f_0(0^+) - c f_1(0^+) = 0$ by the flux balance above. Hence:
$$\boxed{f_1(x) = \frac{p_0}{c}\, f_0(x) \qquad \text{for all } x \in (0, B)}$$

This is the **decoupling property**: in steady state, conditioned on $x \in (0,B)$, the conditional distribution of $n$ given $x$ does not depend on $x$:
$$\Pr(n=1 \mid x) = \frac{p_0}{p_0+c} \qquad \text{for all } x \in (0,B).$$

**Step 2: ODE for $f_0$.** Substituting into (E0):
$$f_0'(x) = \theta\, f_0(x), \qquad \theta = \frac{\mu p_2}{c} - \frac{\lambda}{p_0}$$

The unique solution is $f_0(x) = f_0(0^+)\, e^{\theta x}$.

**Physical meaning of $\theta$.** The quantity $\theta$ is the mean drift of the battery level in the absence of boundary constraints: positive when the station is lightly loaded (battery tends to fill), negative when heavily loaded (battery tends to drain), zero at critical balance. For the single-server model, $\theta = \mu p_2/c - \lambda/p_0$; for the multi-server extension see Section 4.

**Step 3: Point masses.** From (B1) and (B2):
$$\pi_{1,0} = \frac{p_0\, f_0(0^+)}{\mu p_1}, \qquad \pi_{0,B} = \frac{p_0\, f_0(0^+)}{\lambda}\, e^{\theta B}$$

**Step 4: Normalization.** Setting $f_0(0^+) = 1/K$ where:
$$K = \frac{p_0 + c}{c} \cdot \frac{e^{\theta B} - 1}{\theta} + \frac{p_0}{\mu p_1} + \frac{p_0}{\lambda}\, e^{\theta B}$$

(For $\theta = 0$ use $(e^{\theta B}-1)/\theta = B$.)

**Complete stationary distribution:**
$$f_0(x) = \frac{e^{\theta x}}{K}, \qquad f_1(x) = \frac{p_0}{c} \cdot \frac{e^{\theta x}}{K}, \qquad x \in (0, B)$$
$$\pi_{1,0} = \frac{p_0}{\mu p_1 K}, \qquad \pi_{0,B} = \frac{p_0}{\lambda K}\, e^{\theta B}$$



### 3.3 Performance Metrics

**Blocking probability.** By PASTA (Poisson Arrivals See Time Averages), an arriving car sees the stationary distribution. It is blocked iff $n = 1$:

$$P_B = \int_0^B f_1(x)\,dx + \pi_{1,0} = \frac{p_0}{K}\left(\frac{e^{\theta B}-1}{c\,\theta} + \frac{1}{\mu p_1}\right)$$

**Limit $B\to 0$:** $P_B \to \dfrac{\lambda}{\lambda + \mu p_1}$, corresponding to an M/M/1/1 system served always at the slow rate $p_1$.

**Limit $B\to\infty$, $\theta > 0$ (lightly loaded):** $K \sim e^{\theta B}\left(\frac{p_0+c}{c\theta}+\frac{p_0}{\lambda}\right)\to\infty$, so $\pi_{1,0}\to 0$ and $P_B \to \dfrac{\lambda}{\lambda+\mu p_2}$, corresponding to an M/M/1/1 system served always at the fast rate $p_2$.

**Limit $B\to\infty$, $\theta < 0$ (heavily loaded):** $K$ converges to a finite limit $K_\infty = \frac{p_0+c}{c|\theta|}+\frac{p_0}{\mu p_1}$, and $P_B$ converges to a value strictly between $\frac{\lambda}{\lambda+\mu p_2}$ and $\frac{\lambda}{\lambda+\mu p_1}$. Even with an infinite battery, the stationary distribution retains a positive mass at $x=0$ (slow charging) because the battery tends to drain in steady state.

**Mode fractions.** The fraction of time spent in each mode:

$$\rho_{\text{fast}} = \int_0^B f_1(x)\,dx = \frac{p_0}{cK} \cdot \frac{e^{\theta B}-1}{\theta}$$

$$\rho_{\text{slow}} = \pi_{1,0} = \frac{p_0}{\mu p_1 K}$$

$$\rho_{\text{idle}} = 1 - P_B.$$

Note $\rho_{\text{slow}} \to 0$ as $B\to\infty$ when $\theta > 0$: with a large battery in the lightly loaded regime, slow charging disappears. When $\theta < 0$, $\rho_{\text{slow}}$ remains positive even as $B\to\infty$.

**Mean battery level.**

$$\mathbb{E}[x] = \int_0^B x\,[f_0(x) + f_1(x)]\,dx + B\,\pi_{0,B}$$
$$= \frac{p_0 + c}{cK}\left[\frac{B\,e^{\theta B}}{\theta} - \frac{e^{\theta B}-1}{\theta^2}\right] + \frac{p_0 B}{\lambda K}\,e^{\theta B}.$$

**Mean charging time.** By Little's law applied to the loss system: the mean number in service equals $P_B$ (at most one car, busy with probability $P_B$), and the throughput of admitted cars is $\lambda(1-P_B)$. Therefore:

$$\mathbb{E}[T] = \frac{P_B}{\lambda(1-P_B)}$$

with $P_B$ given explicitly above. No integration is needed.

**Normalized mean charging time.** The minimum possible mean charging time — if the battery were always available — would be $\mathbb{E}[Y]/p_2 = 1/(\mu p_2)$. The ratio:

$$\frac{p_2\,\mathbb{E}[T]}{\mathbb{E}[Y]} = \frac{\mu p_2\,P_B}{\lambda(1-P_B)}$$

measures the average slowdown weighted by sojourn time. It equals 1 when all service is at rate $p_2$ (large battery, light load) and approaches $p_2/p_1$ when $B\to 0$ (no battery, all service at rate $p_1$).

Figure 2 shows this metric as a function of $B$ for the same parameters as Figure 1. All curves decrease from $p_2/p_1=3$ at $B=0$ toward their respective limits as $B\to\infty$.

![Normalized mean charging time vs battery capacity](normalized_time.png)

**Conditional charging time and slowdown.** Given admission at battery level $x$ with demand $Y = y$, the total charging time is:

$$T(x, y) = \frac{y}{p_2} + \left(\frac{1}{p_1} - \frac{1}{p_2}\right)\max\!\left(y - \frac{p_2 x}{c},\, 0\right)$$

Define the **critical battery level** $\tilde{x}(y) = cy/p_2$ — the battery level at which the car would exactly deplete the battery during service. If $x \geq \tilde{x}(y)$ the battery is never exhausted and the car is served entirely at rate $p_2$. If $x < \tilde{x}(y)$ the battery runs out and the car finishes at the slower rate $p_1$.

**Mean charging time given $Y = y$.** The admitted battery level has density $g(x) = f_0(x)/(1-P_B)$ on $(0,B)$ plus atom $\pi_{0,B}/(1-P_B)$ at $x=B$. Averaging over the admitted battery distribution:

$$\mathbb{E}[T \mid Y=y] = \frac{y}{p_2} + \left(\frac{1}{p_1} - \frac{1}{p_2}\right) \frac{1}{1-P_B} \int_0^{\min(\tilde{x}(y),\, B)} \left(y - \frac{p_2 x}{c}\right) f_0(x)\,dx$$

plus, when $y \geq p_2 B/c$, the atom contribution $\pi_{0,B}\left(y - p_2 B/c\right)/(1-P_B)$.

**Slowdown factor.** Define $S(y) = \mathbb{E}[T \mid Y=y]\,/\,(y/p_2)$. Substituting $f_0(x) = e^{\theta x}/K$:

**Case 1: $y < p_2 B / c$** (battery may survive, $\tilde{x}(y) < B$):

$$S(y) = 1 + \frac{p_2-p_1}{p_1\, y\, K(1-P_B)\,\theta}\left[\frac{p_2}{c\,\theta}\left(e^{\theta cy/p_2}-1\right) - y\right]$$

**Case 2: $y \geq p_2 B / c$** (battery always depletes, $\tilde{x}(y) \geq B$):

$$S(y) = 1 + \frac{p_2-p_1}{p_1\,y\,K(1-P_B)}\left\{\frac{1}{\theta}\left[\left(y-\frac{p_2 B}{c}\right)e^{\theta B} - y + \frac{p_2}{c\,\theta}(e^{\theta B}-1)\right] + \frac{p_0}{\lambda}\,e^{\theta B}\left(y-\frac{p_2 B}{c}\right)\right\}$$

**Properties of $S(y)$:**
- $S(y) \geq 1$ always, since $p_2 > p_1$.
- $S(y) \to 1$ as $y \to 0$: small demands are served entirely in Phase 1.
- $S(y) \to p_2/p_1$ as $y \to \infty$: large demands are served mostly in Phase 2.
- $S(y)$ is increasing in $y$: larger demands suffer more slowdown.
- $S(y) \to 1$ as $B \to \infty$ for any fixed $y$: a large battery eliminates slowdown entirely.
- $S(y) = p_2/p_1$ for all $y$ when $B \to 0$: no battery means all service at rate $p_1$.

Figure 1 shows $S(y)$ for $p_0=p_1=40$, $p_2=120$, $B=30$, $\mu=1/20$ and $\lambda \in \{1, 3, 5\}$, giving $\theta \in \{0.05, 0, -0.05\}$ and covering the lightly loaded, critical, and heavily loaded regimes.

![Slowdown factor $S(y)$ for three loading regimes](slowdown.png)

---

## 4. Extension: Multi-Server Loss System

### 4.1 Model

We now consider $s \geq 1$ servers sharing a single battery buffer of capacity $B$ and a single grid connection of rate $p_1$. Arriving cars are lost if all $s$ servers are busy.

**Assumptions.** We work under the following conditions:

- $p_1 < s\,p_2$: the grid alone cannot serve all $s$ servers simultaneously — the battery is necessary at full occupancy.
- $p_0 = p_1$: the idle charging rate equals the grid supply rate. We discuss the case $p_0 \neq p_1$ later.

**State space.** $Z(t) = (n(t), x(t))$ with $n(t) \in \{0,1,\ldots,s\}$ and $x(t) \in [0,B]$. The state space $E = \{0,\ldots,s\} \times [0,B]$ is compact.

**Vector field.** When $n$ servers are busy, each car charges at rate $p_2$ and the grid supplies $p_1$ in total. The net battery rate is:

$$r_n = p_1 - n\,p_2$$

This is positive for $n < p_1/p_2$, zero at $n = p_1/p_2$, and negative for $n > p_1/p_2$. Since $p_1 < s\,p_2$, the battery always drains at full occupancy $n=s$. Define the critical occupancy $n^* = \lfloor p_1/p_2 \rfloor \in \{0,\ldots,s-1\}$.

The full vector field is:

$$\dot{x} = \begin{cases} r_n & n \geq 0,\ 0 < x < B,\ r_n \neq 0 \\ 0 & r_n = 0,\ \text{or } x=0 \text{ and } r_n < 0,\ \text{or } x=B \text{ and } r_n > 0 \end{cases}$$

**Jumps.** Arrivals at rate $\lambda$ when $n < s$: $n \to n+1$, $x$ unchanged. When $x > 0$, departures at total rate $n\mu p_2$: $n \to n-1$, $x$ unchanged. When $x = 0$, the grid supplies $p_1$ shared equally among $n$ servers, giving total departure rate $n \cdot \mu \cdot (p_1/n) = \mu p_1$: $n \to n-1$, $x$ unchanged.

**Atoms.** The stationary distribution has:
- Atoms $\pi_n^{(B)}$ at $(n,B)$ for $n \leq n^*$ (flow stops at $x=B$ when charging)
- Atoms $\pi_n^{(0)}$ at $(n,0)$ for $n > n^*$ (flow stops at $x=0$ when draining)

### 4.2 Stationarity Conditions

**Interior equations** for $x \in (0,B)$, $0 \leq n \leq s$:

$$r_n\,f_n'(x) = -(\lambda\,\mathbf{1}_{n<s} + n\mu p_2)\,f_n(x) + \lambda\,f_{n-1}(x)\,\mathbf{1}_{n>0} + (n+1)\mu p_2\,f_{n+1}(x)\,\mathbf{1}_{n<s} \tag{En}$$

### 4.3 Solution

**Telescoping.** Sum (En) from $k=0$ to $k=n$ for $0 \leq n \leq s-1$. The right side telescopes:

$$\sum_{k=0}^n r_k\,f_k'(x) = -\lambda\,f_n(x) + (n+1)\mu p_2\,f_{n+1}(x) \tag{Tn}$$

**Decoupling ansatz.** Motivated by the single-server case, we try $f_n(x) = a_n\,e^{\theta x}$ for all $n$, with $a_0 = 1$ and $\theta$ to be determined. Substituting into (Tn):

$$\theta\,e^{\theta x}\sum_{k=0}^n r_k\,a_k = e^{\theta x}\left[-\lambda\,a_n + (n+1)\mu p_2\,a_{n+1}\right]$$

This gives the **recurrence**:

$$\boxed{a_{n+1} = \frac{\theta S_n + \lambda\,a_n}{(n+1)\mu p_2}} \qquad n = 0, 1, \ldots, s-1 \tag{R}$$

where $S_n = \sum_{k=0}^n r_k\,a_k$, with $S_{-1} = 0$.

For each value of $\theta$, (R) determines $a_1, \ldots, a_s$ uniquely given $a_0 = 1$.

**Special case $\theta = 0$.** The recurrence gives $a_{n+1} = \frac{\lambda\,a_n}{(n+1)\mu p_2} = \frac{\rho}{n+1}\,a_n$ where $\rho = \lambda/(\mu p_2)$, so:

$$a_n = \frac{\rho^n}{n!}$$

This is the **Erlang-B** form — the stationary distribution is uniform in $x$ with Poisson weights in $n$. We will see that $\theta = 0$ corresponds to a specific relationship between parameters.

**Consistency with $s=1$.** For $s=1$ and $p_0=p_1$, the recurrence gives $a_1 = (\theta p_1+\lambda)/(\mu p_2)$. The single-server boundary condition forces $a_1 = p_1/(p_2-p_1) = p_1/c$, which gives $\theta = \mu p_2/c - \lambda/p_1$ — recovering the base model exactly. $\checkmark$

**Determining $\theta$.** The parameter $\theta$ is determined by the boundary conditions. At $x=B$ for each $n \leq n^*$ (where $r_n > 0$), flux balance gives:

$$r_n\,f_n(B^-) = \lambda\,\pi_{n-1}^{(B)}\,\mathbf{1}_{n>0} - \lambda\,\pi_n^{(B)} + (n+1)\mu p_2\,\pi_{n+1}^{(B)}\,\mathbf{1}_{n<n^*}$$

At $x=0$ for each $n > n^*$ (where $r_n < 0$), flux balance gives:

$$|r_n|\,f_n(0^+) = \lambda\,\pi_{n-1}^{(0)}\,\mathbf{1}_{n>0} - (\lambda\,\mathbf{1}_{n<s}+\mu p_1)\,\pi_n^{(0)} + \mu p_1\,\pi_{n+1}^{(0)}\,\mathbf{1}_{n<s}$$

Together with normalization, these conditions determine $\theta$, the atoms, and the normalization constant.

**Connection to the mean drift.** In the single-server case, $\theta = \mu p_2/c - \lambda/p_1$ had a clear physical interpretation as the drift of the unconstrained battery level. The same interpretation extends to the multi-server case.

Consider the system without battery constraints — the occupancy $n(t)$ evolves as an M/M/s/s queue with stationary distribution $\pi_n^{\text{Erlang}} \propto \rho^n/n!$. The mean battery drift in this unconstrained system is:

$$\bar{r} = \sum_{n=0}^s r_n\,\pi_n^{\text{Erlang}} = p_1 - \frac{\lambda}{\mu}(1-P_B^{\text{Erlang}})$$

where $P_B^{\text{Erlang}} = \dfrac{\rho^s/s!}{\sum_{k=0}^s \rho^k/k!}$ is the Erlang-B blocking probability and we used $\sum_{n=0}^s n\pi_n^{\text{Erlang}} = \rho(1-P_B^{\text{Erlang}})$.

The sign of $\bar{r}$ determines the qualitative behavior:

- $\bar{r} > 0$: battery tends to fill in steady state, $\theta > 0$, mass concentrated near $x=B$
- $\bar{r} < 0$: battery tends to drain in steady state, $\theta < 0$, mass concentrated near $x=0$
- $\bar{r} = 0$: critical balance, $\theta = 0$, uniform density — Erlang-B weights

Indeed, $\theta = 0$ corresponds to $\bar{r} = 0$: when $a_n = \rho^n/n!$, one can verify that $S_{s-1}|_{\theta=0} \propto \bar{r}$ evaluated at the Erlang-B distribution. The exact value of $\theta$ for $\bar{r} \neq 0$ is the root of a transcendental equation involving $a_n(\theta)$ from the recurrence (R) and the boundary conditions — it does not have a simple closed form for general $s$, but reduces to $\theta = \mu p_2/c - \lambda/p_1$ for $s=1$.

**Complete stationary distribution:**
$$f_n(x) = C\,a_n(\theta)\,e^{\theta x}, \qquad x \in (0,B),\ n = 0,\ldots,s$$

where $C$ is the normalization constant and $a_n(\theta)$ are given by the recurrence (R).

### 4.4 Discussion

**Blocking probability.** By PASTA:

$$P_B = C\int_0^B a_s\,e^{\theta x}\,dx + \pi_s^{(0)}\,\mathbf{1}_{n^*<s} = C\,a_s\,\frac{e^{\theta B}-1}{\theta} + \pi_s^{(0)}$$

**Comparison with Erlang-B.** When $\theta = 0$ (no battery effect in steady state), $a_n = \rho^n/n!$ and the blocking probability reduces to the classical Erlang-B formula:

$$P_B^{\text{Erlang}} = \frac{\rho^s/s!}{\sum_{n=0}^s \rho^n/n!}$$

For $\theta \neq 0$ the blocking probability deviates from Erlang-B — the battery shifts mass toward high or low battery levels, affecting the service rates and hence the effective throughput.

**The role of $p_0 = p_1$.** This assumption ensures that $r_n = p_1 - np_2$ is independent of whether the system is idle or busy at a given $n$, simplifying the structure. We discuss relaxing this to $p_0 \neq p_1$ next.

### 4.5 The Case $p_0 \neq p_1$

When $p_0 \neq p_1$, the vector field at $n=0$ is $\dot{x} = p_0$ rather than $p_1$. This breaks the uniform structure $r_n = p_1 - np_2$ for $n \geq 1$: the $n=0$ equation has drift $p_0$ while the $n \geq 1$ equations have drift $p_1 - np_2$.

The telescoping trick still works but the summed equation for $n \geq 1$ now has a different $k=0$ term. Specifically, (Tn) for $n \geq 1$ becomes:

$$p_0\,f_0'(x) + \sum_{k=1}^n r_k\,f_k'(x) = -\lambda\,f_n(x) + (n+1)\mu p_2\,f_{n+1}(x)$$

The ansatz $f_n(x) = a_n e^{\theta x}$ still works, but now the $n=0$ drift $p_0$ and the $n \geq 1$ drifts $r_k = p_1 - kp_2$ are distinct. The recurrence becomes:

$$a_{n+1} = \frac{\theta(p_0\,a_0 + \sum_{k=1}^n r_k\,a_k) + \lambda\,a_n}{(n+1)\mu p_2}$$

This is still a first-order recurrence in $a_n$, just with a modified $S_n$. The single-server base model ($s=1$, $p_0 \leq p_1$) is recovered as a special case. $\checkmark$



---

## 5. Extension: The Case of Waiting Vehicles

We now remove the loss assumption and allow cars to queue. We consider the general case of $s \geq 1$ servers sharing a battery buffer of capacity $B$ and a grid connection of rate $p_1$, with $p_0 = p_1$ and $p_1 < sp_2$, as in Section 4 — but now arriving cars that find all servers busy **wait** rather than being lost.

### 5.1 Model

**State space.** The state is $Z(t) = (q(t), x(t))$ where:

- $q(t) \in \{0, 1, 2, \ldots\}$: number of cars in the system (in service or waiting)
- $x(t) \in [0, B]$: battery level

The state space is $E = \{0,1,2,\ldots\} \times [0,B]$, which is **not compact**.

**Vector field.** The number of active servers is $\min(q,s)$. The battery dynamics follow the same convention as Section 4:

$$\dot{x} = r_{\min(q,s)} = p_1 - \min(q,s)\cdot p_2$$

with $\dot{x} = 0$ when $x=0$ and $r_{\min(q,s)} < 0$, or $x=B$ and $r_{\min(q,s)} > 0$.

**Jumps.** Arrivals at rate $\lambda$ increment $q$ by 1. Departures occur at total rate $\min(q,s)\cdot\mu p_2$ when $x>0$, and at rate $\mu p_1$ when $x=0$ (the grid supplies $p_1$ shared equally among $\min(q,s)$ servers, giving total departure rate $\mu p_1$ regardless of occupancy). Departures decrement $q$ by 1.

### 5.2 Stability

**Theorem.** *A stationary distribution exists if and only if $\lambda < \mu p_1$, for any $s \geq 1$.*

**Sufficiency: $\lambda < \mu p_1$.** Take $V(q,x) = q$. For $q \geq 1$, we have $\min(q,s) \geq 1$ and $p(x) \geq p_1$, so:

$$\mathcal{A}V(q,x) = \lambda - \min(q,s)\,\mu p(x) \leq \lambda - \mu p_1 < 0$$

So $\mathcal{A}V \leq -\varepsilon$ outside the compact set $\{q=0\} \times [0,B]$, and the Foster-Lyapunov criterion gives stability. $\square$

**Necessity: $\lambda > \mu p_1$.** Suppose for contradiction that a stationary distribution $\pi$ exists. We derive two equations from stationarity $\int \mathcal{A}V\,d\pi = 0$.

**Flow balance** (taking $V(q,x) = q$):

$$\lambda = \mu p_1\,\Pr_\pi(q\geq 1,\,x=0) + \mu p_2\,\mathbb{E}_\pi[\min(q,s)\,\mathbf{1}_{q\geq 1,\,x>0}] \tag{F}$$

**Energy balance** (taking $V(q,x) = x$; jump terms vanish since $x$ is unchanged at jumps, and $|\mathcal{A}x| \leq sp_2$ is bounded so the application is justified):

$$p_0\,\Pr_\pi(q=0,\,x<B) = \mathbb{E}_\pi[(\min(q,s)\,p_2 - p_1)\,\mathbf{1}_{q\geq 1,\,0<x<B}] \tag{E}$$

Multiplying (E) by $\mu$ and rearranging:

$$\mu p_2\,\mathbb{E}_\pi[\min(q,s)\,\mathbf{1}_{q\geq 1,\,x>0}] = \mu p_1\,\Pr_\pi(q\geq 1,\,0<x<B) + \mu p_0\,\Pr_\pi(q=0,\,x<B)$$

Substituting into (F):

$$\lambda = \mu p_1\,\Pr_\pi(q\geq 1,\,x=0) + \mu p_1\,\Pr_\pi(q\geq 1,\,0<x<B) + \mu p_0\,\Pr_\pi(q=0,\,x<B)$$

$$= \mu p_1\,\Pr_\pi(q\geq 1) + \mu p_0\,\Pr_\pi(q=0,\,x<B)$$

$$\leq \mu p_1\,\Pr_\pi(q\geq 1) + \mu p_1\,\Pr_\pi(q=0) = \mu p_1$$

using $p_0 \leq p_1$. This contradicts $\lambda > \mu p_1$. $\square$

**Remark.** The stability threshold $\lambda < \mu p_1$ is independent of $s$, $B$, and $p_2$ — adding more servers or a larger battery does not improve stability. The proof uses only the flow balance, the energy balance, and $p_0 \leq p_1$. The bottleneck is always the grid rate $p_1$ at $x=0$, which caps total service capacity regardless of the number of servers. The battery improves transient performance but provides no stability benefit.



### 5.3 Stationary Distribution

The stationary distribution does not admit a product-form solution. The interior equations admit an exponential ansatz $f_q(x) = a_q e^{\theta x}$, but the resulting sequence $a_q$ is not non-negative for all $q$, and the atom balance equations at $x=0$ are incompatible with the interior structure. The coupling between the queue length and the battery level through the state-dependent service rate is the fundamental obstruction, and it applies for any $s\geq 1$. Numerical methods are required — for instance via the matrix-analytic approach for fluid-driven queues (Kulkarni 1997, Rogers 1994).

**Bounds on $\mathbb{E}[q]$.** Since $p_1 \leq p(x) \leq p_2$ for all $x$, the system is stochastically dominated by an M/M/1 queue at rate $\mu p_1$ from above and by an M/M/1 at rate $\mu p_2$ from below:

$$\frac{\lambda}{\mu p_2 - \lambda} \leq \mathbb{E}[q] \leq \frac{\lambda}{\mu p_1 - \lambda}$$

The lower bound corresponds to always charging at rate $p_2$ (infinite battery), the upper bound to always charging at rate $p_1$ (no battery).


---

## 6. Optimal Sizing

We return to the loss system of Sections 3 and 4. The operator controls the infrastructure and wishes to set it optimally. The mean demand $1/\mu$, fast charging rate $p_2$, and demand function $\lambda(\nu)$ are all exogenous. We set $p_0 = p_1$ throughout.

### 6.1 Problem Formulation

**Revenue.** The operator charges customers at price $\nu$ per unit energy delivered. The energy delivered per unit time in stationarity is $\lambda(\nu)(1-P_B)/\mu$. Revenue per unit time is $\nu \cdot \lambda(\nu)(1-P_B)/\mu$.

**Costs.** The operator pays:
- $c_g \cdot p_1$ per unit time for the grid connection
- $c_B \cdot B$ per unit time for the battery (amortized capital cost)
- $c_s \cdot s$ per unit time for the chargers (amortized capital + operating cost per server)

**Free variables.** The operator chooses:
- $s \in \mathbb{N}$: number of servers
- $p_1 \in (0, sp_2)$: grid connection rate
- $B \geq 0$: battery capacity
- $\nu \geq 0$: charging price per unit energy

The demand function $\lambda = \lambda(\nu)$ is decreasing in $\nu$ — a higher price reduces arrivals, which reduces $P_B$, coupling pricing and infrastructure decisions.

**Optimization problem:**

$$\max_{s \in \mathbb{N},\, p_1 \in (0,sp_2),\, B \geq 0,\, \nu \geq 0} \quad \nu \cdot \frac{\lambda(\nu)(1-P_B(s, p_1, B, \lambda(\nu)))}{\mu} - c_g\,p_1 - c_B\,B - c_s\,s \tag{OPT}$$

where $P_B$ depends on $s$, $p_1$, $B$, and $\lambda(\nu)$ through the multi-server stationary distribution of Section 4.

### 6.2 Structure of the Problem

**Relaxing $s$ to a continuous variable.** Replace $s \in \mathbb{N}$ with $s > 0$ real. The integer solution can be recovered by comparing $\lfloor s^* \rfloor$ and $\lceil s^* \rceil$.

**Monotonicity.** $P_B$ is decreasing in each of $s$, $p_1$, $B$ and increasing in $\lambda$. The revenue term $\nu\lambda(\nu)(1-P_B)/\mu$ is therefore increasing in $s$, $p_1$, $B$ and has ambiguous dependence on $\nu$ through both the direct effect and the indirect effect via $\lambda(\nu)$.

**Convexity.** The cost terms $c_g p_1 + c_B B + c_s s$ are linear. However, the revenue term is in general **neither convex nor concave**:

- As a function of $B$ or $p_1$ alone: $P_B$ decreases from a positive value toward 0, exhibiting diminishing returns — suggesting concavity of the revenue in these variables, but this is not guaranteed without explicit verification.
- As a function of $\lambda$: $P_B$ increases from 0 to 1 in an S-shaped curve — neither convex nor concave globally.
- As a function of $\nu$: the dependence goes through $\lambda(\nu)$ and then through $P_B(\lambda)$, inheriting the non-convexity.

The problem (OPT) is therefore a **nonlinear non-convex optimization** in general. First-order conditions are necessary but not sufficient for a global optimum.

**Necessary conditions.** At any interior local optimum:

$$\frac{\nu\lambda(\nu)}{\mu}\cdot\left(-\frac{\partial P_B}{\partial B}\right) = c_B, \qquad \frac{\nu\lambda(\nu)}{\mu}\cdot\left(-\frac{\partial P_B}{\partial p_1}\right) = c_g, \qquad \frac{\nu\lambda(\nu)}{\mu}\cdot\left(-\frac{\partial P_B}{\partial s}\right) = c_s$$

Each condition equates the marginal revenue from reduced blocking to the marginal cost of the respective resource. The pricing condition is:

$$\frac{d}{d\nu}\left[\nu\lambda(\nu)(1-P_B)\right] = 0$$

which is a standard monopoly pricing condition modified by the congestion effect.

*[Further analysis — special cases, structural properties, numerical solutions — to be developed.]*





---

## 7. Open Problems

- **Heavy-traffic limits.** In the queue model (Section 5), the stability threshold is $\lambda_c = \mu p_1$. The interesting heavy-traffic regime requires the battery to scale with the traffic intensity: let $\varepsilon = \mu p_1 - \lambda \to 0$ and $B = B(\varepsilon) \to \infty$ at an appropriate rate $B \sim \varepsilon^{-\alpha}$. With the correct scaling, the battery remains non-empty on the diffusion time scale and produces a non-trivial correction to the standard M/M/1 heavy-traffic limit. The candidate limit is a two-dimensional diffusion $(Q(t), X(t))$ where $Q(t)$ is a reflected Brownian motion with drift $-\varepsilon$ and $X(t)$ is driven by $Q(t)$ — charging when $Q=0$, draining at rate $\min(Q,s)\cdot p_2 - p_1$ when $Q \geq 1$. Identifying the correct scaling exponent $\alpha$ and characterizing the limit process is an open problem.
- **General demands and numerical methods.** For non-exponential demand distributions $Y \sim F$, the state space must be augmented with attained service, leading to a Volterra integro-differential equation for the stationary density. Closed-form solutions are available for phase-type distributions via matrix-exponential methods; for general $F$ the equation must be solved numerically via Volterra solvers. The connection to the M/G/1+G abandonment queue literature may yield further tractable cases.

---

## References

- Davis, M. H. A. (1984). Piecewise-deterministic Markov processes: a general class of non-diffusion stochastic models. *Journal of the Royal Statistical Society: Series B*, 46(3), 353–388.
- Davis, M. H. A. (1993). *Markov Models and Optimization*. Chapman & Hall.
- Kulkarni, V. G. (1997). Fluid models for single buffer systems. In *Frontiers in Queueing*, CRC Press.
- Meyn, S. P. and Tweedie, R. L. (1993). *Markov Chains and Stochastic Stability*. Springer.
- Rogers, L. C. G. (1994). Fluid models in queueing theory and Wiener-Hopf factorization of Markov chains. *Annals of Applied Probability*, 4(2), 390–413.
- Sharpe, M. (1988). *General Theory of Markov Processes*. Academic Press.

# Biomineralization

The conceptual model for biomineralization follows the one presented by [@Ebigbo2012].
It accounts for two-phase multi-component reactive transport on the continuum scale, including biofilm and calcite as solid phases.
The reactions considered are pH-dependent dissociation reactions, microbial growth, and decay as well as microbially catalyzed ureolysis and mass-transfer reactions between the different phases.
A mass transfer may occur between both fluid phases by the mutual dissolution of
water and CO<sub>2</sub> in the gas or the aqueous phase. It may also occur between
the aqueous phase and the two *solid* phases,
biofilm and calcite denoted by subscripts ($$f$$) and ($$c$$) respectively, by attachment or
detachment of biomass and precipitation or dissolution of calcite.

The mobile components, denoted by superscripts $$\kappa$$, are water ($$w$$),
dissolved inorganic carbon ($$\mathrm{C_{tot}}$$),
sodium ($$Na$$), chloride ($$Cl$$), calcium ($$Ca$$), urea ($$u$$),
ammonium and ammonia ($$\mathrm{N_{tot}}$$),
substrate ($$s$$), oxygen (O<sub>2</sub>), and suspended biomass ($$b$$).
A substrate is the carbon and energy source of the bacteria
and O<sub>2</sub> is the electron acceptor.


The primary variables solved are the aqueous-phase pressure $$p_\mathrm{w}$$,
mole fractions $$x^\kappa_\mathrm{w}$$ of component $$\kappa$$ in the water phase, and
for the solid phases biofilm and calcite, volume fractions $$\phi_\lambda$$.
Although the previous model by [@Ebigbo2012] was formulated using mass fractions $$X^\kappa_\mathrm{w}$$ as primary variables, this description is formulated using mole fractions  $$x^\kappa_\mathrm{w}$$.

All calcium carbonate is assumed to precipitate as calcite, since experimental
investigations of \citet{Phillips2015a, Mitchell2013, Lauchnor2013, Cuthbert2012}
confirmed by XRD measurements that calcite is the predominant polymorph of calcium carbonate precipitates forming under MICP conditions. In \citet{Phillips2013a}, calcite and possibly vaterite were observed.

However, the CO<sub>2</sub>-phase saturation is used as the primary variable instead of the mole fraction of total inorganic carbon in water $$x^\mathrm{C_{tot}}_\mathrm{w}$$
whenever both fluid phases are present within the same control volume \citep{Class2002}.
All reactive and mass-transfer processes are incorporated in the mass balance equations for the components Eq. [1](#mjx-eqn-eq:mb_transport) and Eq. [2](#mjx-eqn-eq:eq:mb_solid) by component-specific source and sink terms:

$$\begin{equation}
\sum\limits_{\alpha} \left[\frac{\partial}{\partial t}\left(\phi \rho_\mathrm{\alpha,\,mol} x^\kappa_\alpha S_\alpha \right) + \nabla\cdotp \left(\rho_\mathrm{\alpha,\,mol} x^\kappa_\alpha \mathbf{v}_\alpha \right) - \nabla\cdotp \left(\rho_\mathrm{\alpha,\,mol} \mathbf{D}^\kappa_\mathrm{pm,\alpha} \nabla x^\kappa_\alpha \right) \right] = q^\kappa,\:\alpha\in \mathrm{\{n;w\}} .
\label{eq:mb_transport}\tag{1}
\begin{equation}$$

Here, $$t$$ is time, $$\phi$$ porosity, $$\rho_\mathrm{\alpha,\,mol}$$, $$S_\alpha$$,
and $$\mathbf{v}_\alpha$$ the molar density, saturation and the velocity of phase $$\alpha$$ respectively,
$$x^\kappa_\alpha$$ the mole fraction of component $$\kappa$$ in phase $$\alpha$$.
$$\mathbf{D}_\mathrm{pm,\alpha}$$ is the dispersion tensor of phase $$\alpha$$ in the porous medium,and $$q^\kappa$$ is the source term of component $$\kappa$$ due to biochemical reactions. However, all components except water, CO<sub>2</sub> , and O<sub>2</sub> are assumed to be restricted to the water phase.

The mass balances for the solid phases calcite and biofilm contain
only a storage and source term since they are immobile:

$$\begin{equation*}
\frac{\partial}{\partial t} \left(\phi_\lambda \rho_\lambda \right) = q^\lambda,\:\lambda\in \mathrm{\{c;f\}}.
\label{eq:mb_solid}\tag{2}
\end{equation}$$

Here, $$\phi_\lambda$$ and $$\rho_\lambda$$ are volume fraction and mass density of
the solid phase $$\lambda$$, and $$q^\lambda$$ is the source term of phase $$\lambda$$ due to biochemical reactions.
The sources and sinks due to reactions $$q^\kappa$ and $q^\lambda$$ are specific to the components and are discussed in details in the subsequent section.

## Component-specific reactive source and sink terms

The source and sink terms account for the biogeochemical reactions occurring during MICP and the presence of CO<sub>2</sub>:
ureolysis, calcite precipitation, and dissolution, biomass growth under consumption of oxygen and substrate, biomass decay, as well as attachment and detachment of biomass.

## Water, sodium and chloride
Sodium and chloride do not participate in the reactions and water is the solute and is abundant, which is why its consumption by the hydrolysis of urea (Eq. [3](#mjx-eqn-eq:q_w_na_cl)) is considered negligible.
Thus, the reactive source terms for water $$q^\mathrm{w}$$, sodium $$q^\mathrm{Na}$$ and chloride $$q^\mathrm{Cl}$$ are zero:

$$\begin{equation*}
q^\mathrm{w} =   q^\mathrm{Na} = q^\mathrm{Cl} = 0
\label{eq:q_w_na_cl}\tag{3}
\end{equation}$$

## Urea and total nitrogen
The source term for $$\mathrm{N_{tot}}$$, $$q^\mathrm{N_{tot}}$$, and the sink term for urea $$q^\mathrm{u}$$ result from ureolysis (Eq. [5](#mjx-eqn-eq:q_ntot)).
For each mole of urea hydrolyzed, 2 moles of $$\mathrm{N_{tot}}$$ are generated. The $$q^\kappa$$ are thus:

$$\begin{equation*}
q^\mathrm{u}=-r_\mathrm{urea},
\label{eq:q_u}\tag{4}
\end{equation}$$

$$\begin{equation*}
q^\mathrm{N_{tot}}=2r_\mathrm{urea},
\label{eq:q_ntot}\tag{5}
\end{equation}$$

where $$r_\mathrm{urea}$$ is the ureolysis rate from [@Ebigbo2012].
The rate equation is derived from the complex enzymatic rate equation originally published by [@Fidaleo2003] for the ureolysis rate of pure urease extracted from jack bean ($$\textit{Canavalia ensiformis}$$) seeds.

It accounts for enzyme inactivation due to suboptimal pH, inhibition caused by high product (NH<sup>4+</sup>) concentrations and a Michaelis-Menten type reaction rate dependency of the urea concentration:

In batch experiments, [@Lauchnor2015] investigated the influences of urea,
NH<sup>4+</sup>, cell concentration, and pH of the medium on the ureolysis of whole cells of $$\textit{S.~pasteurii}$$.
These new observations and measurements can be accounted for in the MICP model by adjusting the Michaelis-Menten kinetics and parameters:

$$\begin{equation*}
r_\mathrm{urea} = k_\mathrm{urease}
\:k_\mathrm{ub}\:\rho_\mathrm{f}\:\phi_\mathrm{f}
\:\frac{m^\mathrm{u}}{m^\mathrm{u}+K_\mathrm{u}}.
\label{eq:rureol_new}
\end{equation}$$

$$r_\mathrm{urea}$$ represents the revised rate of ureolysis according to [@Lauchnor2015],
$$k_\mathrm{urease}$$ the revised maximum activity of urease adapted from  [@Lauchnor2015],
$$\rho_\mathrm{f}$$ and $$\phi_\mathrm{f}$$ the density and volume fraction of biofilm respectively, $$k_\mathrm{ub}$$ the mass ratio of urease to biofilm as given in \citet{Bachmeier2002},
$$m^\mathrm{u}$$ the molality of urea calculated from the water phase composition,
and $$K_\mathrm{u}$$ is the half saturation constant for urea adapted from [@Lauchnor2015].

To convert the new half-saturation constant given by \citet{Lauchnor2015}
$K = 355\:\unitfrac{mmol}{l}$ to
molalities $\unitfrac{mol}{kg_\mathrm{H_{2}O}}$, we assume that
the concentrations used in the experiments did not affect the density.
Thus, the half-saturation constant $K_\mathrm{u}$ in the
model is set to $0.355\: \unitfrac{mol}{kg_\mathrm{H_{2}O}}$.

The value of the apparent urease activity $k_\mathrm{urease}$ in a whole-cell system
is calculated by dividing the maximum reaction rate $V_\mathrm{max}$ by the urease
content of a bacterial cell $k_\mathrm{ub}$. According to \citet{Bachmeier2002},
the urease content in \textit{S.~pasteurii}
is at most 1\% of the dry cell mass.
As the cells in \citet{Lauchnor2015} are in late exponential phase with
urea-replete conditions,
it is assumed that they produce the maximum amount of urease possible.
Thus, for further analysis, it is assumed that $k_\mathrm{ub}$ is equal to 0.01.
To convert the units of $V_\mathrm{max}$ given in
\citet{Lauchnor2015} from $\unitfrac{mmol}{CFU \cdot h}$
(CFU standing for colony-forming units)
to SI units for compatibility with the conceptual model,
the same cell weight of $2.5\:10^{-16}\: \unitfrac{kg}{CFU}$ is used as in \citet{Ebigbo2012},
originally given in \citet{Norland1987}. The resulting updated value for
$k_\mathrm{urease}$ is $706.7 \:\unitfrac{mol}{kg_\mathrm{biomass} \cdot s}$,
calculated from the value of the maximum reaction rate
$V_\mathrm{max}=6.4\:10^{-9} \:\unitfrac{mmol}{CFU \cdot h}$ given in \citet{Lauchnor2015}.

As the observed ureolysis rates in \citet{Lauchnor2015} are directly proportional
to the cell concentration, the exponent previously fitted in Equation~\eqref{eq:rureol_old}
\citep{Ebigbo2012} is set to $n_\mathrm{ub}=1$ to achieve this linear dependence.

\subsubsection*{Calcium and calcite}\label{sec:source_ca_caco3}
The source terms of calcium $q^\mathrm{Ca}$ and calcite $q^\mathrm{c}$
are determined by the rates of precipitation and dissolution.
When the aqueous phase is oversaturated with respect to calcite, it precipitates according to (Eq.~\eqref{eq:caco3_precip_reaction}).
In the opposite case, calcite dissolves until the solution is saturated or all calcite is already dissolved:

\begin{equation}
q^\mathrm{Ca} = r_\mathrm{diss} - r_\mathrm{prec},
\label{eq:q_ca}
\end{equation}
\begin{equation}
q^\mathrm{c} = - r_\mathrm{diss} + r_\mathrm{prec}.
\label{eq:q_caco3}
\end{equation}

Here, $r_\mathrm{diss}$ is the rate of calcite dissolution and $r_\mathrm{prec}$ the rate of calcite precipitation.
Both reaction rates are calculated as follows, depending on the interfacial area available for the reaction
as well as the saturation index $\Omega$ and, in the case of the dissolution, additionally on the molality of \ce{H+}.
The precipitation rate of calcite is calculated as:

\begin{equation}
r_\mathrm{prec} = k_\mathrm{prec}A_\mathrm{sw}\left(\Omega-1\right)^{n_\mathrm{prec}};\:\text{for}\: \Omega\ge1,
\label{eq:r_prec}
\end{equation}
\begin{equation}
A_\mathrm{sw} = A_\mathrm{sw,0}\left( 1-\frac{\phi_\mathrm{c}}{\phi_0} \right)^\frac{2}{3},
\label{eq:Asw}
\end{equation}
\begin{equation}
\Omega = \frac{m^\mathrm{Ca^{2+}}\gamma^\mathrm{Ca^{2+}} m^\mathrm{CO_3^{2-}}\gamma^\mathrm{CO_3^{2-}}}{K_\mathrm{sp}},
\label{eq:Omega}
\end{equation}

where $k_\mathrm{prec}$ and $n_\mathrm{prec}$ are empirical precipitation parameters from \citet{Zhong1989},
$A_\mathrm{sw}$ and $A_\mathrm{sw,0}$ are the current and initial interfacial areas respectively
between the water phase and the solid phases,
$K_\mathrm{sp}$ the calcite solubility product and
$m^\mathrm{Ca^{2+}}$ and $m^\mathrm{CO_3^{2-}}$ the molalities of calcium and carbonate respectively.
The activity coefficients $\gamma^\kappa$ are calculated using
Pitzer equations \citep{Clegg1995, Millero1984, Wolf1989}.
The dissolution rate of calcite is calculated as:

\begin{equation}
r_\mathrm{diss} = \left(k_\mathrm{diss,1}m^\mathrm{H^{+}}+
k_\mathrm{diss,2}\right)A_\mathrm{cw}\left(\Omega-1\right)^{n_\mathrm{diss}};\:\text{for}\: \Omega<1,
\label{eq:r_diss}
\end{equation}
\begin{equation}
A_\mathrm{cw} = \mathrm{min} \left(A_\mathrm{sw}, a_\mathrm{c}\phi_\mathrm{c}\right),
\label{eq:Acw}
\end{equation}

where $k_\mathrm{diss,1}$, $k_\mathrm{diss,2}$, and $n_\mathrm{diss}$ are dissolution parameters \citep{Chou1989,Compton1989}
and $a_\mathrm{c}$ is the specific surface area and $\phi_\mathrm{c}$ the volume fraction of calcite.

However, in the present simplified chemistry system, dissolution is neglected, $r_\mathrm{diss}=0$,
and instead of calculating the complex geochemistry, e.g. dissociation of inorganic carbon into carbonate and bicarbonate,
it is assumed that the system has reached a steady state, were precipitation rate is equal to the ureolysis rate and thus it results:
\begin{equation}
r_\mathrm{prec} = r_\mathrm{urea}
\end{equation}

\subsubsection*{Dissolved inorganic carbon}\label{sec:source_ctot}
Dissolved inorganic carbon is generated by the hydrolysis of urea (Eq.~\eqref{eq:MICP_chem_reaction})
as well as by the dissolution of calcite while it is consumed by the precipitation of calcite (Eq.~\eqref{eq:caco3_precip_reaction}).
Thus, the source term of dissolved inorganic carbon $q^\mathrm{C_{tot}}$ results in:

\begin{equation}
q^\mathrm{C_{tot}} =  r_\mathrm{urea} + r_\mathrm{diss} - r_\mathrm{prec},
\label{eq:q_ctot}
\end{equation}


\subsubsection*{Suspended and attached biomass} \label{sec:source_bio}

The source and sink terms of suspended and attached biomass (biofilm), $q^\mathrm{b}$ and $q^\mathrm{f}$,
include four reaction rates each, corresponding to the biomass-related processes the model accounts for.
These processes are growth and decay increasing and decreasing the suspended or attached biomass as well as
attachment and detachment describing the transfer of biomass from the suspended to the attached state and vice versa:

\begin{equation}
q^\mathrm{b} = \frac{r^\mathrm{b}_\mathrm{g} - r^\mathrm{b}_\mathrm{b} - r_\mathrm{a} + r_\mathrm{d}}{M^\mathrm{b}},
\label{eq:q_b}
\end{equation}
\begin{equation}
q^\mathrm{f} = \frac{r^\mathrm{f}_\mathrm{g} - r^\mathrm{f}_\mathrm{b} + r_\mathrm{a} - r_\mathrm{d}}{M^\mathrm{f}}
\label{eq:q_f}
\end{equation}

where $r^\mathrm{b}_\mathrm{g}$ is the growth rate and $r^\mathrm{b}_\mathrm{b}$ the decay rate of suspended biomass,
$r_\mathrm{a}$ the attachment rate, $r_\mathrm{d}$ the detachment rate and
$M^\mathrm{b}$ the molar mass of biomass to convert the rates in the units
from mass to moles per volume and time.
Accordingly, $r^\mathrm{f}_\mathrm{g}$ and $r^\mathrm{f}_\mathrm{b}$ are the growth and decay of biofilm
and $M^\mathrm{f}$ is the molar mass of biofilm.
All rates influencing both attached and suspended biomass are assumed to be of a first-order type, where
the rate is calculated by the product of a specific rate and the respective biomass, which is
$C_\mathrm{w}^\mathrm{b}S_\mathrm{w}\phi$ in the case of suspended and
$\phi_\mathrm{f}\rho_\mathrm{f}$ in the case of attached biomass.
Here, $C_\mathrm{w}^\mathrm{b}$ is the suspended biomass mass concentration in the water phase,
$S_\mathrm{w}$ the water phase saturation.

The growth rates of suspended and attached biomass are as follows:

\begin{equation}
r_\mathrm{g}^\mathrm{b} = \mu_\mathrm{g} C_\mathrm{w}^\mathrm{b}S_\mathrm{w}\phi,
  \label{eq:r_g_susp}
\end{equation}
\begin{equation}
  r_\mathrm{g}^\mathrm{f} = \mu_\mathrm{g} \phi_\mathrm{f}\rho_\mathrm{f},
  \label{eq:r_g_f}
\end{equation}
with the specific growth rate $\mu_\mathrm{g}$.
$\mu_\mathrm{g}$ is calculated using double Monod kinetics to reproduce the dependency of
the microbial growth on both substrate and oxygen.

\begin{equation}
 \mu_\mathrm{g} = k_\mathrm{\mu,old} Y
 \frac{C_\mathrm{w}^\mathrm{s}}{K_\mathrm{s} + C_\mathrm{w}^\mathrm{s}}
 \frac{C_\mathrm{w}^\mathrm{O_2}}{K_\mathrm{O_2} + C_\mathrm{w}^\mathrm{O_2}}.
  \label{eq:mu}
\end{equation}

Here, $k_\mathrm{\mu,old}$ is the maximum specific growth rate, which was fitted in \citet{Ebigbo2012},
$Y$ the yield coefficient expressing the ratio of biomass generated to the mass of substrate consumed \citep{Seto1985}.
$C_\mathrm{w}^\mathrm{s}$ and $C_\mathrm{w}^\mathrm{O_2}$
are the mass concentrations of substrate and oxygen in the water phase and
$K_\mathrm{s}$ \citep{Taylor1990} and  $K_\mathrm{O_2}$ \citep{Hao1983}
are the half-saturation coefficients for substrate and oxygen respectively.

The decay rates are calculated similarly to the growth rates:

\begin{equation}
r_\mathrm{b}^\mathrm{b} = k_\mathrm{b}^\mathrm{b} C_\mathrm{w}^\mathrm{b}S_\mathrm{w}\phi,
   \label{eq:r_b_susp}
\end{equation}

\begin{equation}
 r_\mathrm{b}^\mathrm{f} =  k_\mathrm{b}^\mathrm{f}\phi_\mathrm{f}\rho_\mathrm{f},
  \label{eq:r_b_f}
\end{equation}

except that the specific decay rates of suspended and attached biomass,
$k_\mathrm{b}^\mathrm{b}$ and $k_\mathrm{b}^\mathrm{f}$ respectively,
take different processes into account, increasing inactivation.
For suspended biomass, non-optimal, acidic pH is assumed to increase inactivation:

\begin{equation}
k_\mathrm{b}^\mathrm{b} = b_0 \left(1 + \frac{K_\mathrm{pH}}{m_\mathrm{H^{+}} ^2} \right),
   \label{eq:k_b_susp}
\end{equation}

where $b_0$ is the endogenous decay rate \citep{Taylor1990}
and $K_\mathrm{pH}$ is an empirical constant from \citet{Kim2000},
accounting for increased cell inactivation at low pH.
High pH is assumed not to influence the inactivation of suspended biomass, as \textit{S. pasteurii} is alkaliphile.
On the contrary, attached cells are protected from harsh environmental conditions and the presence of \ce{CO2} within the biofilm
by protective mechanisms such as extracellular polymers \citep{Mitchell2008}.
However, as calcite precipitates mainly in or close to the biofilm, cells may be covered with calcite precipitates
or disrupted by crystals inactivating the affected cells \citep{Dupraz2009a, Whiffin2007}. %, Parks2009}.
Consequently, the precipitation rate is assumed to increase the specific decay rate of attached biomass:

\begin{equation}
 k_\mathrm{b}^\mathrm{f} = b_0 + \frac{r_\mathrm{prec} M^\mathrm{c}}{\rho_\mathrm{c}\left(\phi_0 - \phi_\mathrm{c}\right)},
  \label{eq:k_b_f}
\end{equation}

where $\frac{r_\mathrm{prec} M^\mathrm{c}}{\rho_\mathrm{c}}$ is the volumetric calcite precipitation rate,
$M^\mathrm{c}$ being the molar mass and $\rho_\mathrm{c}$ the density of calcite, and
$\phi_0 - \phi_\mathrm{c}$ is the space available for calcite precipitation, which may occur in the biofilm or the pore space.

The attachment rate $r_\mathrm{a}$ quantifies the biomass transfer from
the suspended to the attached state.
As attachment is modeled assuming a first-order kinetic rather than a
sorption-type behavior, it is independent of the amount of attached biomass:

\begin{equation}
r_\mathrm{a}= k_\mathrm{a} C^\mathrm{b}_\mathrm{w} \phi S_\mathrm{w},
\label{eq:r_a}
\end{equation}

where $k_\mathrm{a}$ is the specific attachment rate from \citet{Taylor1990}.
It is considered to
consist of two terms, $c_\mathrm{a,1}\phi_\mathrm{f}$ and $c_\mathrm{a,2}$,
whose coefficient $c_\mathrm{a,1}$ accounts for preferential attachment to existing biofilm
while $c_\mathrm{a,2}$ accounts for the unspecific attachment to arbitrary surfaces:

\begin{equation}
k_\mathrm{a}=c_\mathrm{a,1} \phi_\mathrm{f} + c_\mathrm{a,2}.
\label{eq:k_a}
\end{equation}

Detachment of biomass from biofilm is assumed to be proportional to the shear stress.
Additionally, the growth contributes to the detachment rate, as vigorously growing
biofilm is typically weaker and as such more susceptible to detachment \citep{Rittmann1982}.
As the model of \citet{Ebigbo2012} is defined on the Darcy scale but the shear stress is a micro-scale property,
it is approximated using the absolute value of the water-phase potential gradient:

\begin{equation}
r_\mathrm{d}=k_\mathrm{d} \phi_\mathrm{f} \rho_\mathrm{f},
\label{eq:r_d}
\end{equation}

where $r_\mathrm{d}$ is the rate of detachment and $k_\mathrm{d}$ the detachment coefficient:

\begin{equation}
k_\mathrm{d}=c_\mathrm{d}
\left( \phi S_\mathrm{w} \left|\nabla p_\mathrm{w} -  \rho_\mathrm{w} \mathbf{g} \right| \right)^{0.58}
+ \frac{\phi_\mathrm{f}}{\phi_0-\phi_\mathrm{c}} \mu_\mathrm{g}.
\label{eq:k_d}
\end{equation}

Here, $c_\mathrm{d}$ is a coefficient for the shear-stress-dependent detachment,
$\left|\nabla p_\mathrm{w} -  \rho_\mathrm{w} \mathbf{g} \right|$ the absolute value of
the water-phase potential gradient, and $\phi_0$ the initial porosity.
% $\phi_\mathrm{c}$ the volume fraction of calcite.
% $\mu$ the specific growth rate of biomass, calculated based on double Monod kinetics dependent
% on both substrate and oxygen, see Table~\ref{tab:q_table}.

% Since the detachment rate $r_\mathrm{d}$ is dependent on the potential gradient, it
% increases as the intrinsic permeability $K$ is reduced during MICP
% by bio-clogging and precipitated calcite.



\subsubsection*{Substrate and oxygen} \label{sec:source_sub_o2}
The consumption of substrate and oxygen is linked to the growth of both suspended and attached biomass
by the yield coefficient $Y$ of substrate:

\begin{equation}
q^\mathrm{s} = - \frac{r^\mathrm{b}_\mathrm{g} + r^\mathrm{f}_\mathrm{g}}{M^\mathrm{s}Y}.
\label{eq:q_s}
\end{equation}

In the case of oxygen, the coefficient $F$, which is the ration of oxygen consumed per substrate consumed,
is used to express the biomass yield per oxygen consumed using $Y$:

\begin{equation}
 q^\mathrm{O_2} = F q^\mathrm{s} = - F \frac{ r^\mathrm{b}_\mathrm{g} + r^\mathrm{f}_\mathrm{g}} {M^\mathrm{O_2}Y}.
 \label{eq:q_o2}
\end{equation}

\subsection{Supplementary equations}\label{sec:supplement_eq}

\subsubsection{Permeability and porosity}\label{sec:perm_poro_pw}
The permeability decreases due to biofilm growth and calcite precipitation as already discussed in
Chapter~\ref{ch:introduction}. In the model, the reduction of permeability is calculated based on the
reduction of porosity:

\begin{equation}
\frac{K}{K_\mathrm{0}}=\left(\frac{\phi}{\phi_\mathrm{0}}\right)^3 \left(\frac{1-\phi_0}{1-\phi}\right)^2.
\label{eq:perm}
\end{equation}

Here, $K_\mathrm{0}$ is the initial permeability, %$\phi$ the porosity,
and $\phi_0$ is the initial porosity. The porosity $\phi$ decreases as
the volume fractions of biofilm and calcite increase (Eq.~\eqref{eq:porosity_update_general}:

\begin{equation}
\phi=\phi_\mathrm{0}-\phi_\mathrm{c}-\phi_\mathrm{f}.
\label{eq:poro}
\end{equation}

\subsubsection{Capillary pressure and relative permeability}\label{sec:pc-Sw_kr-Sw_pw}
The capillary pressure saturation relation of Brooks and Corey \citep[e.g.][]{Brooks1964, Corey1994} is used to
calculate the capillary pressure based on the wetting phase saturation (Eq.~\eqref{eq:pc-Sw}), %see Section~\ref{sec:pc},
using an entry pressure $p_\mathrm{d} = 10^4$~Pa
and a pore-size distribution parameter $\lambda = 2$ \citep[see][]{Ebigbo2010}.
The relative permeabilities of the wetting and the non-wetting phase are also
calculated using the relations given by Brooks and Corey, see Section~\ref{sec:kr}, with
identical parameters as for the capillary pressure.

\subsubsection{Fluid properties}\label{sec:fluid_properties_pw}
The density and the viscosity of the \ce{CO2} phase are calculated using
the relation given by \citet{Span1996} and \citet{Fenghour1998} respectively.
In these calculations, the effects of the small amounts of water and oxygen are neglected.
The density and the viscosity of the aqueous phase are calculated
according to \citet{Batzle1992} as a function of salinity.
Sodium, chloride and calcium are considered to contribute to the salinity.


% \subsubsection{Mutual dissolution of \ce{CO2} and brine}\label{sec:phase_partition_pw}
\subsubsection{Phase partitioning of components}\label{sec:phase_partition_pw}
The dissolution of \ce{CO2} in the aqueous phase is calculated according
to \citet{Duan2003} as a function of temperature, pressure and salinity.
In the two-phase case, the concentration of \ce{H2CO3} is assumed to be
equal to the solubility while, in the case of the exclusive presence of the aqueous phase,
the concentration of inorganic carbon %in the aqueous phase
is exclusively dependent on the precipitation, dissolution and ureolysis reactions.
In this case, the solubility represents the maximum possible concentration.
The mass fraction of water in the \ce{CO2} phase is  assumed to be constant as in \citet{Bielinski2006}.
For the solubility of oxygen in the aqueous phase, Henry's law is used
with parameters according to \citet{Sander1999}, now published in \citet{Sander2015}.

\subsection{Diffusion, dispersion and tortuosity}\label{sec:disp_model}
The dispersion of components is accounted for in the model and
calculated according to Equation~\eqref{eq:Dispersion_coeff}.
The tortuosity is calculated according to \citet{Millington1961}:
\begin{equation}
\tau_\alpha=\frac{\left( \phi S_\alpha \right)^{\frac{7}{3}}}{\phi^2}.
\label{eq:tortuosity}
\end{equation}

\section{Microbially induced calcite precipitation}\label{sec:intro_MICP}
In this thesis, the biogeochemical engineering application investigated
is microbially induced calcite precipitation (MICP)
in the context of sealing possible leakage pathways of subsurface gas (\ce{CO2}, \ce{CH4}, \ce{H2}) reservoirs.
The necessary processes are two-phase multi-component reactive transport including precipitation and dissolution
of calcite as well as the biomass-related processes: attachment of biomass to surfaces,
detachment of biomass from a biofilm, and growth and decay of biomass.
Additionally, the reduction in porosity and permeability has to be considered;
this results from the presence of the solid phases biofilm and calcite in the pore space.
 
MICP offers an engineering option that
uses controlled biofilm growth to achieve targeted
calcite precipitation. In subsurface applications, this process is typically
associated with a reduction of porosity and, even more importantly, of permeability.
As an engineering technology, it can be used to alter hydraulic flow conditions
and can be applied, for example, to cut off highly permeable pathways such as
fractures, faults, or behind-casing defects in boreholes within a geological
formation \citep{Phillips2013a, Mitchell2013, Phillips2013c}.

\textit{S.~pasteurii}
expresses the enzyme urease that catalyzes the hydrolysis reaction of
urea (\ce{CO(NH2)2}) into ammonia (\ce{NH_3}) and carbon
dioxide (\ce{CO2}):

\begin{equation}
% \ce{CO(NH_2)_2 + 2H2O -> [{urease}] 2NH3 + H2CO3}.
\mathrm{CO(NH_2)_2} + 2\mathrm{H_2O} \xrightarrow{urease}%\longrightarrow
2\mathrm{NH_{3}} + \mathrm{H_2CO_{3}}.
\label{eq:ureolysis_reaction}
\end{equation}

Aqueous solutions of ammonia become alkaline as \ce{H+} is consumed until the equilibrium of ammonium and ammonia is reached.
Thus, the ureolysis reaction leads to an increase in pH until the pH is equal to the pKa of ammonia:

\begin{equation}
\ce{NH3 + H+ <=> NH4+}.
\label{eq:nh3_dissociation_reaction}
\end{equation}

This shifts the
carbonate balance in an aqueous solution toward higher concentrations of
dissolved carbonate (\ce{CO_3^2-}):

\begin{equation}
\ce{H2CO3 <=> HCO3- + H+},
\label{eq:h2co3_dissociation_reaction}
\end{equation}

\begin{equation}
\ce{HCO3- <=> CO3^2- + H+}.
\label{eq:hco3_dissociation_reaction}
\end{equation}

Adding calcium (\ce{Ca2+}) to the system then results
in the precipitation of calcium carbonate (\ce{CaCO3}):

\begin{equation}
\ce{CO3^2- + Ca^2+ -> CaCO3 v}.
\label{eq:caco3_precip_reaction}
\end{equation}

The resulting overall MICP reaction equation is:
\begin{equation}
\mathrm{CO(NH_2)_2} + 2\mathrm{H_2O} + \mathrm{Ca^{2+}} \xrightarrow{urease}%\longrightarrow
2\mathrm{NH_{4}^{+}} + \mathrm{CaCO_{3}}\downarrow.
\label{eq:MICP_chem_reaction}
\end{equation}

In a porous medium, this process,
which results in the aforementioned impacts on the hydraulic properties,
depends on the interplay between biofilm growth, fluid dynamics,
and reaction rates.
A pore-scale sketch of the most important processes of MICP is shown
in Figure~\ref{fig:pore}.

  \begin{figure}
  \centering % better than center because it does not start new line   
  \includegraphics[width=0.4\textwidth]{Pictures/pore-scale_w_processes_named.eps}
 \caption{Schematic view of relevant processes during MICP.}
 \label{fig:pore}
 \end{figure}
% Accordingly, a biofilm may develop which is characterized by an interaction
% of attachment, nutrient-dependent growth, decay, and detachment with the hydraulic conditions in the
% pore space. The density of a biofilm, i.e. the number of cells within a
% certain volume, is dependent on several factors
% such as shear stress or nutrient supply \citep{Paul2012};
% data for these factors are not yet available for a reliable quantitative description.

A major difficulty for practical engineering applications of
MICP is the predictive planning of its use and impact,
since it involves a number of complex interacting processes.
While the basic chemistry and the flow processes are known, it is the exact
quantitative description of the interactions and, in particular, the influence
of the biofilm and the developing precipitates
that pose challenges to achieving predictability.


However, for the sake of simplicity, the dissociation processes are not modeled in this example and it
is assumed that the system has reached a steady state were the calcite precipitation rate is equal to the rate of ureolysis.












INFORMATION REGARDING THE COLUMN EXPERIMENTAL SETUP


Data from three column experiments are used in the inverse modeling and
determination of the fitting parameters.
One column from the previous study, Column~4 (C4), is used to
ensure reproducibility. 
Two more column studies were performed with modifications to the previous study. 
As these columns were operated as duplicates, they are denoted column D1 and column D2.

All columns were constructed from clear PVC pipe of 2.54~cm diameter
and 61~cm length, filled with 40-mesh quartz sand (0.5~mm effective filtration size,
manufacturer information, Unimin Corp., Emmet, ID), packed under water and vertically positioned. 
They were inoculated with \textit{S.~pasteurii} and subjected to a 16~h attachment period,
after which a 24~h growth period was initiated by continuous addition of growth medium in
an upflow configuration at 10~$\unitfrac{ml} {min}$. 
Following the growth period, calcium additions commenced, consisting of injecting
two pore volumes of a medium containing calcium, followed by a no-flow period to allow reaction and
precipitation to occur. This mineralization period was followed by injection of growth
medium to resuscitate the microorganisms and the procedure was repeated until calcium
had been injected 22 times for column C4 and 30 times for columns D1 and D2.

Differences between the operating parameters for the columns C4, D1 and D2 are
presented in Table~\ref{tab:ureol_exp_setups},
while the composition of injected fluids is given in Table~\ref{tab:ureol_composition}.
The first modification to the columns D1 and D2 was a reduction in the $\mathrm{Ca^{2+}}$ concentration
to 0.33~M, which is equal to the molar concentration of urea in the mineralization medium.
Second, the no-flow period for biomineralization was reduced from 24~h to 4~h,
as successfully implemented by \citep{Phillips2013a} to increase $\mathrm{Ca^{2+}}$-precipitation
efficiency and decrease the time necessary to achieve the desired permeability reduction. 
Finally, the column studies used in the improved model calibration
were constructed with five sampling ports
at locations 10, 20, 30, 40 and 50~cm downstream of the inlet.
Ports were made by drilling 1.252~cm
holes in the PVC and securing rubber stoppers in the holes,
allowing for sampling of pore fluid via a 23~gauge needle and 5~ml syringe. 
The sampling took place along the direction of flow and
over time during the experiments, providing additional measurements for
model calibration.

While columns D1 and D2 were operated, pore fluid was sampled from the ports
at 1-hour intervals during selected no-flow mineralization periods. 
The effluent pH, $\mathrm{NH_{4}^{+}}$ and $\mathrm{Ca^{2+}}$ were monitored at the
end of each mineralization and growth period, in addition to the
$\mathrm{NH_{4}^{+}}$ and $\mathrm{Ca^{2+}}$ measured in sampling ports over time
during mineralization periods.
Viable bacterial counts were performed on effluent samples, ensuring that
\textit{S.~pasteurii} remained the active organism during the experiments.
Finally, the columns were cut into eight sections to measure the precipitated calcite.

Additionally, a 2D radial flow
experiment was conducted in a
radial flow reactor constructed from two Plexiglas plates and a
bicycle rim (BR), which enabled us to
investigate the influence of a radial flow field on MICP.
The radial flow experiment is valuable for investigating
differences between different fluid-dynamic cases and making the
model more broadly applicable to various environments, especially
for field-relevant radial flow conditions.

The radial flow reactor was constructed from  two 31~cm-diameter, $\nicefrac{1}{4}^{\prime\prime}$ Plexiglas plates
that sandwiched a pack of the 40-mesh quartz sand, which was sealed on the
outer diameter by squeezing the Plexiglas onto a metal ring (bicycle rim) with clamps. 
In between the Plexiglas and metal was a bead of Teflon putty which
sealed the edge to prevent leaks when clamped tightly.
Sixteen~$\nicefrac{1}{4}^{\prime\prime}$ NPT threads were inserted around the circumference of the top
Plexiglas plate (outlet) as well as one thread in the center (inlet);
these were all attached by barbed fittings to silicone tubing.
The central inlet was connected via tubing to a peristaltic pump;
the outlet tubes emptied into an elevated effluent collection tub.

The radial flow reactor was inoculated with 3~l of an overnight culture of \textit{S.~pasteurii}
which was allowed to attach for 16.5~hours and then subjected to a 12-hour growth stage. 
After the growth period, a medium containing calcium was injected under constant flow conditions
for 12~hours followed by an overnight precipitation period. 
At the termination of the experiment, the Plexiglas plates were removed and the sand was sampled
at 2.5, 4, 6.5 and 9 inches from the center in eight radial directions.
Relevant experimental design parameters are summarized in Table~\ref{tab:ureol_exp_setups},
while Table~\ref{tab:ureol_composition} gives the
composition of the injected fluids.

Microbially, whole-cell-catalyzed ureolysis is the driving force of
MICP, increasing pH to the pKa of $\mathrm{NH_{3}}$-$\mathrm{NH_{4}^{+}}$
by the production of $\mathrm{NH_{3}}$.
At this pH, a substantial amount of carbonate is present in
the solution (the pKa of $\mathrm{HCO_{3}^{-}}$-$\mathrm{CO_{3}^{2-}}$
is approximately one order of magnitude higher),
which in turn, in the presence of calcium ions, can lead to a supersaturation
of calcium carbonate in the solution,
thereby promoting the precipitation of calcium carbonate.
The dissociation coefficients of $\mathrm{H_2CO_{3}}$ and $\mathrm{NH_{3}}$ are calculated
using relations given by \citet{Millero2007} and \citet{Bell2008} respectively.

Another key parameter is the distribution of biomass within the porous medium,
since the presence of ureolytic microbes is a prerequisite for MICP.
As the distribution of biomass is not exactly known for the MICP experiments,
it is not possible to validate the equations
governing the distribution of biomass directly by comparison with biofilm measurements.
However, the ureolytic activity of the bacteria,
which is monitored by measuring the product $\mathrm{NH_{4}^{+}}$
in recent column experiments, allows for an indirect evaluation of the
ability of the model to predict the biomass distribution.




[@Ebigbo2012]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR011714 "Darcy-scale modeling of microbially induced carbonate mineral precipitation in sand columns"
[@Fidaleo2003]: http://pierre.fkit.hr/hdki/cabeq/pdf/17_4_2003/Cabeq%202003-04_3.pdf "Kinetic study of enzymatic urea hydrolysis in the {pH} range 4-9"
[@Lauchnor2015]: http://dx.doi.org/10.1111/jam.12804 "Whole cell kinetics of ureolysis by Sporosarcina pasteurii"

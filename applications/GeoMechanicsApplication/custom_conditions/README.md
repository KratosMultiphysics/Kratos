# Microclimate Boundary Condition
A microclimate boundary condition refers to the specific set of environmental factors and conditions that influence the climate within a small, localized area, known as a microclimate. Unlike the broader regional or global climate, microclimates are unique climatic conditions that can vary significantly over a relatively small distance.

Microclimate boundary conditions are determined by factors such as:

- Air temperature
- Solar radiation
- Humidity
- Precipitation
- Wind speed

## Radiation
The net radiation $R_n$ $\mathrm{[W/m^2]}$ is calculated by:

$$ R_n = R^s + R^{l,ab} - R^{l,em} $$

where

- $R^s$ = absorbed solar radiation $\mathrm{[W/m^2]}$
- $R^{l,ab}$ = absorbed atmospheric radiation $\mathrm{[W/m^2]}$
- $R^{l,em}$ = emitted black/grey body radiation from the soil $\mathrm{[W/m^2]}$

The incoming short wave radiation is defined by

$$ R^s = \left(1 - \alpha \right) R_{g} $$

where

- $\alpha$ = Albedo cover or vegetation coefficient  $\mathrm{[-]}$
- $R_{g}$ = Short wave radiation  $\mathrm{[W/m^2]}$ 

$R_g$ is a user defined parameter which is usually provided as a "Time-Radiation" table.

The absorbed long wave radiation is:

$$ R^{l,ab}=  \epsilon \sigma \left( T_{at} +273.15 \right)^4 $$

where

- $\epsilon$ = effective emissivity  $\mathrm{[-]}$
- $\sigma$ = Stefan Boltzmann's constant  $\mathrm{[5.67 \cdot 10^{-8}  W/m^{2 \circ}C^4]}$
- $T_{at}$ = atmospheric temperature $\mathrm{[^{\circ}C]}$

And the emitted long wave radiation,

$$ R^{l,em} = \sigma \left( T_{ss} +273.15 \right)^4 $$

where $T_{ss}$ $\mathrm{[^{\circ}C]}$ is land surface temperature.

## Surface Heat Storage
The surface heat storage term $Q_s$ $\mathrm{[W/m^2]}$ is calculated by,

$$ Q_s =  a_1 R_n(t+\Delta t) + a_2 \frac{R_n(t+\Delta t) - R_n(t)}{\Delta t} + a_3 $$

where $a_1$ $\mathrm{[-]}$, $a_2$ $\mathrm{[s]}$ and $a_3$ $\mathrm{[W/m^2]}$ are user defined variables.


## Evaporation Flux
The evaporation flux $E_a$ $\mathrm{[m/s]}$ is calculated by,

$$ E_a = \frac{L_v}{\lambda \rho_w} $$

where

- $L_v$ = latent heat flux
- $\lambda$ = latent heat of vaporization $\mathrm{[2.45 \space MJ/kg]}$
- $\rho_w$  = density of water $\mathrm{[kg/m^3]}$

The latent heat flux is:

$$ L_v  = \frac{e^\prime_{at} \left(R_n + Q_f - Q_s \right) + C_a \rho_a \left(e_{at}^s - e_{at}^a\right) / r_a}{e^\prime_{at} + \gamma (1+r_{s} / r_a)} $$

where

- $r_a$ = atomospheric resistance $\mathrm{[s/m]}$
- $r_s$ = surface resistance constant $\mathrm{[30 \space s/m]}$ 
- $\gamma$ = psychometric constant $\mathrm{[0.67 \space hPa/ ^{\circ}C]}$
- $Q_f$ = build enviromental radiation $\mathrm{[W/m^2]}$ (user defined parameter)
- $C_a$ = specific heat of moist air $\mathrm{[J/kg ^{\circ}C]}$
- $\rho_a$ = air density $\mathrm{[kg/m^3]}$ 
- $e^s$ = saturated vapour pressure $\mathrm{[hPa]}$ 
- $e^\prime$ = slope of the saturation vapor curve $\mathrm{[hPa/ ^{\circ}C]}$
- $e^a$ = actual vapour pressure $\mathrm{[hPa]}$

$$ e^s = 6.11 \exp \left( \frac{17.27 \space T}{T+237.3} \right) $$

$$ e^\prime=  \frac{4098 \space e^s}{\left(T + 237.3\right)^2} $$   

$$ e^a = \frac{\mathrm{RH}}{100} e^s $$

$\mathrm{RH}$ [%] is the relative humidity which is a user defined parameter.

The atomospheric resistance,

$$ r_a = \frac{1}{0.007+0.0056 u} $$

where $u$ $\mathrm{[m/s]}$ is the wind speed which is a user defined parameter. $t$ $\mathrm{[s]}$ and $\Delta t$ $\mathrm{[s]}$ are time and time step, respectively.

## Roughness Layer Temperature
The energy balance equation for the roughness layer reads as,

$$ h_{rl} \frac{\partial T_{rl}}{\partial t} = \frac{T_{ss}-T_{rl}}{r_g} + u f_h a_d^2 \left(T_{at}-T_{rl}\right)  $$

After discretization,

$$ h_{rl} \frac{T_{rl}(t+\Delta t) - T_{rl}(t)}{\Delta t} = \frac{T_{ss}(t)-T_{rl}(t+\Delta t)}{r_g} + u f_h a_d^2 \left[T_{at}(t+\Delta t)-T_{rl}(t+\Delta t) \right] $$

Then the temperature in the roughness layer $T_{rl}$ reads as,

$$ T_{rl}(t+\Delta t) = \frac{r_g \space h_{rl} \space T_{rl}^(t) + \Delta t \space T_{ss}(t) + r_g \space \Delta t \space u \space f_h \space d_d^2 \space T_{at}(t+\Delta t)}{r_g \space h_{rl} + \Delta t + r_g \space \Delta t \space u \space f_h \space a_d^2} $$

where

- $f_h$ = surface roughness factor $\mathrm{[-]}$
- $a_d$ = friction drag coefficient $\mathrm{[-]}$
- $r_g$ = roughness layer resistance $\mathrm{[s/m]}$
- $T_{at}$ = atmospheric temperature $\mathrm{[C]}$

The friction drag coefficient $a_d$ is:

$$ a_d = \frac{\kappa}{\ln\left(z_m / z_0\right)} $$

where

- $\kappa$ = friction drag coefficient $\mathrm{[-]}$
- $z_m$ = measurment height $\mathrm{[10 m]}$
- $z_0$ = roughness length $\mathrm{[1 m]}$

Surface roughness factor $f_h$ for unstable weather conditions where $T_{rl} \geq T_{at}$

$$ f_h = 1 - \frac{15 \space r_{i}}{1+75 \space a_d^2 \sqrt{z_m/z_0} \sqrt{|r_{i}|}} $$

and for stable weather conditions where $T_{rl} < T_{at}$

$$ f_h = \frac{}{1 + 15 \space r_i \sqrt{1+5 \space r_i}} $$

where $r_i$ $\mathrm{[-]}$ is Richardson bulk modulus and is defines as,

$$ r_i = \frac{2 g z_m}{T_{at}+T_{rl}+ 546.3} \frac{T_{at}-T_{rl}}{u^2} $$

$g$ is gravitation constant $\mathrm{[9.81 m/s^2]}$

## Actual Precipitation and Evaporation
For the computation of the actual precipitation and actual evaporation the potential and actual storage is introduced. The influx into the soil follows from the actual storage of water.
The potential storage $S_p^{j+1}$ $\mathrm{[mm]}$ at time step $j+1$ follows from the actual storage $S_a^j$ $\mathrm{[mm]}$ at the previous time step $j$, the time step size $\Delta t$ $\mathrm{[s]}$ and the difference in potential precipitation $P_p^{j+1}$ $\mathrm{[mm/s]}$ and the potential evaporation $E_p^{j+1}$ $\mathrm{[mm/s]}$ during the new time step as:

$$ S_p^{j+1} = S_a^j + \Delta t \left( P_p^{j+1} - E_p^{j+1} \right) $$

If the potential storage is larger than the maximum storage $S_{max}$ $\mathrm{[mm]}$ then the actual evaporation $E_a^{j+1}$ $\mathrm{[mm/s]}$ and the actual precipitation $P_a^{j+1}$ $\mathrm{[mm/s]}$ follow from:

$$ E_a^{j+1} = E_p^{j+1}$$

$$ P_a^{j+1} = \left( S_{max} - S_p^j \right) / \Delta t + E_a^{j+1} $$

If the potential storage is smaller than the minimum storage $S_{min}$ $\mathrm{[mm]}$ then the actual values read:

$$ P_a^{j+1} = P_p^{j+1}$$

$$ E_a^{j+1} = \left( S_p^j - S_{min} \right) / \Delta t + P_a^{j+1} $$

If the potential storage does not exceed the storage limits then the actual fluxes match the potential fluxes:

$$ P_a^{j+1} = P_p^{j+1}$$

$$ E_a^{j+1} = E_p^{j+1}$$

The maximum storage capacity $S_{max}$ is a user defined parameter. The potential precipitation $P_p$ is also given by the user, usually in the form of time-precipitation table in the MPDA file. 
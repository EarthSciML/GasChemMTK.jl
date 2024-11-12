export FastJX

# Effective wavelength in 18 bins covering 177–850 nm
const WL = SA_F32[187, 191, 193, 196, 202, 208, 211, 214, 261, 267, 277, 295, 303, 310, 316, 333, 380, 574]

function interp2_func(x1, x2, T1, T2)
    if x1 ≈ x2
        return (T) -> x1
    end
    function interp2(T)
        T = clamp(T, T1, T2)
        x1 + (T - T1) * (x2 - x1) / (T2 - T1)
    end
end

function interp3_func(x1, x2, x3, T1, T2, T3)
    if x1 ≈ x2 && x2 ≈ x3
        return (T) -> x1
    end
    function interp3(T)
        T = clamp(T, T1, T3)
        ifelse(T < T2,
            x1 + (T - T1) * (x2 - x1) / (T2 - T1),
            x2 + (T - T2) * (x3 - x2) / (T3 - T2),
        )
    end
end

function interp_func(temperatures, cross_sections)
    if length(temperatures) == 2
        return interp2_func(cross_sections[1], cross_sections[2], temperatures[1], temperatures[2])
    elseif length(temperatures) == 3
        return interp3_func(cross_sections[1], cross_sections[2], cross_sections[3],
            temperatures[1], temperatures[2], temperatures[3])
    end
    LinearInterpolation(temperatures, cross_sections, extrapolation_bc=Flat())
end

"""
Create a vector of interpolators to interpolate the cross sections σ (TODO: What are the units?) for different wavelengths (in nm) and temperatures (in K).

We use use linear interpolation with flat extrapolation.
"""
function create_fjx_interp(temperatures::Vector{Float32}, cross_sections::Vector{SVector{18,Float32}})
    [interp_func(temperatures, [x[i] for x ∈ cross_sections]) for i ∈ 1:length(WL)]
    #[linear_interpolation(temperatures, [x[i] for x ∈ cross_sections], extrapolation_bc=Flat()) for i ∈ 1:length(WL)]
end

#   Cross sections and quantum yield from GEOS-CHEM "FJX_spec.dat" for photo-reactive species included in SuperFast:

# Cross sections σ for O3(1D) for different wavelengths(18 bins), temperatures
# O3 -> O2 + O(1D) (superfast: O3 -> 2OH  including O(1D) + H2O -> 2OH)
const ϕ_o31D_jx = 1.0f0
const σ_o31D_interp = create_fjx_interp([200.0f0, 260.0f0, 320.0f0], [
    SA_F32[4.842, 4.922, 5.071, 5.228, 6.040, 6.803, 7.190, 7.549, 9.000, 8.989, 8.929, 9.000, 8.901, 4.130, 8.985*0.1, 6.782*0.1, 0.000, 0.000] * 0.1f0,
    SA_F32[4.843, 4.922, 5.072, 5.229, 6.040, 6.802, 7.189, 7.549, 9.000, 8.989, 8.929, 9.000, 8.916, 4.656, 1.417, 6.995*0.1, 0.000, 0.000] * 0.1f0,
    SA_F32[4.843, 4.922, 5.072, 5.229, 6.040, 6.805, 7.189, 7.550, 9.000, 8.989, 8.929, 9.000, 8.967, 5.852, 2.919, 7.943, 0.000, 0.000] * 0.1f0,
])

# Cross sections σ for H2O2 for different wavelengths(18 bins), temperatures
# H2O2 -> OH + OH (same with superfast)
const ϕ_H2O2_jx = 1.0f0
const σ_H2O2_interp = create_fjx_interp([200.0f0, 300.0f0], [
    SA_F32[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.346*0.1, 4.855*0.1, 3.423*0.1, 8.407*0.01, 5.029*0.01, 3.308*0.01, 2.221*0.01, 8.598*0.001, 1.807*0.0001, 0] * 10.0f0^-19.0f0,
    SA_F32[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.465*0.1, 4.966*0.1, 3.524*0.1, 9.354*0.01, 5.763*0.01, 3.911*0.01, 2.718*0.01, 1.138*0.01, 2.419*0.0001, 0] * 10.0f0^-19.0f0,
])

# Cross sections σ for CH2Oa for different wavelengths(18 bins), temperatures
# CH2O -> H + HO2 + CO (superfast: CH2O -> CO + 2HO2 including H + O2 -> HO2)
const ϕ_CH2Oa_jx = 1.0f0
const σ_CH2Oa_interp = create_fjx_interp([223.0f0, 298.0f0], [
    SA_F32[0.0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.143*0.1, 1.021, 1.269, 2.323, 2.498, 1.133, 2.183, 4.746*0.1, 0.000, 0.000] * 10.0f0^-20.0f0,
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.147*0.1, 1.018, 1.266, 2.315, 2.497, 1.131, 2.189, 4.751*0.1, 0.000, 0.000] * 10.0f0^-20.0f0,
])


# Cross sections σ for CH2Ob for different wavelengths(18 bins), temperatures
# CH2O -> CO + H2 （same with superfast）
const ϕ_CH2Ob_jx = 1.0f0
const σ_CH2Ob_interp = create_fjx_interp([223.0f0, 298.0f0], [
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.642*10^-21, 5.787*10^-21, 5.316*10^-21, 8.181*10^-21, 7.917*10^-21, 4.011*10^-21, 1.081*10^-20, 1.082^-20, 2.088^-22, 0.000],
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.649*10^-21, 5.768*10^-21, 5.305*10^-21, 8.154*10^-21, 7.914*10^-21, 4.002*10^-21, 1.085*10^-20, 1.085*10^-20, 2.081^-22, 0.000],
])

# Cross sections σ for CH3OOH for different wavelengths(18 bins), temperatures
# CH3OOH -> OH + HO2 + CH2O (same with superfast)
const ϕ_CH3OOH_jx = 1.0f0
const σ_CH3OOH = SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 3.120, 2.882, 2.250, 2.716*0.1, 2.740*0.1, 2.143*0.1, 5.624*0.01, 3.520*0.01, 2.403*0.01, 1.697*0.01, 7.230*0.001, 6.973*0.0001, 0.000] * 10.0f0^-19.0f0 # 298K
const σ_CH3OOH_interp = [(T) -> σ_CH3OOH[i] for i in 1:18]

# Cross sections σ for NO2 for different wavelengths(18 bins), temperatures
# NO2 -> NO + O (superfast: NO2 -> NO + O3  including O + O2 -> O3)
const ϕ_NO2_jx = 1.0f0
const σ_NO2_interp = create_fjx_interp([200.0f0, 294.0f0], [
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.835*0.1, 4.693*0.1, 7.705*0.1, 1.078, 1.470, 1.832, 2.181, 3.138, 4.321, 1.386*0.001] * 10.0f0^-19.0f0,
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.313*0.1, 4.694*0.1, 7.553*0.1, 1.063, 1.477, 1.869, 2.295, 3.448, 4.643, 4.345*0.001] * 10.0f0^-19.0f0,
])


const actinic_flux = SA_F32[1.391E+12, 1.627E+12, 1.664E+12, 9.278E+11, 7.842E+12, 4.680E+12, 9.918E+12, 1.219E+13, 6.364E+14, 4.049E+14, 3.150E+14, 5.889E+14, 7.678E+14, 5.045E+14, 8.902E+14, 3.853E+15, 1.547E+16, 2.131E+17]

"""
    cos_solar_zenith_angle(lat, t, long)
This function is to compute the cosine of the solar zenith angle, given the unixtime, latitude and longitude
The input variables: lat=latitude(°), long=longitude(°), t=unixtime(s)
    the cosine of the solar zenith angle (SZA) is given by:                                                                            .
           cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)

           where LAT = the latitude angle,
                 DEC = the solar declination angle,
                 AHR = the hour angle, all in radians.  All in radians
"""
function cos_solar_zenith_angle(t, lat, long)
    ut = Dates.unix2datetime(t)
    DOY = dayofyear(ut)
    hours = Dates.hour(ut) + Dates.minute(ut) / 60 + Dates.second(ut) / 3600
    y = Dates.year(ut)
    yearday_temp = 2 * pi * (DOY - 1 + (hours - 12) / 24)
    γ = ifelse(mod(y, 4) == 0, # the fraction year in radians
        yearday_temp / 366,
        yearday_temp / 365,
    )
    Eot = 229.18 * (0.000075 + 0.001868 * cos(γ) - 0.032077 * sin(γ) - 0.014615 * cos(γ * 2) - 0.040849 * sin(γ * 2))

    timezone = floor(long / 15)
    time_offset = Eot + 4 * (long - 15 * timezone) # in minutes
    dt = floor(long / 15) # in hours
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    tst = Dates.hour(ut_local) + Dates.minute(ut_local) / 60 +
          Dates.second(ut_local) / 3600 + time_offset / 60
    AHR = deg2rad(15) * (tst - 12) # in radians

    LAT = abs(lat * pi / 180) #lat>0, northern hemisphere; lat<0, southern hemisphere
    DEC = asin(sin(deg2rad(-23.44)) * cos(deg2rad(360 / 365.24 * (DOY + 10) + 360 / pi * 0.0167 * sin(deg2rad(360 / 365.24 * (DOY - 2))))))
    CSZA = sin(LAT) * sin(DEC) + cos(LAT) * cos(DEC) * cos(AHR)
    return CSZA
end

@register_symbolic cos_solar_zenith_angle(t, lat, long)

# Dummy function for unit validation. Basically ModelingToolkit
# will call the function with a DynamicQuantities.Quantity or an integer to
# get information about the type and units of the output.
cos_solar_zenith_angle(t::DynamicQuantities.Quantity, lat, long) = 1.0


# calculate actinic flux at the given cosine of the solar zenith angle `csa` and
# maximium actinic flux `max_actinic_flux`
function calc_flux(csa, max_actinic_flux)
    max(zero(max_actinic_flux), max_actinic_flux * csa)
end

# calculate all actinic fluxes
function calc_fluxes(csa)
    [calc_flux(csa, actinic_flux[i]) for i in 1:18]
end

# Symbolic equations for actinic flux
function flux_eqs(csa)
    flux_vals = calc_fluxes(csa)
    flux_vars = []
    @constants c_flux = 1.0 #[unit = u"s^-1", description = "Constant actinic flux (for unit conversion)"]
    for i in 1:18
        wl = WL[i]
        n = Symbol("F_", Int(round(wl)))
        v = @variables $n(t) #[unit = u"s^-1", description = "Actinic flux at $wl nm"]
        push!(flux_vars, only(v))
    end
    flux_vars, (flux_vars .~ flux_vals .* c_flux)
end


"""
Get mean photolysis rates at different times
"""
function j_mean(σ_interp, ϕ, Temperature, fluxes)
    j = zero(Temperature)
    for i in 1:18
        j += fluxes[i] * σ_interp[i](Temperature) * ϕ
    end
    j
end

j_mean_H2O2(T, fluxes) = j_mean(σ_H2O2_interp, ϕ_H2O2_jx, T, fluxes)
j_mean_CH2Oa(T, fluxes) = j_mean(σ_CH2Oa_interp, ϕ_CH2Oa_jx, T, fluxes)
j_mean_CH2Ob(T, fluxes) = j_mean(σ_CH2Ob_interp, ϕ_CH2Ob_jx, T, fluxes)
j_mean_CH3OOH(T, fluxes) = j_mean(σ_CH3OOH_interp, ϕ_CH3OOH_jx, T, fluxes)
j_mean_NO2(T, fluxes) = j_mean(σ_NO2_interp, ϕ_NO2_jx, T, fluxes)
j_mean_o31D(T, fluxes) = j_mean(σ_o31D_interp, ϕ_o31D_jx, T, fluxes)

"""
    adjust_j_o31D(T, P, H2O)

Adjust the photolysis rate of O3 -> O2 + O(1D) to represent the effective rate for O3 -> 2OH.
This adjustment is based on the fraction of O(1D) that reacts with H2O to produce 2 OH.
"""
function adjust_j_o31D(T, P, H2O)
    @constants(
        T_unit = 1, [unit = u"K", description = "inverse unit of Temperature"],
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
        ppb_inv = 1, [unit = u"ppb^-1"],
    )
    num_density_unitless = A*P/(R*T)*num_density_unit_inv

    # Define species concentrations value in unit of molec/cm3
    C_H2O = H2O*ppb_inv*1e-9*num_density_unitless # convert value of H2O concentration in unit of ppb to unit of molec/cm3, but here is unitless
    C_O2 = 0.2095 * num_density_unitless
    C_N2 = 0.7808 * num_density_unitless
    C_H2 = 0.5e-6 * num_density_unitless

    # Define rate constants for reactions involving O(1D)
    RO1DplH2O = 1.63e-10 * exp(60.0*T_unit / T) * C_H2O
    RO1DplH2  = 1.2e-10 * C_H2
    RO1DplN2  = 2.15e-11 * exp(110.0*T_unit / T) * C_N2
    RO1DplO2  = 3.30e-11 * exp(55.0*T_unit / T) * C_O2

    # Total rate constant for O(1D)
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2

    # Prevent division by zero
    return RO1DplH2O / RO1D
end

struct FastJXCoupler
    sys
end

"""
Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme
(Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

# Example
Build Fast-JX model:
``` julia
    fj = FastJX()
```
"""
function FastJX(; name=:FastJX)
    @constants T_unit = 1.0 [unit = u"K", description = "Unit temperature (for unit conversion)"]
    @parameters T = 298.0 [unit = u"K", description = "Temperature"]
    @parameters lat = 40.0 [description = "Latitude (Degrees)"]
    @parameters long = -97.0 [description = "Longitude (Degrees)"]
    @parameters P = 101325 [unit = u"Pa", description = "Pressure"]
    @parameters H2O = 450 [unit = u"ppb"]

    @variables j_h2o2(t) = 1.0097 * 10.0^-5 #[unit = u"s^-1"]
    @variables j_CH2Oa(t) = 0.00014 #[unit = u"s^-1"]
    @variables j_CH2Ob(t) = 0.00014 #[unit = u"s^-1"]
    @variables j_o31D(t) = 4.0 * 10.0^-3 #[unit = u"s^-1"]
    @variables j_o32OH(t) = 2.27e-4 #[unit = u"s^-1"]
    @variables j_CH3OOH(t) = 8.9573 * 10.0^-6 #[unit = u"s^-1"]
    @variables j_NO2(t) = 0.0149 #[unit = u"s^-1"]
    @variables cosSZA(t) [description = "Cosine of the solar zenith angle"]
    # TODO(JL): What's difference between the two photolysis reactions of CH2O, do we really need both?
    # (@variables j_CH2Ob(t) = 0.00014 [unit = u"s^-1"]) (j_CH2Ob ~ mean_J_CH2Ob(t,lat,T)*j_unit)

    flux_vars, fluxeqs = flux_eqs(cosSZA)
    eqs = [
        cosSZA ~ cos_solar_zenith_angle(t, lat, long);
        fluxeqs;
        j_h2o2 ~ j_mean_H2O2(T/T_unit, flux_vars);
        j_CH2Oa ~ j_mean_CH2Oa(T/T_unit, flux_vars);
        j_CH2Ob ~ j_mean_CH2Ob(T/T_unit, flux_vars);
        j_o31D ~ j_mean_o31D(T/T_unit, flux_vars)*1e-17;
        j_o32OH ~ j_o31D*adjust_j_o31D(T, P, H2O);
        j_CH3OOH ~ j_mean_CH3OOH(T/T_unit, flux_vars);
        j_NO2 ~ j_mean_NO2(T/T_unit, flux_vars)
    ]

    ODESystem(eqs, t, [j_h2o2, j_CH2Oa, j_CH2Ob, j_o32OH, j_o31D, j_CH3OOH, j_NO2, cosSZA, flux_vars...],
        [lat, long, T, P, H2O]; name=name, metadata=Dict(:coupletype => FastJXCoupler))
end
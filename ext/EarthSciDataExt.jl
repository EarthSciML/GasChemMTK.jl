module EarthSciDataExt
using GasChem, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities
@register_unit ppb 1

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, e::EarthSciData.NEI2016MonthlyEmisCoupler)
    c, e = c.sys, e.sys

    @constants(
        MW_NO2 = 46.0055*10^-3, [unit = u"kg/mol", description="NO2 molar mass"],
        MW_NO = 30.01*10^-3, [unit = u"kg/mol", description="NO molar mass"],
        MW_FORM = 30.0260*10^-3, [unit = u"kg/mol", description="Formaldehyde molar mass"],
        MW_CH4 = 16.0425*10^-3, [unit = u"kg/mol", description="Methane molar mass"],
        MW_CO = 28.0101*10^-3, [unit = u"kg/mol", description="Carbon monoxide molar mass"],
        MW_SO2 = 64.0638*10^-3, [unit = u"kg/mol", description="Sulfur dioxide molar mass"],
        MW_ISOP = 68.12*10^-3, [unit = u"kg/mol", description="Isoprene molar mass"],

        nmolpermol = 1e9, [unit = u"ppb", description="noml/mol, Conversion factor from mol to nmol"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
    )

    # Emissions are in units of "kg/m3/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # To do this we need to convert kg of emissions to nmol of emissions,
    # and we need to convert m3 of air to mol of air.
    # nmol_emissions = kg_emissions * gperkg / MW_emission * nmolpermol = kg / kg/mol * nmol/mol = nmol
    # mol_air = m3_air / R / T * P = m3 / (m3*Pa/mol/K) / K * Pa = mol
    # So, the overall conversion is:
    # nmol_emissions / mol_air = (kg_emissions / MW_emission * nmolpermol) / (m3_air / R / T * P)
    uconv =  nmolpermol * R * c.T / c.P # Conversion factor with MW factored out.
    operator_compose(convert(ODESystem, c), e, Dict(
        c.NO2 => e.NO2 => uconv / MW_NO2,
        c.NO => e.NO => uconv / MW_NO,
        c.CH2O => e.FORM => uconv / MW_FORM,
        c.CO => e.CO => uconv / MW_CO,
        #c.SO2 => e.SULF => uconv / MW_SO2, #TODO Need to find a way to make c.SO2 => e.SULF + e.SO2 => uconv / MW_SO2
        c.SO2 => e.SO2 => uconv / MW_SO2,
        c.ISOP => e.ISOP => uconv / MW_ISOP,
    ))
end

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, g::EarthSciData.GEOSFPCoupler)
    c, g = c.sys, g.sys

    @constants(
        MW_H2O = 18.015*10^-3, [unit = u"kg/mol", description="Water molar mass"],
        MW_air = 28.97*10^-3, [unit = u"kg/mol", description = "Air molar mass"],
        ppbpersecond = 1e9, [unit = u"ppb/s"],
    )
    uconv = MW_air/MW_H2O*ppbpersecond
    operator_compose(convert(ODESystem, c), g, Dict(
        c.H2O => g.I3₊QV => uconv,
    ))

    c = param_to_var(c, :T, :P)
    ConnectorSystem([
            c.T ~ g.I3₊T,
            c.P ~ g.P,
        ], c, g)
end

function EarthSciMLBase.couple2(f::GasChem.FastJXCoupler, g::EarthSciData.GEOSFPCoupler)
    f, g = f.sys, g.sys

    f = param_to_var(f, :T, :lat, :long)
    ConnectorSystem([
            f.T ~ g.I3₊T,
            f.lat ~ g.lat,
            f.long ~ g.lon,
        ], f, g)
end

end

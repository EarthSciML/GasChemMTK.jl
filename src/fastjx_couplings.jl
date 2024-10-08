
function EarthSciMLBase.couple2(c::GEOSChemGasPhaseCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(c, :j_3, :j_9, :j_11, :j_7, :j_8, :j_10)
    @constants uconv = 1 [unit = u"s"]
    @constants c_fixme1 = 10^(-21) [unit = u"s", description="Suspicious constant"] # FIXME: Suspicious constant
    ConnectorSystem([
            c.j_9 ~ uconv * p.j_h2o2
            c.j_7 ~ uconv * p.j_CH2Oa
            c.j_8 ~ uconv * p.j_CH2Ob
            c.j_10 ~ uconv * p.j_CH3OOH
            c.j_11 ~ uconv * p.j_NO2
            c.j_3 ~ c_fixme1 * p.j_o31D
        ], c, p)
end

function EarthSciMLBase.couple2(c::SuperFastCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(convert(ODESystem, c), :jO31D, :jH2O2, :jNO2, :jCH2Oa, :jCH2Ob, :jCH3OOH)
    ConnectorSystem([
            c.jH2O2 ~ p.j_h2o2
            c.jCH2Oa ~ p.j_CH2Oa
            c.jCH2Ob ~ p.j_CH2Ob
            c.jCH3OOH ~ p.j_CH3OOH
            c.jNO2 ~ p.j_NO2
            c.jO31D ~ p.j_o31D
        ], c, p)
end

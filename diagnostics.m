function diagnostics(sat)
    clc
    answ = input("Show model data (y/n)\n\n","s");

    if answ == "y" || answ == "Y"
        answ = input("\n 1.) All \n 2.) Power \n 3.) Propulsion \n 4.) Thermals \n 5.) Communications \n 6.) Orbital Mechanics \n 7.) Mass \n 8.) Exit \n\n","s");
        while answ == "1" || answ == "2" || answ == "3" || answ == "4" || answ == "5" || answ == "6" || answ == "7"
            
            switch answ

            case "1"
                clc
                power(sat)
                propulsion(sat)
                thermals(sat)
                comms(sat)
                orbital_mechanics(sat)
                mass(sat)
            case "2"
                clc
                power(sat)
            case "3"
                clc
                propulsion(sat)  
            case "4"
                clc
                thermals(sat)
            case "5"
                clc
                comms(sat)
            case "6"
                clc
                orbital_mechanics(sat)
            case "7"
                clc
                mass(sat)
            otherwise
                fprintf("No valid option selected, exiting... \n")
                return
            end
            
            fprintf("******** Repeating ********\n")
            
            answ = input("\n 1.) All \n 2.) Power \n 3.) Propulsion \n 4.) Thermals \n 5.) Communications \n 6.) Orbital Mechanics \n 7.) Mass \n 8.) Exit \n\n","s");
               
        end
    else
        fprintf("No valid option selected, exiting... \n")
        return
    end
end

function orbital_mechanics(sat)

    fprintf("                 ORBITAL MECHANICS                 \n\n")
    fprintf("Semimajor Axis:                "+sat.orbit.sma+" m\n")
    fprintf("Eccentricity:                  "+sat.orbit.e+"\n")
    fprintf("Inclination:                   "+sat.orbit.i+" deg\n")
    fprintf("RAAN:                          "+sat.orbit.raan+" deg\n")
    fprintf("Argument of Perigee:           "+sat.orbit.ap+" deg\n")
    fprintf("Minimum Altitude:              "+sat.orbit.alt_min+" m\n")
    fprintf("Maximum Altitude:              "+sat.orbit.alt_max+" m\n")
    fprintf("Average Altitude:              "+sat.orbit.pertubations.average_altitude*1000+" m\n")
    fprintf("Average Inclination:           "+sat.orbit.pertubations.average_inclination +" deg\n")
    fprintf("\n")
    return
end

function mass(sat)
    fprintf("                 MASS BUDGET                 \n\n")
    fprintf("Wet/Launch Mass:               "+sat.massbudget.wet+" kg\n")
    fprintf("Dry Mass:                      "+sat.massbudget.dry+" kg\n")
    fprintf("Propellant Mass:               "+sat.massbudget.propellant+" kg\n")
    fprintf("Excess Mass:                   "+sat.massbudget.excess+" kg\n")
    fprintf("Maximum:                       "+sat.massbudget.max_subsystem+" kg\n")
    fprintf("Power:                         "+sat.massbudget.power+" kg\n")
    fprintf("Thermals:                      "+sat.massbudget.thermal+" kg\n")
    fprintf("Communications:                "+sat.massbudget.ttc+" kg\n")
    fprintf("Propulsion:                    "+sat.massbudget.propulsion+"  kg\n")
    fprintf("ADCS:                          "+sat.massbudget.adcs+" kg\n")
    fprintf("Structures:                    "+sat.massbudget.structure+" kg\n")
    fprintf("Computing:                     "+sat.massbudget.computing+" kg\n")
    fprintf("Other:                         "+sat.massbudget.other+" kg\n")
    fprintf("Payload Margin:                "+sat.massbudget.payload_margin+" kg\n")
    fprintf("Structure Margin:              "+sat.massbudget.structure_margin+" kg\n")
    fprintf("\n")
    return
end

function thermals(sat)
    fprintf("                 THERMALS                 \n\n")
    fprintf("Absorbed Solar Heat:           "+sat.thermals.absorbed_heat_solar+" W\n")
    fprintf("Reflected Solar Heat:          "+sat.thermals.absorbed_heat_reflected+" W\n")
    fprintf("Absorbed Infrared Heat:        "+sat.thermals.absorbed_heat_infrared+" W\n")
    fprintf("Absorbed Heat:                 "+sat.thermals.absorbed_heat+" W\n")
    fprintf("Radiated Heat:                 "+sat.thermals.max_radiated_heat+" W\n")
    fprintf("Equilibrium Temperature:       "+sat.thermals.max_equilibrium_temperature+" K\n")
    
    fprintf("\n")
    return
end

function power(sat)
    fprintf("                 POWER BUDGET                 \n\n")
    fprintf("Duty Cycles\n\n")
    fprintf("Comms Duty Cycle Eclipse:      "+sat.powerbudget.ttc_cyc_eclipse+" \n")
    fprintf("Comms Duty Cycle Light:        "+sat.powerbudget.ttc_cyc_light+" \n")
    fprintf("Comms Duty Cycle Eclipse:      "+sat.powerbudget.ttc_cyc_eclipse+" \n")
    fprintf("ADCS Duty Cycle Light:         "+sat.powerbudget.acs_cyc_light+" \n")
    fprintf("ADCS Duty Cycle Eclipse:       "+sat.powerbudget.acs_cyc_eclipse+" \n")
    fprintf("Payload Duty Cycle Light:      "+sat.powerbudget.payload_cyc_light+" \n")
    fprintf("Payload Duty Cycle Eclipse:    "+sat.powerbudget.payload_cyc_eclipse+" \n")
    fprintf("Thermals Duty Cycle Light:     "+sat.powerbudget.thermals_cyc_light+" \n")
    fprintf("Thermals Duty Cycle Eclipse:   "+sat.powerbudget.thermals_cyc_eclipse+" \n")
    fprintf("Processing Duty Cycle Light:   "+sat.powerbudget.processing_cyc_light+" \n")
    fprintf("Processing Duty Cycle Eclipse: "+sat.powerbudget.processing_cyc_eclipse+" \n")
    fprintf("Propulsion Duty Cycle Light:   "+sat.powerbudget.prop_cyc_light+" \n")
    fprintf("Propulsion Duty Cycle Eclipse: "+sat.powerbudget.prop_cyc_eclipse+" \n")
    fprintf("\n")
    fprintf("Payload Power:                 "+sat.powerbudget.payload+" W\n")
    fprintf("Nominal Power:                 "+sat.powerbudget.nominal+" W\n")
    fprintf("Structure Power:               "+sat.powerbudget.structure+" W\n")
    fprintf("Thermals Power:                "+sat.powerbudget.thermal+" W\n")
    fprintf("Power Power:                   "+sat.powerbudget.power+" W\n")
    fprintf("Comms Power:                   "+sat.powerbudget.ttc+" W\n")
    fprintf("Processing Power:              "+sat.powerbudget.processing+" W\n")
    fprintf("ADCSPower:                     "+sat.powerbudget.adcs+" W\n")
    fprintf("Propulsion Power:              "+sat.powerbudget.propulsion+" W\n")
    fprintf("\n")
    fprintf("Max Power Generated:           "+sat.powerbudget.generated_max+" W\n")
    fprintf("Max EOL Power Generated:       "+sat.powerbudget.generated_eol_max+" W\n")
    fprintf("Max EOL Efficiency:            "+sat.powerbudget.solar_eol_eff+" W\n")
    fprintf("Max Power Generated:           "+sat.powerbudget.generated_max+" W\n")
    fprintf("\n")
    return
end

function comms(sat)
    fprintf("                 COMMUNICATIONS                 \n")
    fprintf("Frequency:                     "+sat.comms.frequency+" Hz\n")
    fprintf("Wavelength:                    "+sat.comms.wavelength+" m\n")
    fprintf("Radius:                        "+sat.comms.altitude+" m\n")
    fprintf("Power:                         "+sat.comms.power+" W\n")
    fprintf("TX Area:                       "+sat.comms.tx_effective_area+" m^2\n")
    fprintf("TX Efficiency:                 "+sat.comms.tx_efficiency+" \n")
    fprintf("TX Gain:                       "+sat.comms.tx_gain+"\n")
    fprintf("EIRP:                          "+sat.comms.eirp+" W\n")
    fprintf("RX Area:                       "+sat.comms.rx_effective_area+" m^2\n")
    fprintf("RX Efficiency:                 "+sat.comms.rx_efficiency+" \n")
    fprintf("RX Gain:                       "+sat.comms.rx_gain+"\n")
    fprintf("RX Temperature:                "+sat.comms.rx_temp+" K\n")
    fprintf("Bandwdth:                      "+sat.comms.bandwidth+" Hz\n")
    fprintf("Noise Power:                   "+sat.comms.noise_power+" W\n")
    fprintf("Signal to Noise Ratio:         "+sat.comms.sn+"\n")
    fprintf("Free Space Loss:               "+sat.comms.free_space_loss+"\n")

    return
end

function propulsion(sat)
    fprintf("                 PROPULSION                 \n")
    fprintf("Maintainence Delta V:          "+sat.propulsion.total_maintainence_dV+" ms^-1\n")
    fprintf("Commission Delta V:            "+sat.propulsion.commission_dV+" ms^-1\n")
    fprintf("Total Delta V (No SF):         "+sat.propulsion.total_dV+" ms^-1\n")
    fprintf("Total Delta V (SF):            "+sat.propulsion.total_dV_sf+" ms^-1\n")
    fprintf("Reserve Delta V:               "+sat.propulsion.reserve_dV+" ms^-1\n")
    fprintf("ISP:                           "+sat.propulsion.isp+" s\n")
    return
end

% function printer(fileID, statement)
%     fprintf(statement);
%     fprintf(fileID, statement);
%     return
% end
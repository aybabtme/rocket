// walkthrough of "How to design, build and test small liquid-fuel rocket engines (Leroy J. Krzycki, 1967, China Lake, CA)"

mod rocket {
    pub fn oxidizer_flow_rate(total_flow_rate: f64, mixture_ratio: f64) -> f64 {
        return (total_flow_rate * mixture_ratio) / (mixture_ratio + 1.0);
    }

    pub fn fuel_flow_rate(total_flow_rate: f64, mixture_ratio: f64) -> f64 {
        return total_flow_rate / (mixture_ratio + 1.0);
    }

    pub fn total_propellant_flow_rate(thrust: f64, specific_impulse: f64) -> f64 {
        return thrust / specific_impulse;
    }

    pub fn gas_constant(universal_gas_constant: f64, molecular_weigth: f64) -> f64 {
        return universal_gas_constant / molecular_weigth;
    }

    pub fn temperature_at_throat(
        combustion_chamber_flame_temperature: f64, // degrees Rankine
        ratio_gas_specific_heat: f64,
    ) -> f64 {
        let Tc = combustion_chamber_flame_temperature;
        let ɣ = ratio_gas_specific_heat;
        return Tc * (1.0 / (1.0 + (ɣ - 1.0) / 2.0));
    }

    pub fn pressure_at_throat(
        combustion_chamber_pressure: f64,
        ratio_gas_specific_heat: f64,
    ) -> f64 {
        let Pc = combustion_chamber_pressure;
        let ɣ = ratio_gas_specific_heat;
        return Pc * (1.0 + (ɣ - 1.0) / 2.0).powf(-ɣ / (ɣ - 1.0));
    }

    pub fn squared_mach_number_at_exit(
        combustion_chamber_pressure: f64,
        atmospheric_pressure: f64,
        ratio_gas_specific_heat: f64,
    ) -> f64 {
        let P_c = combustion_chamber_pressure;
        let P_atm = atmospheric_pressure;
        let ɣ = ratio_gas_specific_heat;

        let bracket_factor = (2.0 / (ɣ - 1.0));

        let P_c_over_P_atm = P_c / P_atm;
        let P_c_over_P_atm_exp = (ɣ - 1.0) / ɣ;
        let bracket = P_c_over_P_atm.powf(P_c_over_P_atm_exp) - 1.0;
        return bracket_factor * bracket;
    }

    pub fn nozzle_throat_cross_sectional_area(
        total_propellant_flow_rate: f64,
        pressure_at_throat: f64,
        gas_constant: f64,
        temperature_at_throat: f64,
        ratio_gas_specific_heat: f64,
        acceleration_at_earth_surface: f64,
    ) -> f64 {
        let ω_t = total_propellant_flow_rate;
        let P_t = pressure_at_throat;
        let R = gas_constant;
        let T_t = temperature_at_throat;
        let ɣ = ratio_gas_specific_heat;
        let g_c = acceleration_at_earth_surface;

        let a = ω_t / P_t;
        let b = (R * T_t) / (ɣ * g_c);
        return a * sqrt(b);
    }

    pub fn nozzle_exit_cross_sectional_area(
        nozzle_throat_cross_sectional_area: f64,
        squared_mach_number_at_exit: f64,
        ratio_gas_specific_heat: f64,
    ) -> f64 {
        let A_t = nozzle_throat_cross_sectional_area;
        let M2_e = squared_mach_number_at_exit;
        let ɣ = ratio_gas_specific_heat;

        let M_e = squared_mach_number_at_exit.sqrt();

        let eqn_factor = (A_t / M_e);
        let eqn_num = 1.0 + ((ɣ - 1.0) / 2.0) * M2_e;
        let eqn_denom = (ɣ + 1.0) / 2.0;
        let eqn_exponent = (ɣ + 1.0) / (2.0 * (ɣ - 1.0));
        return eqn_factor * (eqn_num / eqn_denom).powf(eqn_exponent);
    }

    fn sqrt(v: f64) -> f64 {
        return v.sqrt();
    }
}

#[cfg(test)]
mod tests {
    use crate::rocket::{
        nozzle_throat_cross_sectional_area, pressure_at_throat, temperature_at_throat,
        total_propellant_flow_rate,
    };

    use super::*;

    #[test]
    fn test_oxidizer_flow_rate() {
        let want = 0.293;
        let got = rocket::oxidizer_flow_rate(0.41, 2.5);
        assert_in_range(want, got, 0.001)
    }

    #[test]
    fn test_fuel_flow_rate() {
        let want = 0.117;
        let got = rocket::fuel_flow_rate(0.41, 2.5);
        assert_in_range(want, got, 0.001)
    }

    #[test]
    fn test_total_propellant_flow_rate() {
        let want = 0.41;
        let got = rocket::total_propellant_flow_rate(100.0, 244.0);
        assert_in_range(want, got, 0.001)
    }

    #[test]
    fn test_gas_constant() {
        const UNIVERSAL_GAS_CONSTANT: f64 = 1542.32; // ft-lb/lb*R
        let want = 65.0;
        let got = rocket::gas_constant(UNIVERSAL_GAS_CONSTANT, 23.7);
        assert_in_range(want, got, 0.1);
    }

    #[test]
    fn test_temperature_at_throat() {
        let want = 0.909;
        let got = rocket::temperature_at_throat(1.0, 1.2);
        assert_in_range(want, got, 0.001);
    }

    #[test]
    fn test_pressure_at_throat() {
        let want = 0.564;
        let got = rocket::pressure_at_throat(1.0, 1.2);
        assert_in_range(want, got, 0.001);
    }

    #[test]
    fn test_mach_number_at_exit() {
        let combustion_chamber_pressure = 300.0; // psi
        let atmospheric_pressure = 14.7; // psi
        let ratio_gas_specific_heat = 1.2;
        let want = 6.531;
        let got = rocket::squared_mach_number_at_exit(
            combustion_chamber_pressure,
            atmospheric_pressure,
            ratio_gas_specific_heat,
        );
        assert_in_range(want, got, 0.001);
    }

    #[test]
    fn test_nozzle_throat_cross_sectional_area() {
        let total_propellant_flow_rate = 0.077;
        let gas_constant = 65.0;
        let pressure_at_throat = 169.0;
        let temperature_at_throat = 5650.0;
        let ratio_gas_specific_heat = 1.2;
        let acceleration_at_earth_surface = 32.2; // ft/sec^2

        let want = 0.044; // in^2
        let got = rocket::nozzle_throat_cross_sectional_area(
            total_propellant_flow_rate,
            pressure_at_throat,
            gas_constant,
            temperature_at_throat,
            ratio_gas_specific_heat,
            acceleration_at_earth_surface,
        );
        assert_in_range(want, got, 0.01);
    }

    #[test]
    fn test_nozzle_overall_parameters() {
        let universal_gas_constant = 1542.32; // ft-lb/lb*R
        let specific_impulse = 244.0; // for gas.oxygen/gasoline @ 200PSI
        let oxy_gasoline_comb_molecular_weight = 23.7; // "molecular weight of the hot gaseous products of combustion of gaseous oxygen/hydrocarbon fuel is about 24"
        let ratio_gas_specific_heat = 1.2; // for gas oxygen/hydrocarbon

        let pressure_at_earth_surface = 14.7; // psi
        let acceleration_at_earth_surface = 32.2; // ft/sec^2

        struct Condition {
            thrust: f64,
            P_c: f64,

            Isp: f64,
            comb_product_molecular_weight: f64, // molecular weight of hot-gaseous product from combustion
            ɣ: f64,
            P_atm: f64,
            T_c: f64,

            M_e: f64,
            A_e_over_A_t: f64,
            T_e_over_T_c: f64,
        }
        let tests: [Condition; 5] = [
            Condition {
                thrust: 100.0,
                P_c: 100.0,

                Isp: specific_impulse,
                comb_product_molecular_weight: oxy_gasoline_comb_molecular_weight,
                ɣ: ratio_gas_specific_heat,
                P_atm: pressure_at_earth_surface,
                T_c: 5500.0 + 460.0,

                M_e: 1.95,
                A_e_over_A_t: 1.79,
                T_e_over_T_c: 0.725,
            },
            Condition {
                thrust: 100.0,
                P_c: 200.0,

                Isp: specific_impulse,
                comb_product_molecular_weight: oxy_gasoline_comb_molecular_weight,
                ɣ: ratio_gas_specific_heat,
                P_atm: pressure_at_earth_surface,
                T_c: 5700.0 + 460.0,

                M_e: 2.34,
                A_e_over_A_t: 2.74,
                T_e_over_T_c: 0.65,
            },
            Condition {
                thrust: 100.0,
                P_c: 300.0,

                Isp: specific_impulse,
                comb_product_molecular_weight: oxy_gasoline_comb_molecular_weight,
                ɣ: ratio_gas_specific_heat,
                P_atm: pressure_at_earth_surface,
                T_c: 5742.0 + 460.0,

                M_e: 2.55,
                A_e_over_A_t: 3.65,
                T_e_over_T_c: 0.606,
            },
            Condition {
                thrust: 100.0,
                P_c: 400.0,

                Isp: specific_impulse,
                comb_product_molecular_weight: oxy_gasoline_comb_molecular_weight,
                ɣ: ratio_gas_specific_heat,
                P_atm: pressure_at_earth_surface,
                T_c: 5812.0 + 460.0,

                M_e: 2.73,
                A_e_over_A_t: 4.6,
                T_e_over_T_c: 0.574,
            },
            Condition {
                thrust: 100.0,
                P_c: 500.0,

                Isp: specific_impulse,
                comb_product_molecular_weight: oxy_gasoline_comb_molecular_weight,
                ɣ: ratio_gas_specific_heat,
                P_atm: pressure_at_earth_surface,
                T_c: 5862.0 + 460.0,

                M_e: 2.83,
                A_e_over_A_t: 5.28,
                T_e_over_T_c: 0.55,
            },
        ];

        for tt in tests {
            let thrust = tt.thrust;
            let Isp = tt.Isp;
            let ɣ = tt.ɣ;
            let P_atm = tt.P_atm;
            let P_c = tt.P_c;
            let T_c = tt.T_c;

            let ω_t = rocket::total_propellant_flow_rate(thrust, Isp);
            let P_t = rocket::pressure_at_throat(P_c, ɣ);
            let T_t = rocket::temperature_at_throat(T_c, ɣ);
            let R = rocket::gas_constant(universal_gas_constant, tt.comb_product_molecular_weight);

            let A_t = rocket::nozzle_throat_cross_sectional_area(
                ω_t,
                P_t,
                R,
                T_t,
                ɣ,
                acceleration_at_earth_surface,
            );

            let M2_e = rocket::squared_mach_number_at_exit(P_c, P_atm, ɣ);
            let M_e = M2_e.sqrt();

            let A_e = rocket::nozzle_exit_cross_sectional_area(A_t, M2_e, ɣ);

            let want_M_e = tt.M_e;
            let want_A_e_over_A_t = tt.A_e_over_A_t;
            let want_T_e_over_T_c = tt.T_e_over_T_c;

            let got_M_e = M2_e.sqrt();
            let got_A_e_over_A_t = A_e / A_t;
            let got_T_e_over_T_c = 1.0; // how to compute T_e?

            assert_in_range(want_M_e, got_M_e, 0.01);
            assert_in_range(want_A_e_over_A_t, got_A_e_over_A_t, 0.001);
        }
    }

    fn assert_in_range(want: f64, got: f64, epsilon: f64) {
        if (want - got).abs() < epsilon {
            return;
        }
        panic!(
            "assertion failed: `(left ~= right) not within epsilon {}
        left: `{}`
       right: `{}",
            epsilon, want, got
        )
    }
}

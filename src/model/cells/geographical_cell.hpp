// Created by binybrion - 06/29/20
// Modified by Glenn    - 02/07/20
// And by Eric          - Summer/2021

#ifndef PANDEMIC_HOYA_2002_ZHONG_CELL_HPP
#define PANDEMIC_HOYA_2002_ZHONG_CELL_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <cadmium/celldevs/cell/cell.hpp>
#include <iomanip>
#include "vicinity.hpp"
#include "sevirds.hpp"
#include "simulation_config.hpp"
#include "AgeData.hpp"
#include "../Helpers/Assert.hpp"

using namespace std;
using namespace cadmium::celldevs;
using namespace Assert;

// Indices for a vector used in local_compute(), compute_vaccinated(),
// and compute_EIRD()
unsigned int const NVAC = 0;
unsigned int const VAC1 = 1;
unsigned int const VAC2 = 2;
unsigned int const BOOS = 3;

template <typename T>
class geographical_cell : public cell<T, string, sevirds, vicinity>
{
    public:
        template <typename X>
        using cell_unordered = unordered_map<string, X>;

        using cell<T, string, sevirds, vicinity>::simulation_clock;
        using cell<T, string, sevirds, vicinity>::state;
        using cell<T, string, sevirds, vicinity>::neighbors;
        using cell<T, string, sevirds, vicinity>::cell_id;

        using config_type = simulation_config;

        using phase_rates = vector<     // The age sub_division
                            vecDouble>; // The stage of infection

        phase_rates virulence_rates;
        phase_rates incubationD1_rates;
        phase_rates incubationD2_rates;
        phase_rates incubation_rates;
        phase_rates recovery_rates;
        phase_rates recoveryD1_rates;
        phase_rates recoveryD2_rates;
        phase_rates mobility_rates;
        phase_rates fatality_rates;
        phase_rates fatalityD1_rates;
        phase_rates fatalityD2_rates;
        phase_rates vac1_rates;
        phase_rates vac2_rates;

        // To make the parameters of the correction_factors variable more obvious
        using infection_threshold        = float;
        using mobility_correction_factor = array<float, 2>;  // array<mobility correction factor, hysteresis factor>;

        bool reSusceptibility, is_vaccination;

        unsigned int age_segments;

        geographical_cell() : cell<T, string, sevirds, vicinity>() {}

        geographical_cell(string const& cell_id, cell_unordered<vicinity> const& neighborhood,
                            sevirds const& initial_state, string const& delay_id, simulation_config config) :
            cell<T, string, sevirds, vicinity>(cell_id, neighborhood, initial_state, delay_id)
        {
            for (const auto& i : neighborhood)
                state.current_state.hysteresis_factors.insert({i.first, hysteresis_factor{}});

            // Set whether or not vaccines are being modeled
            // to be used in the getters found in sevirds.hpp
            // and later in this file
            is_vaccination               = config.is_vaccination;
            state.current_state.vaccines = is_vaccination;

            // Set the precision divider in the sevirds object
            state.current_state.prec_divider          = (double)config.prec_divider;
            state.current_state.one_over_prec_divider = 1.0 / (double)config.prec_divider;

            virulence_rates  = move(config.virulence_rates);
            incubation_rates = move(config.incubation_rates);
            recovery_rates   = move(config.recovery_rates);
            mobility_rates   = move(config.mobility_rates);
            fatality_rates   = move(config.fatality_rates);

            // Multiplication is always faster then division so set this up to be 1/prec_divider to be multiplied later
            reSusceptibility  = config.reSusceptibility;
            age_segments = initial_state.get_num_age_segments();

            if (is_vaccination)
            {
                vac1_rates = move(config.vac1_rates);
                vac2_rates = move(config.vac2_rates);

                incubationD1_rates = move(config.incubationD1_rates);
                incubationD2_rates = move(config.incubationD2_rates);

                recoveryD1_rates = move(config.recovery_ratesD1);
                recoveryD2_rates = move(config.recovery_ratesD2);

                fatalityD1_rates = move(config.fatality_ratesD1);
                fatalityD2_rates = move(config.fatality_ratesD2);
            }
        }

        /**
         * @brief This is the 'main' function for the class
         * and is where all the equations for the the current cell
         * and on the current day are computed for each age group
         * 
         * @return sevirds
        */
        sevirds local_computation() const override
        {
            // Can't be a reference since it would need to be
            // const and then we wouldn't be allowed to change its values
            sevirds res = state.current_state;

            // Number of AgeData objects needed
            // One for non-vac, dose1, dose2, and any booster shot populations
            int size = 1;
            if (is_vaccination)
                size += 2;

            // Initialize them in a vector for easy moving around the functions
            vector<unique_ptr<AgeData>> datas;
            for (int i = 0; i < size; ++i)
                datas.push_back(unique_ptr<AgeData>(nullptr));

            // Global new susceptible variable as the other equations
            // remove their proportions from this one leaving it with
            // the remaning susceptible proportion
            double new_s;

            // Calculate the next new sevirds variables for each age group
            for (unsigned int age_segment_index = 0; age_segment_index < age_segments; ++age_segment_index)
            {
                // Reset for susceptible equation
                new_s = 1;

                // Init the non-vac object for the current age group
                datas.at(NVAC).reset(new AgeData(age_segment_index, res.susceptible, res.exposed, res.infected,
                                                res.recovered, incubation_rates, recovery_rates, fatality_rates));

                if (is_vaccination)
                {
                    // Init the vac object for the current age group
                    datas.at(VAC1).reset(new AgeData(age_segment_index, res.vaccinatedD1, res.exposedD1, res.infectedD1,
                                                    res.recoveredD1, incubationD1_rates, recoveryD1_rates,
                                                    fatalityD1_rates, vac1_rates.at(age_segment_index),
                                                    res.immunityD1_rate.at(age_segment_index), AgeData::PopType::DOSE1));
                    datas.at(VAC2).reset(new AgeData(age_segment_index, res.vaccinatedD2, res.exposedD2, res.infectedD2,
                                                    res.recoveredD2, incubationD2_rates, recoveryD2_rates,
                                                    fatalityD2_rates, vac2_rates.at(age_segment_index),
                                                    res.immunityD2_rate.at(age_segment_index), AgeData::PopType::DOSE2));

                    // Equations for Vaccinated population (eg. EV1, RV2...)
                    sanity_check(res.get_total_susceptible(true, age_segment_index), __LINE__);
                    compute_vaccinated(datas, res);

                    // S = 1 - V1 - V2
                    new_s -= datas.at(VAC1).get()->GetTotalSusceptible(); // 1e
                    sanity_check(new_s, __LINE__);
                    new_s -= datas.at(VAC2).get()->GetTotalSusceptible(); // 2d
                    sanity_check(new_s, __LINE__);
                }

                // Compute the Exposed, Infected, Recovered, and Fatalities equations
                // for all population types
                compute_EIRD(datas, res);

                // S = 1 - E - I - R - F
                for (unique_ptr<AgeData>& data : datas)
                {
                    new_s -= data.get()->GetTotalExposed();
                    sanity_check(new_s, __LINE__);

                    new_s -= data.get()->GetTotalInfected();
                    sanity_check(new_s, __LINE__);

                    new_s -= data.get()->GetTotalRecovered();
                    sanity_check(new_s, __LINE__);

                    // Using the total is currently broken but is the ideal as it removes the need
                    // to loop every cycle
                    //res.fatalities.at(age_segment_index) += data.get()->GetTotalFatalities();
                    // For some reason doing this loop gives the correct result when the Total does not
                    for (unsigned int q = 0; q < data.get()->GetInfectedPhase(); ++q)
                        res.fatalities.at(age_segment_index) += data.get()->GetNewFatalities(q);
                    sanity_check(res.fatalities.at(age_segment_index), __LINE__);
                }

                new_s -= res.fatalities.at(age_segment_index);
                sanity_check(new_s, __LINE__);

                res.susceptible.at(age_segment_index).front() = new_s;
            } //for(age_groups)

            return res;
        } //local_computation()

        // It returns the delay to communicate cell's new state.
        // It looks useless but it is extremely important. Do NOT delete!
        T output_delay(sevirds const& cell_state) const override { return 1; }

        /**
         * @brief Vaccinated Dose 1 - Equation 1a
         * 
         * @param datas Vector containing the three population types and their data
         * @param res State machine object that holds simulation config data
         * @return double
         */
        double new_vaccinated1(vector<unique_ptr<AgeData>>& datas, sevirds const& res) const
        {
            // Vaccination rate with those who are susceptible
            // vd1 * S
            double new_vac1 = datas.at(VAC1).get()->GetVaccinationRate(0)  // vd1
                            * datas.at(NVAC).get()->GetOrigSusceptible(0); // * S

            // And those who are in the recovery phase
            double sum = 0;
            for (unsigned int q = datas.at(NVAC).get()->GetRecoveredPhase() - 1; q > res.min_interval_recovery_to_vaccine; --q)
            {
                // Remember these values in the non-vac object as
                // they are removed from the susceptible group
                // in increment_recoveries(). Only do math once!!
                datas.at(NVAC).get()->SetVacFromRec(q - 1,
                                                    datas.at(NVAC).get()->GetOrigRecovered(q - 1) // R(q)
                                                    * datas.at(VAC1).get()->GetVaccinationRate(0) // vd1
                );

                sum += datas.at(NVAC).get()->GetVacFromRec(q - 1);
            }

            return new_vac1 + sum;
        }

        /**
         * @brief Vaccinated Dose 2 - Equation 2a
         * 
         * @param datas Vector containing the three population types with their respective data
         * @param res Current state of the cell
         * @return double
         */
        double new_vaccinated2(vector<unique_ptr<AgeData>>& datas, sevirds& res, vecDouble const& earlyVac2) const
        {
            AgeData& age_data_vac1 = *(datas.at(VAC1)).get();
            AgeData& age_data_vac2 = *(datas.at(VAC2)).get();

            // Everybody on the last day of dose 1 is moved to dose 2
            double vac2 = age_data_vac1.GetOrigSusceptibleBack(); // V1(td1)

            // Some people are eligible to receive their second dose sooner
            // and this was already computed ealier in compute_vaccinated()
            // qϵ{mtd1...td1 - 1}
            vac2 += accumulate(earlyVac2.begin(), earlyVac2.end(), 0.0);

            // Some people are eligible to receive their second dose sooner from the dose 1 recovery pop
            // qϵ{mtd1...Tr}
            for (unsigned int q = age_data_vac1.GetRecoveredPhase(); q > res.min_interval_recovery_to_vaccine; --q)
            {
                // Remember these values for when they are removed from the
                // vac1 susceptible group in increment_recoveries()
                age_data_vac1.SetVacFromRec(q - 1,
                                            age_data_vac2.GetVaccinationRate(q - 1 - res.min_interval_recovery_to_vaccine) // v(q)
                                                * age_data_vac1.GetOrigRecovered(q - 1)                                    // RV1(q)
                );

                vac2 += age_data_vac1.GetVacFromRec(q - 1);
            }

            // - V1(td1) * sum(1...k and 1...Ti))
            return vac2 - new_exposed(res, *(datas.at(VAC1).get()), age_data_vac1.GetSusceptiblePhase());
        }

        /**
         * @brief Calculates proportion of new exposures from either non-vac or vac (dose 1 or 2) population.
         * 1b, 1c, 1d, 1e, 1f, 2b, 2c, 2d, 2e, 3a, 3b and 3c use this
         * 
         * @param res State machine object that holds simulation config data
         * @param age_data Reference to current simulation data
         * @param q Index to compute equation
         * @return double
        */
        double new_exposed(sevirds& res, AgeData& age_data, int q=0) const
        {
            double expos = 0, sum = 0, inner_sum, inner_sumV1, inner_sumV2;

            // Calculate the correction factor of the current cell.
            // The current cell must be part of its own neighborhood for this to work!
            vicinity self_vicinity = state.neighbors_vicinity.at(cell_id);
            double current_cell_correction_factor = res.disobedient
                                                    + (1 - res.disobedient)
                                                    * movement_correction_factor(self_vicinity.correction_factors,
                                                                                state.neighbors_state.at(cell_id).get_total_infections(),
                                                                                res.hysteresis_factors.at(cell_id));

            double neighbor_correction;

            // jϵ{1...k}
            for (string neighbor : neighbors)
            {
                sevirds const& nstate = state.neighbors_state.at(neighbor);       // Cell j's state
                vicinity const& v     = state.neighbors_vicinity.at(neighbor);    // Holds cij and a correction factor used in kij

                // Disobedient people have a correction factor of 1. The rest of the population is affected by the movement_correction_factor
                neighbor_correction = nstate.disobedient
                                        + (1 - nstate.disobedient)
                                        * movement_correction_factor(v.correction_factors,
                                                                    nstate.get_total_infections(),
                                                                    res.hysteresis_factors.at(neighbor));

                // Logically makes sense to require neighboring cells to follow the movement restriction that is currently
                // in place in the current cell if the current cell has a more restrictive movement.
                neighbor_correction = min(current_cell_correction_factor, neighbor_correction);

                // Reset the inner sum for the next neighbor
                inner_sum = 0; inner_sumV1 = 0; inner_sumV2 = 0;

                // bϵ{1...A}
                for (unsigned int age_group = 0; age_group < nstate.num_age_groups; ++age_group)
                {

                    // nϵ{1...Ti}
                    for (unsigned int n = 0; n < nstate.infected.at(age_group).size(); ++n)
                    {
                        inner_sum +=
                            mobility_rates.at(age_group).at(n)    // μ(n)
                            * virulence_rates.at(age_group).at(n) // λ(n)
                            * nstate.infected.at(age_group).at(n) // I(n)
                            ;
                    }

                    if (is_vaccination)
                    {
                        // nϵ{1...Ti,V1}
                        for (unsigned int n = 0; n < nstate.infectedD1.at(age_group).size(); ++n)
                        {
                            inner_sumV1 +=
                                mobility_rates.at(age_group).at(n)      // μ(n)
                                * virulence_rates.at(age_group).at(n)   // λ(n)
                                * nstate.infectedD1.at(age_group).at(n) // IV1(n)
                                ;
                        }

                        // nϵ{1...Ti,V2}
                        for (unsigned int n = 0; n < nstate.infectedD2.at(age_group).size(); ++n)
                        {
                            inner_sumV2 +=
                                mobility_rates.at(age_group).at(n)      // μ(n)
                                * virulence_rates.at(age_group).at(n)   // λ(n)
                                * nstate.infectedD2.at(age_group).at(n) // IV2(n)
                                ;
                        }
                    }

                    sum += v.correlation                                // cij
                           * neighbor_correction                        // kij
                           * (inner_sum + inner_sumV1 + inner_sumV2)    // sum(1...Ti)
                           * nstate.age_group_proportions.at(age_group) // Njb / Nj
                        ;
                }
            }

            expos = age_data.GetOrigSusceptible(q) * sum; // S * sum(1...k)

            if (age_data.GetType() != AgeData::PopType::NVAC)
                expos *= 1.0 - age_data.GetImmunityRate( int((q - 1) * 0.14f) ); // 1 - i(q)

            sanity_check(expos, __LINE__);
            return expos;
        } //new_exposed()

        /**
         * @brief Exposed: E(q), EV1(q), EV2(q)
         *  Advance all exposed forward a day, with some proportion leaving exposed(q-1) and entering infected(1)
         * 
         * @param age_data Reference to current simulation data for current age group (nvac, dose1, dose2) and age group
        */
        void increment_exposed(AgeData& age_data) const
        {
            double curr_expos;

            // qϵ{2...Te}
            for (unsigned int q = age_data.GetExposedPhase(); q > 0; --q)
            {
                // Moves each proportion group in the phase to the next day and removes
                // those who become infected earlier via the incubation rate
                curr_expos = (1 - age_data.GetIncubationRate(q - 1)) // 1 - ε(q - 1)
                             * age_data.GetOrigExposed(q - 1)        // * E(q - 1)
                    ;

                sanity_check(curr_expos, __LINE__);
                age_data.SetExposed(q, curr_expos);
            }
        }

        /**
         * @brief Infection: I(1), IV1(1), or IV2(1)
         *  Calculates proportion of new infections from either non-vac or vac (dose 1 or 2) population
         * 
         * @param age_data Reference to current simulation data
         * @return double
        */
        double new_infections(AgeData& age_data) const
        {
            double inf = 0;

            /* Scan through all exposed days and calculate exposed.at(age).at(q)
            *   Incubation Rate on Te must be 1.0
            *   Note: age_data.get()->GetOrigExposed(i) == exposed.at(age).at(q) 
            *   and at timestep t not t+1
            *   qϵ{1...Te-1}
            */
            for (unsigned int q = 1; q <= age_data.GetExposedPhase(); ++q)
            {
                // Calculates those who move early to the infected phase
                // and automatically moves those on the last day to the infected phase
                inf += age_data.GetIncubationRate(q) // ε(q), εV1(q), or εV2(q)
                       * age_data.GetOrigExposed(q)  // E(q), EV1(q), or EV2(q)
                    ;
            }

            sanity_check(inf, __LINE__);
            return inf;
        }

        /**
         * @brief Infectd: I(q), IV1(q), IV2(q)
         *  Advances all infected forward a day, with some already moved to fatalities or recovered prior
         * 
         * @param age_data Reference to current simulation data for current age group (nvac, dose1, dose2) and age group
         * @param recovered Vector of new recoveries from each day
        */
        void increment_infections(AgeData& age_data) const
        {
            double curr_inf;

            // qϵ{2...Ti}
            for (unsigned int q = age_data.GetInfectedPhase(); q > 0; --q)
            {
                // The previous day of infections minus those
                // who have died and those who have recovered
                curr_inf = age_data.GetOrigInfected(q - 1)    // I(q - 1)
                           - age_data.GetNewFatalities(q - 1) // - D(q - 1)
                           - age_data.GetNewRecovered(q - 1)  // - R(q - 1)
                    ;

                sanity_check(curr_inf, __LINE__);
                age_data.SetInfected(q, curr_inf);
            }
        }

        /**
         * @brief Recovered: R(1), RV1(1), or RV2(1)
         *  Calculates proportion of new new recoveries from either non-vac or vac (dose 1 or 2) population
         * 
         * @param age_data Reference to simulation data for the current age group and population type
         * @return double
        */
        double new_recoveries(AgeData& age_data) const
        {
            // Assume that any individuals that are not fatalities on the last stage of infection recover
            double recoveries = age_data.GetOrigInfectedBack()    // I(q)
                                - age_data.GetNewFatalitiesBack() // - D(q)
                ;

            sanity_check(recoveries, __LINE__);
            age_data.SetNewRecovered(age_data.GetInfectedPhase(), recoveries);

            // qϵ{1...Ti - 1}
            double sum;
            for (unsigned int q = 0; q <= age_data.GetInfectedPhase() - 1; ++q)
            {
                // Calculate all of the new recoveries for every day that a population is infected, some recover
                sum = age_data.GetRecoveryRate(q)   // γ(q)
                    * age_data.GetOrigInfected(q) // I(q)
                    ;

                recoveries += sum;
                age_data.SetNewRecovered(q, sum);
            }

            sanity_check(recoveries, __LINE__);
            return recoveries;
        }

        /**
         * @brief Infectd: R(q), RV1(q), RV2(q)
         *  Advances all infected forward a day, with some already moved to fatalities or recovered prior
         * 
         * @param age_data Reference to current simulation data for current age group (nvac, dose1, dose2) and age group
         * @param recovered_index If res-susc is turned on this will avoid processing the population on the last day
         * @param age_data_vac Pointer to a vaccinated age_data object that is used for R(q) and RV1(q)
         * @param res Used to get the minimum interval between doses needed in RV1(q)
        */
        void increment_recoveries(AgeData& age_data) const
        {
            double curr_rec;

            // qϵ{2...Tr}
            for (unsigned int q = age_data.GetRecoveredPhase(); q > 0; --q)
            {
                curr_rec = 0;

                // When resusceptibility is off then those who are recovered stay in that phase
                if (!reSusceptibility && q == age_data.GetRecoveredPhase())
                    curr_rec += age_data.GetRecoveredBack();

                // Each day of the recovered phase is the value of the previous day. The population on the last day is
                // now susceptible (assuming a re-susceptible model); this is implicitly done already as the susceptible value was set to 1.0 and the
                // population on the last day of recovery is never subtracted from the susceptible value.
                // 5d, 5e, 5f
                curr_rec += age_data.GetOrigRecovered(q - 1) - age_data.GetVacFromRec(q - 1); // R(q - 1) * (1 - vd(q - 1))

                sanity_check(curr_rec, __LINE__);
                age_data.SetRecovered(q, curr_rec);
            }
        }

        /**
         * @brief New fatalities for one population group (Non-Vaccinated/Vaccinated Dose 1/Vaccinvated Dose 2)
         *  These are calculated individually as the vector of new fatalities is used by the functions
         *  that follow this one (see compute_vaccinated() and compute_not_vaccinated()). If we calculated
         *  the global number of fatalities later when we do something like dose1.infected().at(q) - fatalities.at(q)
         *  to get the maximum number of possible recoveries the fatality number will be including those from dose 2 and
         *  not vaccinated when in reality the maximum number of recoveries for dose1 is limited to those who are infected
         *  with dose 1 and still alive so we only want to remove those who are dose 1 fatality.
         *
         * @param res State of the geographical cell (holds some global data)
         * @param age_data Contains the data of the proportion. In this function the infections proportion as well as
         *                  the fatality rates are used from here
         * @return double
        */
        double new_fatalities(sevirds const& res, AgeData& age_data) const
        {
            double new_f = 0.0, sum;

            // Calculate all those who have died during an infection stage.
            // qϵ{1...Ti}
            for (unsigned int q = 0; q <= age_data.GetInfectedPhase(); ++q)
            {
                // fa(q) * I(q)
                sum = age_data.GetFatalityRate(q) * age_data.GetOrigInfected(q);

                // Amplify fatality rate if the hospitals are full
                if (res.get_total_infections() > res.hospital_capacity)
                    sum *= res.fatality_modifier;

                new_f += sum;
                age_data.SetNewFatalities(q, sum);
            }

            sanity_check(new_f, __LINE__);
            return new_f;
        }

        double movement_correction_factor(const map<infection_threshold, mobility_correction_factor>& mobility_correction_factors,
                                        double infectious_population, hysteresis_factor& hysteresisFactor) const
        {
            // For example, assume a correction factor of "0.4": [0.2, 0.1]. If the infection goes above 0.4, then the
            // correction factor of 0.2 will now be applied to total infection values above 0.3, no longer 0.4 as the
            // hysteresis is in effect.
            if (infectious_population > hysteresisFactor.infections_higher_bound)
                hysteresisFactor.in_effect = false;

            // This is uses the comparison '>', not '>=' ; otherwise if the lower bound is 0 there is no way for the hysteresis
            // to disappear as the infections can never go below 0
            if (hysteresisFactor.in_effect && infectious_population > hysteresisFactor.infections_lower_bound)
                return hysteresisFactor.mobility_correction_factor;

            hysteresisFactor.in_effect = false;

            double correction = 1.0;
            for (auto const& pair: mobility_correction_factors)
            {
                if (infectious_population >= pair.first)
                {
                    correction = pair.second.front();

                    // A hysteresis factor will be in effect until the total infection goes below the hysteresis factor;
                    // until that happens the information required to return a movement factor must be kept in above variables.

                    // Get the threshold of the next correction factor; otherwise the current correction factor can remain in
                    // effect if the total infections never goes below the lower bound hysteresis factor, but also if it goes
                    // above the original total infection threshold!
                    auto next_pair_iterator = find(mobility_correction_factors.begin(), mobility_correction_factors.end(), pair);
                    AssertLong(next_pair_iterator != mobility_correction_factors.end(), __FILE__, __LINE__);

                    // If there is a next correction factor (for a higher total infection), then use it's total infection threshold
                    if ((long unsigned int) distance(mobility_correction_factors.begin(), next_pair_iterator) != mobility_correction_factors.size() - 1)
                        ++next_pair_iterator;

                    hysteresisFactor.in_effect                  = true;
                    hysteresisFactor.infections_higher_bound    = next_pair_iterator->first;
                    hysteresisFactor.infections_lower_bound     = pair.first - pair.second.back();
                    hysteresisFactor.mobility_correction_factor = pair.second.front();
                } else
                    break;
            }

            return correction;
        } //movement_correction_factor()

        /**
         * @brief Computes all the equations specific to the vaccinated population
         * 
         * @param datas Vector of AgeData objects containing current age group data
         * @param res The current state of the geographical cell
        */
        void compute_vaccinated(vector<unique_ptr<AgeData>>& datas, sevirds& res) const
        {
            double curr_vac1 = 0.0, curr_vac2 = 0.0;

            AgeData& age_data_vac1 = *(datas.at(VAC1).get());
            AgeData& age_data_vac2 = *(datas.at(VAC2).get());

            // Holds those who get their second dose earlier from the susceptible dose 1 group
            // This is not the same as vacFromRec in AgeData.hpp
            vecDouble earlyVac2(age_data_vac1.GetSusceptiblePhase(), 0.0);

            // <VACCINATED DOSE 1>
                // Calculate the number of new vaccinated dose 1
                double new_vac1 = new_vaccinated1(datas, res); // 1a

                // qϵ{2...td1}
                for (unsigned int q = age_data_vac1.GetSusceptiblePhase(); q > 0; --q)
                {
                    // 1b & 1d
                    curr_vac1 = age_data_vac1.GetOrigSusceptible(q - 1); // V1(q - 1)

                    age_data_vac1.SetNewExposed(q, new_exposed(res, age_data_vac1, q - 1));
                    curr_vac1 -= age_data_vac1.GetNewExposed(q); // - ( V1(q - 1) * (1 - iv1(q - 1)) * sum(1..k and 1...Ti) )

                    // Early dose 2
                    if (q > res.min_interval_doses)
                    {
                        // 1d
                        if (q > res.min_interval_recovery_to_vaccine)
                            earlyVac2.at(q - 1) = age_data_vac2.GetVaccinationRate(q - 1 - res.min_interval_recovery_to_vaccine) // vd2(q - 1)
                                                * age_data_vac1.GetOrigSusceptible(q - 1)                                        // * V1(q - 1)
                            ;
                        // 1c substracts early dose2 vaccinations from 1b
                        else
                            earlyVac2.at(q - 1) = age_data_vac2.GetVaccinationRate(q - 1 - res.min_interval_doses) // vd2(q - 1)
                                                * age_data_vac1.GetOrigSusceptible(q - 1)                          // * V1(q - 1)
                            ;

                        curr_vac1 -= earlyVac2.at(q - 1);
                    }

                    sanity_check(curr_vac1, __LINE__);

                    // Update the current day with the modified exposed from yesterday
                    age_data_vac1.SetSusceptible(q, curr_vac1);
                }

                // 1d
                if (reSusceptibility)
                {
                    double susc_from_rec = age_data_vac1.GetOrigRecoveredBack()                                                                             // RV1(Tr)
                                        * (1 - age_data_vac2.GetVaccinationRate(age_data_vac1.GetRecoveredPhase() - res.min_interval_recovery_to_vaccine)); // * (1 - vd2(Tr))
                    age_data_vac1.AddSusceptibleBack(susc_from_rec);
                }

                // Set the new dose1 proportion to the beginning of the phase
                age_data_vac1.SetSusceptible(0, new_vac1);
                sanity_check(age_data_vac1.GetTotalSusceptible(), __LINE__);
            // </VACCINATED DOSE 1>

            // <VACCINATED DOSE 2>
                // Calculate the number of new vaccinated dose 2
                double new_vac2 = new_vaccinated2(datas, res, earlyVac2);
                sanity_check(new_vac2, __LINE__);

                // qϵ{2...td2 - 1}
                for (unsigned int q = age_data_vac2.GetSusceptiblePhase() - 1; q > 0; --q)
                {
                    // 2b
                    curr_vac2 = age_data_vac2.GetOrigSusceptible(q - 1); // V2(q - 1)

                    age_data_vac2.SetNewExposed(q, new_exposed(res, age_data_vac2, q - 1));
                    curr_vac2 -= age_data_vac2.GetNewExposed(q); // - V2(q - 1) * (1 - iv2(q - 1)) * sum( jϵ{1…k}(cij * kij * sum(bϵ{1...A} and nϵ{1...Ti}[...])) )

                    sanity_check(curr_vac2, __LINE__);
                    age_data_vac2.SetSusceptible(q, curr_vac2);
                }

                // 2c
                double end = age_data_vac2.GetOrigSusceptible(age_data_vac2.GetSusceptiblePhase() - 1) // V2(td2 - 1)
                      + age_data_vac2.GetOrigSusceptibleBack();                                        // V2(td2)

                age_data_vac2.SetNewExposed(age_data_vac2.GetSusceptiblePhase() - 1, new_exposed(res, age_data_vac2, age_data_vac2.GetSusceptiblePhase() - 1));
                age_data_vac2.SetNewExposed(age_data_vac2.GetSusceptiblePhase(), new_exposed(res, age_data_vac2, age_data_vac2.GetSusceptiblePhase()));
                end -= age_data_vac2.GetNewExposed(age_data_vac2.GetSusceptiblePhase() - 1); // - V2(td2 - 1) * (1 - iV2(td2 - 1)) * sum( jϵ{1...k}(cij * kij * sum(bϵ{1...A} and nϵ{1...Ti}[...]) )
                end -= age_data_vac2.GetNewExposed(age_data_vac2.GetSusceptiblePhase());     // - V2(td2) * (1 - iV2(td2)) * sum( jϵ{1...k}(cij * kij * sum(bϵ{1...A} and nϵ{1...Ti}[...]) )

                if (reSusceptibility)
                    end += age_data_vac2.GetOrigRecoveredBack(); // + RV2(Tr)

                sanity_check(end, __LINE__);
                age_data_vac2.SetSusceptible(age_data_vac2.GetSusceptiblePhase(), end);
                age_data_vac2.SetSusceptible(0, new_vac2); // Set the first day of the phase
                sanity_check(age_data_vac2.GetTotalSusceptible(), __LINE__);
            // </VACCINATED DOSE 2>
        }

        /**
         * @brief Computes the exposed, infected, recovered, and dead equations for all population types
         * Setup the the datas vector to hold all the population types and they'll be looped through
         * 
         * @param datas Vector of pointers holding the population states (i.e., NVac, Dose1, Dose2)
         * @param res Current cell data
         */
        void compute_EIRD(vector<unique_ptr<AgeData>>& datas, sevirds& res) const
        {
            double new_expos, new_inf, new_rec;

            for (unique_ptr<AgeData>& age_data_ptr : datas)
            {
                AgeData& age_data = *(age_data_ptr.get());
                // <FATALITIES>
                    // Calculates the new fatalities on each day of the infected phase
                    // for easy use and less repetive code later
                    age_data.SetTotalFatalities(new_fatalities(res, age_data));
                    sanity_check(age_data.GetTotalFatalities(), __LINE__);
                // </FATALITIES>

                // <RECOVERIES>
                    // Calculates the new recoveries on each day of the infected phase
                    new_rec = new_recoveries(age_data);
                // </RECOVERIES>

                // <EXPOSED>
                    new_expos = 0.0;

                    // qϵ{1...Td2}
                    for (unsigned int q = 0; q <= age_data.GetSusceptiblePhase(); ++q)
                    {
                        if (age_data.GetType() != AgeData::PopType::NVAC)
                            new_expos += age_data.GetNewExposed(q);
                        else
                            new_expos += new_exposed(res, age_data, q);
                    }

                    increment_exposed(age_data);

                    age_data.SetExposed(0, new_expos);
                // </EXPOSED>

                // <INFECTED>
                    new_inf = new_infections(age_data);

                    increment_infections(age_data);

                    age_data.SetInfected(0, new_inf);
                // </INFECTED>

                // <RECOVERED>
                    increment_recoveries(age_data);

                    // The people on the first day of recovery are those that were on the last stage of infection (minus those who died;
                    // already accounted for) in the previous time step plus those that recovered early during an infection stage.
                    age_data.SetRecovered(0, new_rec);
                // </RECOVERED>
            }
        }

        /**
         * @brief Basic check that the proportion is not
         * less then 0 or bigger then 1
         * 
         * @param value Proportion to check
         * @param line  Line the function is called from (use __LINE__)
         */
        void sanity_check(double value, unsigned int line) const
        {
            sevirds const& res = state.current_state;

            // Can't be bigger then 1 or less then 0
            if (value < (0 - res.one_over_prec_divider) || value > (1 + res.one_over_prec_divider))
            {
                value = res.precision_divider(value);
                    AssertLong(value >= 0 && value <= 1,
                                __FILE__, line,
                                to_string(value) + " is \033[33m" + (value < 0 ? "less then zero" : "bigger then one") + "\033[31m on day " + to_string((int)simulation_clock));
            }
        }
}; //class geographical_cell{}

#endif //PANDEMIC_HOYA_2002_ZHONG_CELL_HPP
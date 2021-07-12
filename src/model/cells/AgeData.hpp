// Created by Eric - Jun/2021

#ifndef AGE_DATA_HPP
#define AGE_DATA_HPP

#include <vector>
#include "sevirds.hpp"

using namespace std;
using vecDouble = vector<double>;
using vecVecDouble = vector<vecDouble>;

static vecDouble EMPTY_VEC; // Used as a null

/**
 * Wrapper class that holds important simulation data
 * at each age segment index during local_compute()
*/
class AgeData
{
    public:
        enum PopType
        {
            NVAC,
            DOSE1,
            DOSE2
        };

    private:
        // The current age the data is referencing
        unsigned int m_currAge;

        // Proportion Vectors
        // These will be at a current age segment index so only one vector of doubles
        vecDouble& m_susceptible;
        vecDouble& m_exposed;
        vecDouble& m_infected;
        vecDouble& m_recovered;
        vecDouble  m_fatalities;

        // Config Vectors
        vecDouble const& m_incubRates;
        vecDouble const& m_recovRates;
        vecDouble const& m_fatalRates;
        vecDouble const& m_mobilityRates;
        vecDouble const& m_vacRates;
        vecDouble const& m_virulRates;
        vecDouble const& m_immuneRates;

        PopType m_popType;

    public:
        AgeData(unsigned int age, vecVecDouble& susc, vecVecDouble& exp, vecVecDouble& inf,
                vecVecDouble& rec, vecVecDouble const& incub_r, vecVecDouble const& rec_r,
                vecVecDouble const& fat_r, vecDouble const& vac_r, vecVecDouble const& mob_r,
                vecVecDouble const& vir_r, vecDouble const& immu_r, PopType type=PopType::NVAC) :
            m_currAge(age),
            m_susceptible(susc.at(age)),
            m_exposed(exp.at(age)),
            m_infected(inf.at(age)),
            m_recovered(rec.at(age)),
            m_incubRates(incub_r.at(age)),
            m_recovRates(rec_r.at(age)),
            m_fatalRates(fat_r.at(age)),
            m_mobilityRates(mob_r.at(age)),
            m_vacRates(vac_r),      // Don't .at() this one since it may be EMPTY_VEC
            m_virulRates(vir_r.at(age)),
            m_immuneRates(immu_r),   // This one too may be EMPTY_VEC
            m_popType(type)
        { }

        // Non-Vaccinated
        //  so no vaccination rate or immunity rate
        AgeData(unsigned int age, vecVecDouble& susc, vecVecDouble& exp, vecVecDouble& inf,
            vecVecDouble& rec, vecVecDouble const& incub_r, vecVecDouble const& rec_r,
            vecVecDouble const& fat_r, vecVecDouble const& mob_r, vecVecDouble const& vir_r) :
            AgeData(age, susc, exp, inf, rec, incub_r, rec_r, fat_r, EMPTY_VEC, mob_r, vir_r, EMPTY_VEC)
        { }

        vecDouble&  GetSusceptible()    { return m_susceptible; }
        vecDouble&  GetExposed()        { return m_exposed;     }
        vecDouble&  GetInfected()       { return m_infected;    }
        vecDouble&  GetRecovered()      { return m_recovered;   }
        vecDouble&  GetFatalities()     { return m_fatalities;  }

        double& GetSusceptible(int index)    { return m_susceptible.at(index);   }
        double& GetExposed(int index)        { return m_exposed.at(index);       }
        double& GetInfected(int index)       { return m_infected.at(index);      }
        double& GetRecovered(int index)      { return m_recovered.at(index);     }
        double& GetFatalities(int index)     { return m_fatalities.at(index);    }

        vecDouble const& GetIncubationRates()   { return m_incubRates;      }
        vecDouble const& GetRecoveryRates()     { return m_recovRates;      }
        vecDouble const& GetFatalityRates()     { return m_fatalRates;      }
        vecDouble const& GetMobilityRates()     { return m_mobilityRates;   }
        vecDouble const& GetVaccinationRates()  { return m_vacRates;        }
        vecDouble const& GetVirulenceRates()    { return m_virulRates;      }
        vecDouble const& GetImmunityRates()     { return m_immuneRates;     }

        double GetIncubationRate(int index)    { return m_incubRates.at(index);    }
        double GetRecoveryRate(int index)      { return m_recovRates.at(index);    }
        double GetFatalityRate(int index)      { return m_fatalRates.at(index);    }
        double GetMobilityRate(int index)      { return m_mobilityRates.at(index); }
        double GetVaccinationRate(int index)   { return m_vacRates.at(index);      }
        double GetVirulenceRate(int index)     { return m_virulRates.at(index);    }
        double GetImmunityRate(int index)      { return m_immuneRates.at(index);   }

        PopType& GetType() { return m_popType; }
        int      GetAge()  { return m_currAge; }
};

#endif // AGE_DATA_HPP
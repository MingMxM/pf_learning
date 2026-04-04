//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineH2.h"
#include "BrineFluidProperties.h"
#include "HystreSinglePhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"
#include "MathUtils.h"
#include "Conversion.h"
#include "libmesh/utility.h"

registerMooseObject("PF_learningApp", PorousFlowBrineH2);

//defineLegacyParams(PorousFlowBrineH2)

InputParameters
PorousFlowBrineH2::validParams()
{
  InputParameters params = PorousFlowFluidStateMultiComponentBase::validParams();
  params.addRequiredParam<UserObjectName>("brine_fp", "The name of the user object for brine");
  params.addRequiredParam<UserObjectName>("h2_fp", "The name of the user object for H2");
  params.addParam<unsigned int>("salt_component", 2, "The component number of salt");
  params.addClassDescription("Fluid state class for brine and H2");
  return params;
}

PorousFlowBrineH2::PorousFlowBrineH2(const InputParameters & parameters)
  : PorousFlowFluidStateMultiComponentBase(parameters),
    _salt_component(getParam<unsigned int>("salt_component")),
    _brine_fp(getUserObject<BrineFluidProperties>("brine_fp")),
    _h2o_fp(_brine_fp.getComponent(BrineFluidProperties::WATER)),
    _h2_fp(getUserObject<HystreSinglePhaseFluidProperties>("h2_fp")),
    _Mh2o(_brine_fp.molarMassH2O()),
    _invMh2o(1.0 / _Mh2o),
    _ln_invMh2o(std::log(_invMh2o)),
    _Mh2(_h2_fp.molarMass()),
    _Mnacl(_brine_fp.molarMassNaCl()),
    _Zmin(1.0e-4)
{
  // Check that the correct FluidProperties UserObjects have been provided
  if (_h2_fp.fluidName() != "hydrogen")
    paramError("h2_fp", "A valid H2 FluidProperties UserObject must be provided");

  if (_brine_fp.fluidName() != "brine")
    paramError("brine_fp", "A valid Brine FluidProperties UserObject must be provided");

  // Set the number of phases and components, and their indexes
  _num_phases = 2;
  _num_components = 3;
  _gas_phase_number = 1 - _aqueous_phase_number;
  _gas_fluid_component = 3 - _aqueous_fluid_component - _salt_component;

  // Check that _aqueous_phase_number is <= total number of phases
  if (_aqueous_phase_number >= _num_phases)
    paramError("liquid_phase_number",
               "This value is larger than the possible number of phases ",
               _num_phases);

  // Check that _aqueous_fluid_component is <= total number of fluid components
  if (_aqueous_fluid_component >= _num_components)
    paramError("liquid_fluid_component",
               "This value is larger than the possible number of fluid components",
               _num_components);

  // Check that the salt component index is not identical to the liquid fluid component
  if (_salt_component == _aqueous_fluid_component)
    paramError(
        "salt_component",
        "The value provided must be different from the value entered in liquid_fluid_component");

  // Check that _salt_component is <= total number of fluid components
  if (_salt_component >= _num_components)
    paramError("salt_component",
               "The value provided is larger than the possible number of fluid components",
               _num_components);

  _empty_fsp = FluidStateProperties(_num_components);
}

std::string
PorousFlowBrineH2::fluidStateName() const
{
  return "brine-h2";
}

void
PorousFlowBrineH2::thermophysicalProperties(Real pressure,
                                             Real temperature,
                                             Real Xnacl,
                                             Real Z,
                                             unsigned int qp,
                                             std::vector<FluidStateProperties> & fsp) const
{
  // Make AD versions of primary variables then call AD thermophysicalProperties()
  ADReal p = pressure;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);
  ADReal T = temperature;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);
  ADReal Zh2 = Z;
  Moose::derivInsert(Zh2.derivatives(), _Zidx, 1.0);
  ADReal X = Xnacl;
  Moose::derivInsert(X.derivatives(), _Xidx, 1.0);

  thermophysicalProperties(p, T, X, Zh2, qp, fsp);
}

void
PorousFlowBrineH2::thermophysicalProperties(const ADReal & pressure,
                                             const ADReal & temperature,
                                             const ADReal & Xnacl,
                                             const ADReal & Z,
                                             unsigned int qp,
                                             std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(pressure.value(), temperature.value());

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(pressure, temperature, Xnacl, Z, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      // Set the gas saturations
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(pressure, temperature, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      // Calculate the liquid properties
      const ADReal liquid_pressure = pressure - _pc.capillaryPressure(1.0, qp);
      liquidProperties(liquid_pressure, temperature, Xnacl, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas and liquid properties in the two phase region
      twoPhaseProperties(pressure, temperature, Xnacl, Z, qp, fsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;

  // Save pressures to FluidStateProperties object
  gas.pressure = pressure;
  liquid.pressure = pressure - _pc.capillaryPressure(liquid.saturation, qp);
}

void
PorousFlowBrineH2::massFractions(const ADReal & pressure,
                                  const ADReal & temperature,
                                  const ADReal & Xnacl,
                                  const ADReal & Z,
                                  FluidStatePhaseEnum & phase_state,
                                  std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  ADReal Xh2 = 0.0;
  ADReal Yh2o = 0.0;
  ADReal Yh2 = 0.0;

  // If the amount of H2 is less than the smallest solubility, then all H2 will
  // be dissolved, and the equilibrium mass fractions do not need to be computed
  if (Z.value() < _Zmin)
    phase_state = FluidStatePhaseEnum::LIQUID;

  else
  {
    // Equilibrium mass fraction of H2 in liquid and H2O in gas phases
    equilibriumMassFractions(pressure, temperature, Xnacl, Xh2, Yh2o);

    Yh2 = 1.0 - Yh2o;

    // Determine which phases are present based on the value of z
    phaseState(Z.value(), Xh2.value(), Yh2.value(), phase_state);
  }

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z
  ADReal Xh2o = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xh2 = Z;
      Yh2 = 0.0;
      Xh2o = 1.0 - Z;
      Yh2o = 0.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xh2 = 0.0;
      Yh2 = Z;
      Yh2o = 1.0 - Z;
      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      Xh2o = 1.0 - Xh2;
      break;
    }
  }

  // Save the mass fractions in the FluidStateProperties object
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xh2;
  liquid.mass_fraction[_salt_component] = Xnacl;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yh2;
}

void
PorousFlowBrineH2::gasProperties(const ADReal & pressure,
                                  const ADReal & temperature,
                                  std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Gas density, viscosity and enthalpy are approximated with pure H2 - no correction due
  // to the small amount of water vapor is made
  ADReal h2_density, h2_viscosity;
  _h2_fp.rho_mu_from_p_T(pressure, temperature, h2_density, h2_viscosity);

  ADReal h2_enthalpy = _h2_fp.h_from_p_T(pressure, temperature);

  // Save the values to the FluidStateProperties object. Note that derivatives wrt z are 0
  gas.density = h2_density;
  gas.viscosity = h2_viscosity;
  gas.enthalpy = h2_enthalpy;

  mooseAssert(gas.density.value() > 0.0, "Gas density must be greater than zero");
  gas.internal_energy = gas.enthalpy - pressure / gas.density;
}

void
PorousFlowBrineH2::liquidProperties(const ADReal & pressure,
                                     const ADReal & temperature,
                                     const ADReal & Xnacl,
                                     std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved H2
  const ADReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of H2 in liquid phase
  const ADReal Xh2 = liquid.mass_fraction[_gas_fluid_component];

  // The partial density of H2
  const ADReal h2_partial_density = _Mh2 / molarVolume(pressure, temperature);

  // The liquid density
  const ADReal liquid_density = 1.0 / (Xh2 / h2_partial_density + (1.0 - Xh2) / brine_density);

  // Assume that liquid viscosity is just the brine viscosity
  const ADReal liquid_viscosity = _brine_fp.mu_from_p_T_X(pressure, temperature, Xnacl);

  // Liquid enthalpy (including contribution due to the enthalpy of dissolution)
  const ADReal brine_enthalpy = _brine_fp.h_from_p_T_X(pressure, temperature, Xnacl);

  // Enthalpy of H2
  const ADReal h2_enthalpy = _h2_fp.h_from_p_T(pressure, temperature);

  // Enthalpy of liquid
  const ADReal liquid_enthalpy = (1.0 - Xh2) * brine_enthalpy + Xh2 * (h2_enthalpy);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;
  liquid.enthalpy = liquid_enthalpy;

  mooseAssert(liquid.density.value() > 0.0, "Liquid density must be greater than zero");
  liquid.internal_energy = liquid.enthalpy - pressure / liquid.density;
}

ADReal
PorousFlowBrineH2::saturation(const ADReal & pressure,
                               const ADReal & temperature,
                               const ADReal & Xnacl,
                               const ADReal & Z,
                               std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
  auto & liquid = fsp[_aqueous_fluid_component];
  // auto & liquid = fsp[_aqueous_phase_number];

  // Approximate liquid density as saturation isn't known yet, by using the gas
  // pressure rather than the liquid pressure. This does result in a small error
  // in the calculated saturation, but this is below the error associated with
  // the correlations. A more accurate saturation could be found iteraviely,
  // at the cost of increased computational expense

  // Gas density
  const ADReal gas_density = _h2_fp.rho_from_p_T(pressure, temperature);

  // Approximate liquid density as saturation isn't known yet
  const ADReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of H2 in liquid phase
  const ADReal Xh2 = liquid.mass_fraction[_gas_fluid_component];

  // The liquid density
  const ADReal liquid_density = brine_density;

  const ADReal Yh2 = gas.mass_fraction[_gas_fluid_component];

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  const ADReal K0 = Yh2 / Xh2;
  const ADReal K1 = (1.0 - Yh2) / (1.0 - Xh2);
  const ADReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const ADReal saturation = vapor_mass_fraction * liquid_density /
                            (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  return saturation;
}

void
PorousFlowBrineH2::twoPhaseProperties(const ADReal & pressure,
                                       const ADReal & temperature,
                                       const ADReal & Xnacl,
                                       const ADReal & Z,
                                       unsigned int qp,
                                       std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];

  // Calculate all of the gas phase properties, as these don't depend on saturation
  gasProperties(pressure, temperature, fsp);

  // The gas saturation in the two phase case
  gas.saturation = saturation(pressure, temperature, Xnacl, Z, fsp);

  // The liquid pressure and properties can now be calculated
  const ADReal liquid_pressure = pressure - _pc.capillaryPressure(1.0 - gas.saturation, qp);
  liquidProperties(liquid_pressure, temperature, Xnacl, fsp);
}

void
PorousFlowBrineH2::equilibriumMassFractions(const ADReal & pressure,
                                             const ADReal & temperature,
                                             const ADReal & Xnacl,
                                             ADReal & Xh2,
                                             ADReal & Yh2o) const
{
  // Mole fractions at equilibrium
  ADReal xh2, yh2o;
  equilibriumMoleFractions(pressure, temperature, Xnacl, xh2, yh2o);

  // The mass fraction of H2O in gas (assume no salt in gas phase) and derivatives
  // wrt p, T, and X
  Yh2o = yh2o * _Mh2o / (yh2o * _Mh2o + (1.0 - yh2o) * _Mh2);

  // NaCl molality (mol/kg)
  const ADReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // The molality of H2 in 1kg of H2O
  const ADReal mh2 = xh2 * (2.0 * mnacl + _invMh2o) / (1.0 - xh2);
  // The mass fraction of H2 in brine is then
  const ADReal denominator = (1.0 + mnacl * _Mnacl + mh2 * _Mh2);
  Xh2 = mh2 * _Mh2 / denominator;
}

void
PorousFlowBrineH2::equilibriumMoleFractions(const ADReal & pressure,
                                             const ADReal & temperature,
                                             const ADReal & Xnacl,
                                             ADReal & xh2,
                                             ADReal & yh2o) const
{
  // To calculate equilibrium mole fractions, need NaCl molality
  const ADReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // Molality of H2 in liquid phas
  const ADReal mh2 = equilibriumMolality(pressure, temperature, mnacl);

  // Mole fraction of H2 in liquid phase, this part may have problem, as the corect equation
  // should be xh2 = mh2 / (_invMH2O + mh2 + 2*mnacl)
  xh2 = mh2 / (_invMh2o + mh2);

  // Mole fraction of H2O in gas phase
  yh2o = this->yh2o(pressure, temperature);

}

ADReal
PorousFlowBrineH2::equilibriumMolality(const ADReal & pressure,
                                       const ADReal & temperature,
                                       const ADReal & mnacl) const
{
  using std::log, std::exp;

  const ADReal lnyH2 = log(1.0 - yh2o(pressure, temperature));
  const ADReal lnphi = log(fugacityCoefficient(pressure, temperature));
  const ADReal lnKh = log(henryConstant(temperature));
  const ADReal PF = poyntingFactor(pressure, temperature);
  const ADReal lngamma = log(activityCoefficient(temperature, mnacl));

  const ADReal lnmH2 =
      lnyH2 + log(pressure * 1.0e-6) + lnphi - lnKh - PF - lngamma + _ln_invMh2o;

  return exp(lnmH2);
}

ADReal
PorousFlowBrineH2::yh2o(const ADReal & pressure, const ADReal & temperature) const
{
  // mooseAssert(pressure > 0.0, "PorousFlowBrineH2::yH2O(): pressure must be greater than zero");

  return _h2o_fp.vaporPressure(temperature) / pressure;
}

ADReal
PorousFlowBrineH2::fugacityCoefficient(const ADReal & pressure,
                                       const ADReal & temperature) const
{
  return _h2_fp.psi_from_p_T(pressure, temperature);
}

ADReal
PorousFlowBrineH2::henryConstant(const ADReal & temperature) const
{
  using std::exp;
  // Henry's constant for dissolution in water
  const ADReal lnKh_H2O = _a[0] * Utility::pow<2>(temperature) + _a[1] * temperature + _a[2] +
                            _a[3] / temperature + _a[4] / Utility::pow<2>(temperature);

  // Note: ln(Kh) provided in Li et al
  return exp(lnKh_H2O);
}

ADReal
PorousFlowBrineH2::poyntingFactor(const ADReal & pressure, const ADReal & temperature) const
{
  // Pressure in MPa
  const ADReal p = pressure * 1.0e-6;

  const ADReal PF = _b[0] * p / temperature + _b[1] * p + _b[2] * p * temperature +
                      _b[3] * Utility::pow<2>(p) / temperature;

  return PF;
}

ADReal
PorousFlowBrineH2::activityCoefficient(const ADReal & temperature, const ADReal & mnacl) const
{
  using std::exp;
  // Note: ln(gamma) provided in Li et al
  return exp((0.64485 - 0.00142 * temperature) * mnacl);
}

ADReal
PorousFlowBrineH2::molarVolume(const ADReal & pressure, const ADReal & temperature) const
{
  // Pressure in MPa
  const ADReal p = pressure * 1.0e-6;
  const ADReal V =
      51.1904 - 0.208062 * temperature + 3.4427e-4 * Utility::pow<2>(temperature) - 0.022 * p;

  // V is in cm^3/mol, return m^3/mol
  return V * 1.0e-6;
}

void
PorousFlowBrineH2::checkVariables(Real pressure, Real temperature) const
{
  // The calculation of mass fractions is valid from 0C <= T <= 100C, and
  // pressure less than 60 MPa
  if (temperature < 273.15 || temperature > 373.15)
    mooseException(name() + ": temperature " + Moose::stringify(temperature) +
                   " is outside range 273.15 K <= T <= 373.15 K");

  if (pressure > 60.0e7)
    mooseException(name() + ": pressure " + Moose::stringify(pressure) +
                   " must be less than 60 MPa");
}

Real
PorousFlowBrineH2::totalMassFraction(
    Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const
{
  // Check whether the input pressure and temperature are within the region of validity
  checkVariables(pressure, temperature);

  // As we do not require derivatives, we can simply ignore their initialisation
  const ADReal p = pressure;
  const ADReal T = temperature;
  const ADReal X = Xnacl;

  // FluidStateProperties data structure
  std::vector<FluidStateProperties> fsp(_num_phases, FluidStateProperties(_num_components));
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Calculate equilibrium mass fractions in the two-phase state
  ADReal Xh2, Yh2o;
  equilibriumMassFractions(p, T, X, Xh2, Yh2o);

  // Save the mass fractions in the FluidStateMassFractions object
  const ADReal Yh2 = 1.0 - Yh2o;
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - Xh2;
  liquid.mass_fraction[_gas_fluid_component] = Xh2;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yh2;

  // Gas properties
  gasProperties(pressure, temperature, fsp);

  // Liquid properties
  const ADReal liquid_saturation = 1.0 - saturation;
  const ADReal liquid_pressure = p - _pc.capillaryPressure(liquid_saturation, qp);
  liquidProperties(liquid_pressure, T, X, fsp);

  // The total mass fraction of ncg (z) can now be calculated
  const ADReal Z = (saturation * gas.density * Yh2 + liquid_saturation * liquid.density * Xh2) /
                     (saturation * gas.density + liquid_saturation * liquid.density);

  return Z.value();
}

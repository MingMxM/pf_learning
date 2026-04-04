//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowFluidStateMultiComponentBase.h"

class BrineFluidProperties;
class HystreSinglePhaseFluidProperties;
class SinglePhaseFluidProperties;
class Water97FluidProperties;
class PorousFlowBrineH2;

/**
 * Specialized class for brine and H2 base on
 * Li, Beyer, and Bauer, A unified phase equilibrium model for hydrogen solubility
 * and solution density. International Journal of hydrogen energy, 43, 512-529 (2018).
 * 
 * This EOS is valid for temperaute up to 100C (373.15K), pressure 
 * between 1 MPa and 50 Mpa, and salinity up to 5 mol/kg NaCl
 * 
 * Notation convention
 * Throughout this class, both mole fractions and mass fractions will be used.
 * The following notation will be used:
 * yk: mole fraction of component k in the gas phase
 * xk: mole fraction of component k in the liquid phase
 * Yk: mass fraction of component k in the gas phase
 * Xk: mass fraction of component k in the liquid phase
 */
class PorousFlowBrineH2 : public PorousFlowFluidStateMultiComponentBase
{
public:
  static InputParameters validParams();

  PorousFlowBrineH2(const InputParameters & parameters);

  virtual std::string fluidStateName() const override;

  void thermophysicalProperties(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real Z,
                                unsigned int qp,
                                std::vector<FluidStateProperties> & fsp) const override;

  void thermophysicalProperties(const ADReal & pressure,
                                const ADReal & temperature,
                                const ADReal & Xnacl,
                                const ADReal & Z,
                                unsigned int qp,
                                std::vector<FluidStateProperties> & fsp) const override;

  /**
   * Mole fractions of H2 in brine and water vapor in H2 at equilibrium
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] xh2 mole fraction of H2 in liquid
   * @param[out] yh2o mole fraction of H2O in gas
   */
  void equilibriumMoleFractions(const ADReal & pressure,
                                const ADReal & temperature,
                                const ADReal & Xnacl,
                                ADReal & xh2,
                                ADReal & yh2o) const;

  /**
   * Mass fractions of H2 in brine and water vapor in H2 at equilibrium
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] Xh2 mass fraction of h2 in liquid (kg/kg)
   * @param[out] Yh2o mass fraction of H2O in gas (kg/kg)
   */
  void equilibriumMassFractions(const ADReal & pressure,
                                const ADReal & temperature,
                                const ADReal & Xnacl,
                                ADReal & Xh2,
                                ADReal & Yh2o) const;

  /**
   * Mass fractions of H2 and H2O in both phases, as well as derivatives wrt
   * PorousFlow variables. Values depend on the phase state (liquid, gas or two phase)
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of H2 component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateProperties data structure
   */
  void massFractions(const ADReal & pressure,
                     const ADReal & temperature,
                     const ADReal & Xnacl,
                     const ADReal & Z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the gaseous state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param[out] FluidStateProperties data structure
   */
  void gasProperties(const ADReal & pressure,
                     const ADReal & temperature,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the liquid state
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateProperties data structure
   */
  void liquidProperties(const ADReal & pressure,
                        const ADReal & temperature,
                        const ADReal & Xnacl,
                        std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas saturation in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of H2 component
   * @param FluidStateProperties data structure
   * @return gas saturation (-)
   */
  ADReal saturation(const ADReal & pressure,
                    const ADReal & temperature,
                    const ADReal & Xnacl,
                    const ADReal & Z,
                    std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas and liquid properties in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of NCG component
   * @param qp quadpoint for capillary presssure
   * @param[out] FluidStateProperties data structure
   */
  void twoPhaseProperties(const ADReal & pressure,
                          const ADReal & temperature,
                          const ADReal & Xnacl,
                          const ADReal & Z,
                          unsigned int qp,
                          std::vector<FluidStateProperties> & fsp) const;

   /**
   * Henry's constant of dissolution of gas phase H2 in brine. (Eq.(13) of Li et al)
   *
   * @param temperature fluid temperature (K)
   * @return Henry's constant (Pa)
   */
  ADReal henryConstant(const ADReal & temperature) const;

  /**
   * Apparent molar volume of H2 in H2O (Eq. (14) of Li et al)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @return molar volume (m^3/mol)
   */
  ADReal molarVolume(const ADReal & pressure, const ADReal & temperature) const;

  /**
   * Poynting factor (Eq. (16) of Li et al)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @return Poynting factor
   */
  ADReal poyntingFactor(const ADReal & pressure, const ADReal & temperature) const;

  /**
   * Activity coefficient of H2 in brine (Eq. (17) of Li et al)
   * Note: this is gamma, not ln(gamma) as in Eq. (17)
   *
   * @param temperature fluid temperature (K)
   * @param mnacl molality of NaCl in solution (mol/kg)
   * @return gamma the activity coefficient
   */
  ADReal activityCoefficient(const ADReal & temperature, const ADReal & mnacl) const;

  /**
   * Mole fraction of H2O in gas phase (yH2O)  (Dalton's law)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @return yH2O the mole fraction of H2O in gas phase
   */
  ADReal yh2o(const ADReal & pressure, const ADReal & temperature) const;

  /**
   * Fugactiy coefficient of H2  (Eq. (8) of Li et al)
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @return fugacity coefficient (-)
   */
  ADReal fugacityCoefficient(const ADReal & pressure, const ADReal & temperature) const;

  /**
   * Molality (mol/kg) of H2 in liquid phase (Eq. (6) of Li et al)
   * Note: this is molality, not ln(molality) as in Eq. (6)
   * Note: There is a typo in Eq. (6) from Li et al. The last term should be '+ 4.0166'
   * instead of '- 4.0166'
   *
   * @param pressure fluid pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param mnacl molality of NaCl in solution (mol/kg)
   * @return molality of H2 in liquid phase (mol/kg)
   */
  ADReal equilibriumMolality(const ADReal & pressure,
                               const ADReal & temperature,
                               const ADReal & mnacl) const;

  /**
   * The index of the salt component
   * @return salt component number
   */
  unsigned int saltComponentIndex() const { return _salt_component; };

  virtual Real totalMassFraction(
      Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const override;

protected:
  /// Check the input variables
  virtual void checkVariables(Real pressure, Real temperature) const;
  /// Salt component index
  const unsigned int _salt_component;
  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for the H2O
  const SinglePhaseFluidProperties & _h2o_fp;
  /// Fluid properties UserObject for the H2
  const HystreSinglePhaseFluidProperties & _h2_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Inverse of molar mass of H2O (mol/kg)
  const Real _invMh2o;
  /// Natural logarithm of inverse of molar mass of H2O (mol/kg)
  const Real _ln_invMh2o;
  /// Molar mass of H2 (kg/mol)
  const Real _Mh2;
  /// Molar mass of NaCL
  const Real _Mnacl;
  /// Minimum Z - below this value all H2 will be dissolved. This reduces the
  /// computational burden when small values of Z are present
  const Real _Zmin;
  /// Henry constant coefficients
  const std::array<Real, 5> _a{{2.68721e-5, -0.05121, 33.55196, -3411.0432, -31258.74683}};
  /// Poynting factor coefficients
  const std::array<Real, 4> _b{{6.156755, -2.502396e-2, 4.140593e-5, -1.322988e-3}};
};

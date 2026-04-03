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
   * @param[out] xh2 mole fraction of CO2 in liquid
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
   * @param[out] Xh2 mass fraction of CO2 in liquid (kg/kg)
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
   * @param Z total mass fraction of CO2 component
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
   * The index of the salt component
   * @return salt component number
   */
  unsigned int saltComponentIndex() const { return _salt_component; };



  /**
   * Enthalpy of dissolution of gas phase CO2 in brine calculated using Henry's constant
   * From Himmelblau, Partial molal heats and entropies of solution for gases dissolved
   * in water from the freezing to the near critical point, J. Phys. Chem. 63 (1959).
   * Correction due to salinity from Battistelli et al, A fluid property module for the
   * TOUGH2 simulator for saline brines with non-condensible gas, Proc. Eighteenth Workshop
   * on Geothermal Reservoir Engineering (1993).
   *
   * @param temperature fluid temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @return enthalpy of dissolution (J/kg)
   */
  ADReal enthalpyOfDissolutionGas(const ADReal & temperature, const ADReal & Xnacl) const;

  /**
   * Enthalpy of dissolution of CO2 in brine calculated using linear fit to model of
   * Duan and Sun, An improved model calculating CO2 solubility in pure water and aqueous NaCl
   * solutions from 273 to 533 K and from 0 to 2000 bar, Chemical geology, 193, 257--271 (2003).
   *
   * In the region of interest, the more complicated model given in Eq. (8) of Duan and Sun
   * is well approximated by a simple linear fit (R^2 = 0.9997).
   *
   * Note: as the effect of salt mass fraction is small, it is not included in this function.
   *
   * @param temperature fluid temperature (K)
   * @return enthalpy of dissolution (J/kg)
   */
  ADReal enthalpyOfDissolution(const ADReal & temperature) const;

  /**
   * Mole fractions of CO2 in brine and water vapor in CO2 at equilibrium in the low
   * temperature regime (T <= 99C).
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] xco2 mole fraction of CO2 in liquid
   * @param[out] yh2o mass fraction of mole in gas
   */
  void equilibriumMoleFractionsLowTemp(const ADReal & pressure,
                                       const ADReal & temperature,
                                       const ADReal & Xnacl,
                                       ADReal & xco2,
                                       ADReal & yh2o) const;

  /**
   * Function to solve for yh2o and xco2 iteratively in the elevated temperature regime (T > 100C)
   *
   * @param pressure gas pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param co2_density CO2 density (kg/m^3)
   * @param[out] xco2 mole fraction of CO2 in liquid phase (-)
   * @param[out] yh2o mole fraction of H2O in gas phase (-)
   */
  void solveEquilibriumMoleFractionHighTemp(Real pressure,
                                            Real temperature,
                                            Real Xnacl,
                                            Real co2_density,
                                            Real & xco2,
                                            Real & yh2o) const;

protected:
  /// Check the input variables
  virtual void checkVariables(Real pressure, Real temperature) const;

  /**
   * Binary interaction coefficient of H2, N2, and CH4 (Eq.(26) of Li et al)
   * 
   * @param temperature fluid temperature (K)
   * @param a array of constans for evaluation of Eq.(26)
   * @return binary interaction coefficient (m^3/mol)
   */
  ADReal binaryInteractionCoeff(const ADReal & temperature,
                                const std::array<Real, 4> & a) const;

  /// Salt component index
  const unsigned int _salt_component;
  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for the H2O
  const SinglePhaseFluidProperties & _H2O_fp;
  /// Fluid properties UserObject for the H2
  const HystreSinglePhaseFluidProperties & _h2_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Inverse of molar mass of H2O (mol/kg)
  const Real _invMh2o;
  
  /// Molar mass of CO2 (kg/mol)
  const Real _Mco2;
  /// Molar mass of NaCL
  const Real _Mnacl;
  /// Molar gas constant in bar cm^3 /(K mol)
  const Real _Rbar;
  /// Temperature below which the Spycher, Pruess & Ennis-King (2003) model is used (K)
  const Real _Tlower;
  /// Temperature above which the Spycher & Pruess (2010) model is used (K)
  const Real _Tupper;
  /// Minimum Z - below this value all CO2 will be dissolved. This reduces the
  /// computational burden when small values of Z are present
  const Real _Zmin;
  /// Henry's coefficeients for CO2
  const std::vector<Real> _co2_henry;
};

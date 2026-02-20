//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "PF_learningTestApp.h"
#include "PF_learningApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
PF_learningTestApp::validParams()
{
  InputParameters params = PF_learningApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

PF_learningTestApp::PF_learningTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  PF_learningTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

PF_learningTestApp::~PF_learningTestApp() {}

void
PF_learningTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  PF_learningApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"PF_learningTestApp"});
    Registry::registerActionsTo(af, {"PF_learningTestApp"});
  }
}

void
PF_learningTestApp::registerApps()
{
  registerApp(PF_learningApp);
  registerApp(PF_learningTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
PF_learningTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PF_learningTestApp::registerAll(f, af, s);
}
extern "C" void
PF_learningTestApp__registerApps()
{
  PF_learningTestApp::registerApps();
}

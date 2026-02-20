#include "PF_learningApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
PF_learningApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

PF_learningApp::PF_learningApp(const InputParameters & parameters) : MooseApp(parameters)
{
  PF_learningApp::registerAll(_factory, _action_factory, _syntax);
}

PF_learningApp::~PF_learningApp() {}

void
PF_learningApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<PF_learningApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"PF_learningApp"});
  Registry::registerActionsTo(af, {"PF_learningApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
PF_learningApp::registerApps()
{
  registerApp(PF_learningApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
PF_learningApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PF_learningApp::registerAll(f, af, s);
}
extern "C" void
PF_learningApp__registerApps()
{
  PF_learningApp::registerApps();
}

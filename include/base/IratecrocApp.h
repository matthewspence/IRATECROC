#ifndef IRATECROCAPP_H
#define IRATECROCAPP_H

#include "MooseApp.h"

class IratecrocApp;

template<>
InputParameters validParams<IratecrocApp>();

class IratecrocApp : public MooseApp
{
public:
  IratecrocApp(InputParameters parameters);
  virtual ~IratecrocApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* IRATECROCAPP_H */

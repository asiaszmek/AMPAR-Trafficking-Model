# AMPAR-Trafficking-Model

AMPAR-Trafficking-Model contains the Python code used in the paper "The biophysical basis underlying the maintenance of early phase long-term potentiation". 

* A package "ampartrafficking" is provided containing four modules:
  ** rate_model: Contains functions and classes for the rate model of AMPAR-trafficking (inlcuding the mean-field approximation of cooperative binding/unbinding rates).
  ** stochastic_model: Contains functions and classes for the stochastic, spatial model of cooperative receptor binding.
  ** parameter_sampling: Contains functions for random parameter sampling and subsequent simulation of the model. Depends on rate_model.
  ** frap: Contains functions to carry out fluorescence recovery after photobleaching (FRAP) simulations. Depends on stochastic_model.
* The folders "Mean-Field-Model" and "Stochastic-Model" further contain python scripts to reproduce the figures shown in the paper. Each script is named after the respective figure as appearing in the Paper.
* Documentation of the ampartrafficking package and its modules is provided in the docs folder (docs/build/index.html).

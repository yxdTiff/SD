import os

import pyomo_sdlib

import baa99_util
import baa99_pyomo

reference_model = baa99_pyomo.model
reference_model.name = 'Scenario1'
scenario_tree_model = \
    baa99_util.generate_scenario_tree_model(1)

pyomo_sdlib.setup('baa99_pyomo',
                  reference_model,
                  scenario_tree_model)

import baa99_exotic_pyomo
reference_model = baa99_exotic_pyomo.model
reference_model.name = 'Scenario1'
scenario_tree_model = \
    baa99_util.generate_scenario_tree_model(1)

pyomo_sdlib.setup('baa99_exotic_pyomo',
                  reference_model,
                  scenario_tree_model)

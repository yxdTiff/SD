import os
import operator

from pyutilib.services import TempfileManager

from pyomo.pysp.scenariotree import ScenarioTree
from pyomo.core import *
from pyomo.core.base.piecewise import _PiecewiseData
from pyomo.repn import LinearCanonicalRepn

from six import iteritems

# NOTES:
#  - Do we handle SOS constraints in the LP file?
#  - Constants in the objective?

def map_constraint_stages(reference_model, scenario_tree, LP_symbol_map):

    reference_scenario = scenario_tree.get_scenario(reference_model.name)

    assert len(scenario_tree._stages) == 2
    stage1 = scenario_tree._stages[0]
    stage2 = scenario_tree._stages[1]

    StageToConstraintMap = {}
    StageToConstraintMap[stage1._name] = []
    StageToConstraintMap[stage2._name] = []

    LP_byObject = LP_symbol_map.byObject
    # deal with the fact that the LP writer prepends constraint
    # names with things like 'c_e_', 'c_l_', etc depending on the
    # constraint bound type and will even split a constraint into
    # two constraints if it has two bounds
    LP_reverse_alias = dict()
    for symbol in LP_symbol_map.bySymbol:
        LP_reverse_alias[symbol] = []
    for alias, obj_weakref in iteritems(LP_symbol_map.aliases):
        LP_reverse_alias[LP_byObject[id(obj_weakref())]].append(alias)
    for block in reference_model.all_blocks(active=True):
        canonical_repn = getattr(block,"canonical_repn",None)
        if canonical_repn is None:
            raise ValueError("Unable to find canonical_repn ComponentMap "
                             "on block %s" % (block.cname(True)))
        piecewise_stage = None
        if isinstance(block, (Piecewise, _PiecewiseData)):
            piecewise_stage = stage1
            for vardata in block.referenced_variables():
                variable_node = scenario_tree.variableNode(vardata, reference_model)
                if variable_node._stage == stage2:
                    piecewise_stage = stage2
                else:
                    assert variable_node._stage == stage1

        for constraint_data in block.active_component_data.itervalues(SOSConstraint, descend_into=False):
            raise TypeError("SOSConstraints are not handled by the SD interface: %s"
                            % (constraint_data.cname(True)))

        for constraint_data in block.active_component_data.itervalues(Constraint, descend_into=False):
            LP_name = LP_byObject[id(constraint_data)]
            # if it is a range constraint this will account for
            # that fact and hold and alias for each bound
            LP_aliases = LP_reverse_alias[LP_name]
            assert len(LP_aliases) > 0
            if piecewise_stage is None:
                constraint_node = reference_scenario.constraintNode(constraint_data,
                                                                    repn=canonical_repn,
                                                                    instance=reference_model)
                constraint_stage = constraint_node._stage
            else:
                constraint_stage = piecewise_stage

            StageToConstraintMap[constraint_stage._name].append((LP_aliases, constraint_data))

    assert sorted(StageToConstraintMap.keys()) == \
        sorted([stage1._name, stage2._name])

    # sort each by name
    for key in StageToConstraintMap:
        StageToConstraintMap[key].sort(key=operator.itemgetter(0))

    return StageToConstraintMap

def map_variable_stages(reference_model, scenario_tree, LP_symbol_map):

    FirstStageVars = {}
    SecondStageVars = {}

    all_vars_cnt = 0
    for block in reference_model.all_blocks(active=True):
        all_vars_cnt += len(list(block.all_component_data.itervalues(Var,
                                                                     descend_into=False)))

    rootnode = scenario_tree.findRootNode()
    assert len(scenario_tree._stages) == 2
    stageone = scenario_tree._stages[0]
    stagetwo = scenario_tree._stages[1]
    anode = stagetwo._tree_nodes[0]
    firststage_blended_variables = rootnode._standard_variable_ids
    LP_byObject = LP_symbol_map.byObject
    all_vars_on_tree = []
    for scenario_tree_id, vardata in \
          iteritems(reference_model.\
          _ScenarioTreeSymbolMap.bySymbol):
        try:
            LP_name = LP_byObject[id(vardata)]
        except:
            print(("FAILED ON VAR DATA= "+vardata.cname(True)))
            foobar
        if LP_name == "RHS":
            raise RuntimeError("Congratulations! You have hit an edge case. The SD input "
                               "format forbids variables having name 'RHS'. Please rename it")
        if scenario_tree_id in firststage_blended_variables:
            FirstStageVars[LP_name] = (vardata, scenario_tree_id)
        elif (scenario_tree_id in rootnode._derived_variable_ids) or \
             (scenario_tree_id in anode._variable_ids):
            SecondStageVars[LP_name] = (vardata, scenario_tree_id)
        else:
            # More than two stages?
            assert False
        all_vars_on_tree.append(LP_name)

    for stage in scenario_tree._stages:
        cost_variable_name, cost_variable_index = \
            stage._cost_variable
        stage_cost_component = \
            reference_model.\
            find_component(cost_variable_name)
        if stage_cost_component.type() is not Expression:
            raise RuntimeError("The StageCostVariables must be declared "
                               "as Expression objects when using this tool")
            """
            stage_cost_vardata = stage_cost_component[cost_variable_index]
            LP_name = LP_byObject[id(stage_cost_vardata)]
            assert LP_name not in FirstStageVars
            if LP_name not in all_vars_on_tree:
                # This would imply the user has declared something as
                # a StageCostVariable and a StageVariable (and they need
                # to fix)
                assert LP_name not in SecondStageVars
                SecondStageVars[LP_name] = (stage_cost_vardata, None)
                all_vars_on_tree.append(LP_name)
            """
    # Make sure every variable on the model has been
    # declared on the scenario tree
    if len(all_vars_on_tree) != all_vars_cnt:
        print("**** THERE IS A PROBLEM ****")
        print("Not all model variables are on the scenario tree. Investigating...")
        all_vars = set()
        for block in reference_model.all_blocks(active=True):
            all_vars.update(vardata.cname(True) \
                            for vardata in block.all_component_data.itervalues(Var,
                                                                               descend_into=False))
        tree_vars = set()
        for scenario_tree_id, vardata in \
            iteritems(reference_model.\
                      _ScenarioTreeSymbolMap.bySymbol):
            tree_vars.add(vardata.cname(True))
        cost_vars = set()
        for stage in scenario_tree._stages:
            cost_variable_name, cost_variable_index = \
                stage._cost_variable
            stage_cost_component = \
                reference_model.\
                find_component(cost_variable_name)
            if stage_cost_component.type() is not Expression:
                cost_vars.add(stage_cost_component[cost_variable_index].cname(True))
        print(("Number of Scenario Tree Variables (found ddsip LP file): "+str(len(tree_vars))))
        print(("Number of Scenario Tree Cost Variables (found ddsip LP file): "+str(len(cost_vars))))
        print(("Number of Variables Found on Model: "+str(len(all_vars))))
        print(("Variables Missing from Scenario Tree (or LP file):"+str(all_vars-tree_vars-cost_vars)))

    # A necessary but not sufficient sanity check to make sure the
    # second stage variable sets are the same for all
    # scenarios. This is not required by pysp, but I think this
    # assumption is made in the rest of the code here
    for tree_node in stagetwo._tree_nodes:
        assert len(anode._variable_ids) == \
            len(tree_node._variable_ids)

    assert len(scenario_tree._stages) == 2

    StageToVariableMap = {}
    StageToVariableMap[stageone._name] = \
        [(LP_name, FirstStageVars[LP_name][0], FirstStageVars[LP_name][1])
         for LP_name in sorted(FirstStageVars)]
    StageToVariableMap[stagetwo._name] = \
        [(LP_name, SecondStageVars[LP_name][0], SecondStageVars[LP_name][1])
         for LP_name in sorted(SecondStageVars)]

    return StageToVariableMap

def setup(output_prefix, reference_model, scenario_tree_model):

    import pyomo.environ
    import pyomo.repn.plugins.cpxlp

    stoch = reference_model.find_component('PySP_StochasticParameters')
    for param_data in stoch:
        if not param_data.parent_component()._mutable:
            raise RuntimeError("Stochastic parameters must be mutable")
        if value(param_data) == 0:
            raise RuntimeError("Stochastic parameters should be initialized "
                               "with nonzero values to avoid issues due to sparsity "
                               "of Pyomo expressions. Please update the value of %s"
                               % (param_data.cname(True)))
    try:
        os.mkdir(output_prefix)
    except OSError:
        pass

    scenario_tree = ScenarioTree(scenariotreeinstance=scenario_tree_model,
                                 scenariobundlelist=[reference_model.name])
    scenario_tree.linkInInstances({reference_model.name: reference_model})

    reference_scenario = scenario_tree.get_scenario(reference_model.name)
    reference_model.preprocess()

    with pyomo.repn.plugins.cpxlp.ProblemWriter_cpxlp() as writer:
        lp_filename = os.path.join(output_prefix, output_prefix+".setup.lp")
        io_options = {'symbolic_solver_labels': True}
        output_filename, LP_symbol_map = writer(reference_model,
                                                lp_filename,
                                                lambda x: True,
                                                io_options)
        assert output_filename == lp_filename

    StageToVariableMap = \
        map_variable_stages(reference_model, scenario_tree, LP_symbol_map)

    StageToConstraintMap = \
        map_constraint_stages(reference_model, scenario_tree, LP_symbol_map)

    assert len(scenario_tree._stages) == 2
    firststage = scenario_tree._stages[0]
    secondstage = scenario_tree._stages[1]

    #
    # Make sure the objective references all first stage variables.
    # We do this by directly modifying the canonical_repn of the objective
    # which the LP writer will reference next time we call it. In addition,
    # make that that the first second-stage variable in our column
    # ordering also appears in the objective so that ONE_VAR_CONSTANT
    # does not get identified as the first second-stage variable.
    # ** Just do NOT preprocess again until we call the LP writer **
    #
    canonical_repn = reference_model.canonical_repn
    obj_repn = canonical_repn[reference_scenario._instance_objective]
    first_stage_varname_list = [item[0] for item in StageToVariableMap[firststage._name]]
    if isinstance(obj_repn, LinearCanonicalRepn) and \
       (obj_repn.linear is not None):
        referenced_var_names = [LP_symbol_map.byObject[id(vardata)]
                                for vardata in obj_repn.variables]
        update_vars = []
        # add the first-stage variables (if not present)
        for LP_name in first_stage_varname_list:
            if LP_name not in referenced_var_names:
                update_vars.append(LP_symbol_map.bySymbol[LP_name])
        # add the first second-stage variable (if not present)
        if StageToVariableMap[secondstage._name][0][0] not in \
           referenced_var_names:
            update_vars.append(StageToVariableMap[secondstage._name][0][1])
        obj_repn.variables = list(obj_repn.variables) + \
                             update_vars
        obj_repn.linear = list(obj_repn.linear) + \
                          [0.0 for vardata in update_vars]
    else:
        raise RuntimeError("Unexpected objective representation")

    #
    # Create column (variable) ordering maps for LP files
    #
    column_order = ComponentMap()
    # first-stage variables
    for column_index, (LP_name, vardata, scenario_tree_id) \
          in enumerate(StageToVariableMap[firststage._name]):
        column_order[vardata] = column_index
    # second-stage variables
    for column_index, (LP_name, vardata, scenario_tree_id) \
        in enumerate(StageToVariableMap[secondstage._name],
                     len(column_order)):
        column_order[vardata] = column_index

    #
    # Create row (constraint) ordering maps for LP files
    #
    row_order = ComponentMap()
    # first-stage constraints
    for row_index, (LP_names, condata) \
          in enumerate(StageToConstraintMap[firststage._name]):
        row_order[condata] = row_index
    # second-stage constraints
    for row_index, (LP_names, condata) \
        in enumerate(StageToConstraintMap[secondstage._name],
                     len(row_order)):
        row_order[condata] = row_index

    #
    # Write the ordered LP file
    #
    with pyomo.repn.plugins.cpxlp.ProblemWriter_cpxlp() as writer:
        lp_filename = os.path.join(output_prefix, output_prefix+".lp")
        io_options = {'symbolic_solver_labels': True,
                      'column_order': column_order,
                      'row_order': row_order}
        output_filename, LP_symbol_map = writer(reference_model,
                                                lp_filename,
                                                lambda x: True,
                                                io_options)
        assert output_filename == lp_filename

    with open(os.path.join(output_prefix, output_prefix+".tim"),'w') as f:
        f.write('TIME\n')
        f.write('PERIODS\tLP\n')
        f.write('\t'+str(StageToVariableMap[firststage._name][0][0])+
                '\t'+str(LP_symbol_map.byObject[id(reference_scenario._instance_objective)])+
                '\tTIME1\n')
        LP_names = StageToConstraintMap[secondstage._name][0][0]
        if len(LP_names) == 1:
            # equality constraint
            assert (LP_names[0].startswith('c_e_') or \
                    LP_names[0].startswith('c_l_') or \
                    LP_names[0].startswith('c_u_'))
            stage2_row_start = LP_names[0]
        else:
            # range constraint (assumed the LP writer outputs
            # the lower range constraint first)
            LP_names = sorted(LP_names)
            assert (LP_names[0].startswith('r_l_') or \
                    LP_names[0].startswith('r_u_'))
            stage2_row_start = LP_names[0]

        f.write('\t'+str(StageToVariableMap[secondstage._name][0][0])+
                '\t'+stage2_row_start+
                '\tTIME2\n')
        f.write('ENDATA\n')

    constraint_LP_names = ComponentMap()
    for stage_name in StageToConstraintMap:
        for LP_names, constraint_data in StageToConstraintMap[stage_name]:
            constraint_LP_names[constraint_data] = LP_names

    with open(os.path.join(output_prefix, output_prefix+".sto"),'w') as f:
        f.write('STOCH retail\n')
        f.write('INDEP\tDISCRETE\n')
        for param, table in reference_model.PySP_StochasticParameters.items():
            nominal_value = value(param)
            for param_value, probability in table:
                rhs_template = "RHS\t%s\t%+.17g\t%.17g\n"
                if param in reference_model.PySP_StochasticRHS:
                    constraint_data = reference_model.PySP_StochasticRHS[param]
                    constraint_repn = constraint_data.parent_block().canonical_repn[constraint_data]
                    body_constant = constraint_repn.constant
                    if body_constant is None:
                        body_constant = 0.0
                    assert isinstance(constraint_repn, LinearCanonicalRepn)
                    #
                    # **NOTE: In the code that follows we assume the LP writer
                    #         always moves constraint body constants to the rhs
                    #
                    LP_names = constraint_LP_names[constraint_data]
                    if constraint_data.equality:
                        assert len(LP_names) == 1
                        assert LP_names[0].startswith('c_e_')
                        assert constraint_data.lower is param
                        # equality constraint
                        f.write(rhs_template % (LP_names[0], param_value-body_constant, probability))
                    elif ((constraint_data.lower is not None) and \
                          (constraint_data.upper is not None)):
                        # range constraint

                        assert len(LP_names) == 2
                        assert constraint_data.lower is param
                        assert constraint_data.upper is param

                        for LP_name in LP_names:
                            if LP_name.startswith('r_l_'):
                                # lower bound
                                f.write(rhs_template % (LP_names, param_value-body_constant, probability))
                            else:
                                assert LP_name.startswith('r_u_')
                                # upper_bound
                                f.write(rhs_template % (LP_names, param_value-body_constant, probability))
                    elif constraint_data.lower is not None:
                        # lower bound
                        assert len(LP_names) == 1
                        assert LP_names[0].startswith('c_l_')
                        assert constraint_data.lower is param
                        f.write(rhs_template % (LP_names[0], param_value-body_constant, probability))
                    else:
                        # upper bound
                        assert constraint_data.upper is not None
                        assert len(LP_names) == 1
                        assert LP_names[0].startswith('c_u_')
                        assert constraint_data.upper is param
                        f.write(rhs_template % (LP_names[0], param_value-body_constant, probability))
        f.write('ENDATA\n')

    print("Output saved to: "+output_prefix)

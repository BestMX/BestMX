# LIBTBX_SET_DISPATCHER_NAME bestmx

from __future__ import absolute_import, division
from dials.array_family import flex
import iotbx.phil

help_message = '''

'''

phil_scope = iotbx.phil.parse("""\
  force_static = False
    .type = bool
    .help = "For a scan varying model, force static prediction"

  buffer_size = 0
    .type = int(value_min=0)
    .help = "Calculate predictions within a buffer zone of n images either"
            "size of the scan"

  d_min = None
    .type = float(value_min=0)
    .help = "Minimum d-spacing of predicted reflections"

""", process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json | experiments.json | image_*.cbf" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments)==0 or len(reflections)==0:
    parser.print_help()
    exit()

  crystal = experiments.crystals()[0]
  space_group = crystal.get_space_group()
  unit_cell = crystal.get_unit_cell()
  print space_group.info()
  print unit_cell

  expt = experiments[0]

  import copy
  import math
  # Hack to make the predicter predict reflections outside of the range
  # of the scan
  expt_input = copy.deepcopy(expt)
  scan = expt.scan
  image_range = scan.get_image_range()
  oscillation = scan.get_oscillation()
  scan.set_image_range((1, int(math.ceil(360/oscillation[1]))))
  scan.set_oscillation((0, oscillation[1]))
  print scan

  # Populate the reflection table with predictions
  predicted = flex.reflection_table.from_predictions(
    expt,
    force_static=params.force_static,
    dmin=params.d_min)
  predicted['id'] = flex.int(len(predicted), 0)

  print len(predicted)

  from cctbx import crystal, miller
  space_group = space_group.build_derived_reflection_intensity_group(anomalous_flag=True)
  cs = crystal.symmetry(unit_cell=unit_cell, space_group=space_group)

  ms = miller.set(cs, indices=predicted['miller_index'], anomalous_flag=True)
  ma = miller.array(ms, data=flex.double(ms.size(),1),
                    sigmas=flex.double(ms.size(), 1))

  d_star_sq = ma.d_star_sq().data()
  n_bins = 1
  binner = ma.setup_binner_d_star_sq_step(
    d_star_sq_step=(flex.max(d_star_sq)-flex.min(d_star_sq)+1e-8)/n_bins)
  image_number = predicted['xyzcal.px'].parts()[2]
  print flex.min(image_number)
  print flex.max(image_number)
  #dose = flex.size_t(list(flex.floor(image_number).iround()))
  angle_deg = predicted['xyzcal.mm'].parts()[2] * 180/math.pi
  dose = flex.size_t(list(flex.floor(angle_deg).iround()))
  range_width = 1
  range_min = flex.min(dose) - range_width
  range_max = flex.max(dose)
  n_steps = 2 + int((range_max - range_min) - range_width)

  binner_non_anom = ma.as_non_anomalous_array().use_binning(
    binner)
  n_complete = flex.size_t(binner_non_anom.counts_complete()[1:-1])


    
  range_for_50 = []
  range_for_99 = []
  from xia2.Modules.PyChef2 import ChefStatistics
  step = 1
  for i in range(0, 360, step):
    sel = dose < step
    dose.set_selected(sel, dose.select(sel) + 360)
    dose -= flex.min(dose)
    chef_stats = ChefStatistics(
      ma.indices(), ma.data(), ma.sigmas(),
      ma.d_star_sq().data(), dose, n_complete, binner,
      ma.space_group(), ma.anomalous_flag(), n_steps)

  
    ieither_completeness = chef_stats.ieither_completeness()
    iboth_completeness = chef_stats.iboth_completeness()


    range_for_50.append(min((ieither_completeness > 0.5).iselection()))
    range_for_99.append(min((ieither_completeness > 0.99).iselection()))
    
  from matplotlib import pyplot
  pyplot.figure()
  pyplot.plot(range_for_50)
  pyplot.plot(range_for_99)  
  
#  pyplot.scatter(range(iboth_completeness.size()), iboth_completeness)
  pyplot.show()
  

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

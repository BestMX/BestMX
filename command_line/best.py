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

  d_min = None
    .type = float(value_min=0)
    .help = "Minimum d-spacing of predicted reflections"

  space_group = None
    .type = space_group
    .help = "Optionally override the space group."

  dose_rate = None
    .type = float(value_min=0)
    .help = " Dose rate for characterization images" 
""", process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] experiments.json | integrated.pickle" %(
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
  reflections = reflections[0]
  cryst = experiments.crystals()[0]
  unit_cell = cryst.get_unit_cell()
  if params.space_group is not None:
    space_group = params.space_group.group()
    assert space_group.is_compatible_unit_cell(unit_cell), unit_cell
  else:
    space_group = cryst.get_space_group()
  print space_group.info()
  print unit_cell

  expt = experiments[0]
  from cctbx import miller, crystal
  from mmtbx.scaling.data_statistics import wilson_scaling
  sel = reflections.get_flags(reflections.flags.integrated_sum)
  reflections = reflections.select(sel)
  cs = crystal.symmetry(unit_cell=unit_cell, space_group=space_group)
  ms = miller.set(cs, indices=reflections['miller_index'], anomalous_flag=True)  
  intensities = miller.array( ms, data=reflections['intensity.sum.value'],
                              sigmas=flex.sqrt(reflections['intensity.sum.variance']))
  intensities.set_observation_type_xray_intensity()
  d_star_sq = intensities.d_star_sq().data()
  n_bins = 20
#  binner = intensities.setup_binner_d_star_sq_step(
#    d_star_sq_step=(flex.max(d_star_sq)-flex.min(d_star_sq)+1e-8)/n_bins)
  binner = intensities.setup_binner_counting_sorted(n_bins=n_bins)
  # wilson = intensities.wilson_plot(use_binning=True)
  # wilson.show()
  # from matplotlib import pyplot
  # pyplot.figure()
  # pyplot.scatter(wilson.binner.bin_centers(2), wilson.data[1:-1])
  # pyplot.show()
  intensities = intensities.merge_equivalents().array()
  wilson=wilson_scaling(intensities, n_residues=200)
  wilson.iso_scale_and_b.show()

  from matplotlib import pyplot
  pyplot.figure()
  pyplot.scatter(wilson.d_star_sq, wilson.mean_I_obs_data, label='Data')
  pyplot.plot(wilson.d_star_sq, wilson.mean_I_obs_theory, label='theory')
  pyplot.plot(wilson.d_star_sq, wilson.mean_I_normalisation, label='smoothed' )
  pyplot.yscale('log')
  pyplot.legend()
  pyplot.show()
  
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


    
  ranges_dict = {}
  completeness_levels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99]
  for c in completeness_levels:
          ranges_dict[c]=[]
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

    for c in completeness_levels:
      ranges_dict[c].append(min((ieither_completeness > (c/100)).iselection()))



  from matplotlib import pyplot
  pyplot.figure()
  for c in completeness_levels:
    pyplot.plot(ranges_dict[c], label=str(c))
#  pyplot.plot(range_for_50)
#  pyplot.plot(range_for_99)  
  
#  pyplot.scatter(range(iboth_completeness.size()), iboth_completeness)
  pyplot.legend()
  pyplot.show()



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

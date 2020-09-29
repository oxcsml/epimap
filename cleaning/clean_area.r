source('cleaning/epiclean.r')

option_list = list(
  make_option(c("-t", "--task_id"), type="integer", default=0, help="Task ID for Slurm usage. Maps to area_index.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

area_index = opt$task_id

epiclean(area_index)


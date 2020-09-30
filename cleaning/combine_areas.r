source('cleaning/epiclean.r')

epiclean_options = epiclean_options()

option_list = list(
  make_option(c("-c", "--clean_directory"), type="character", default=epiclean_options$clean_directory, help="Directory to put cleaned results in.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

for (o in names(opt)) {
    epiclean_options[o] = opt[o]
}

epiclean_combine(epiclean_options)

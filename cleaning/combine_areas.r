source('cleaning/epiclean.r')

opt = epiclean_get_cmdline_options()
epiclean_combine(opt)


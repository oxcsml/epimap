source('cleaning/epiclean.r')

opt = epiclean_get_cmdline_options()
area_index = opt$task_id
print(area_index)
epiclean(area_index, opt)


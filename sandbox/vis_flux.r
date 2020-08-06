areanames = uk_cases[1:N,2]

for (i in 1:N) {
  js = order(-flux[i,])[1:80]
  print(areanames[i])
  print(data.frame(area=areanames[js],flux=cumsum(flux[i,js]),flux_sd=flux_sd[i,js],pop=population[js],dist=geodist[i,js]))
  invisible(readline())
}

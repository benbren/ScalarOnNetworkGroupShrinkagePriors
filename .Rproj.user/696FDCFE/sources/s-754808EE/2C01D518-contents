library(profvis)

profvis({
  get_posteriors(5,1000, 1000,1000, debug = F)
})


Rprof(tmp <- tempfile())
get_posteriors(5,1000, 1000,1000, debug = F)
Rprof()
summaryRprof(tmp)
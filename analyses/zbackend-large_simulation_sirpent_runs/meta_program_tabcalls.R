library(tidyverse)

#............................................................
# expand out metaprogram lvls
#...........................................................
fctrowslvls <- tidyr::expand_grid(
  modname = c("degreedist", "modularity", "unity", "cluster", "NEdynamicity"),
  beta = c(0.005, 0.01, 0.05, 0.075, 0.1),
  durI = seq(3, 15, by = 3))



#............................................................
# metaprogram function
#...........................................................
meta_make_tbas <- function(modname, beta, durI) {
  header <- paste("####", "B:", beta, ";", "D:", durI, sep = " ")
  fxncall <- paste("plotSimRet(df = ogsimrets, modname = ", paste0("\'", modname, "\'"), ", beta = ",
  beta, ", durI = ", paste0(durI, ")"))
  out <- capture.output( cat( header,"\n","```{r}","\n",fxncall,"\n","```", "\n","\n") )
  return(out)
}


#............................................................
# run out
#...........................................................
out <- fctrowslvls %>%
  dplyr::mutate(mycall = purrr::pmap(., meta_make_tbas))
readr::write_lines(unlist(out$mycall), file = "analyses/zbackend-large_simulation_sirpent_runs/meta_program_tabcalls_out.txt")

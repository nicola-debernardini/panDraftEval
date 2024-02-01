# get metabolic by-products
get.produced.metabolites <- function(mod) {
  
  # get MTF solution
  sol.mtf <- optimizeProb(mod, algorithm = "mtf")
  dt.mtf  <- data.table(ex = mod@react_id,
                        mtf.flux = sol.mtf@fluxdist@fluxes[1:mod@react_num])
  dt.mtf.tmp <- copy(dt.mtf[grepl("^EX_", ex)])
  
  # this following two lines are there to prevent the case that a nutrient (e.g. L-Lactate)
  # from the environment is taken up, and thus enables the production of D-Lactate.
  dt.mtf.tmp[mtf.flux > 0, mtf.flux := 0]
  model.tmp <- changeBounds(mod, react = dt.mtf.tmp$ex, lb = dt.mtf.tmp$mtf.flux)
  
  
  # get FV solution
  sol.fv <- fluxVar(model.tmp, react = mod@react_id[grep("^EX_", mod@react_id)])
  
  dt <- data.table(ex = rep(mod@react_id[grep("^EX_", mod@react_id)],2),
                   rxn.name = rep(mod@react_name[grep("^EX_", mod@react_id)],2),
                   dir = c(rep("l",length(grep("^EX_", mod@react_id))),rep("u",length(grep("^EX_", mod@react_id)))),
                   fv = sol.fv@lp_obj)
  dt <- dcast(dt, ex + rxn.name ~ dir, value.var = "fv")[(u>1e-4 & l >= 0) | (u > 1)]
  
  dt <- merge(dt, dt.mtf, by = "ex")
  
  return(dt[order(-u)])
}
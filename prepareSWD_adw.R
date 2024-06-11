
# función del paquete 'SDMtune' con un bug corregido

prepareSWD_adw <-
function (species, env, p = NULL, a = NULL, categorical = NULL) 
{
  df_coords <- data.frame(X = numeric(), Y = numeric())
  df_data <- data.frame()
  dfs <- list(p, a)
  text <- c("presence", "absence/background")
  values <- c(1, 0)
  pa <- c()
  for (i in 1:2) {
    if (!is.null(dfs[[i]])) {
      coords <- as.data.frame(dfs[[i]])
      colnames(coords) <- c("X", "Y")
      message("Extracting predictor information for ", 
              text[i], " locations...")
      data <- as.data.frame(raster::extract(env, coords))
      aa <- names(data)
      index <- stats::complete.cases(data)
      discarded <- nrow(data) - sum(index)
      if (discarded > 0) {
        data <- as.data.frame(data[index, ])
        names(data) <- aa
        coords <- coords[index, ]
        message("Info: ", discarded, " ", text[i], ifelse(discarded == 
                                                            1, " location is", " locations are"), " NA for some environmental variables, ", 
                ifelse(discarded == 1, "it is ", "they are "), 
                "discarded!")
      }
      df_coords <- rbind(df_coords, coords)
      df_data <- rbind(df_data, data)
      pa <- c(pa, rep(values[i], nrow(coords)))
    }
  }
  
  if (!is.null(categorical)) {
    for (i in categorical) {
      df_data[, i] <- as.factor(df_data[, i])
    }
  }
  rownames(df_coords) <- NULL
  rownames(df_data) <- NULL
  swd <- SWD(species = species, coords = df_coords, data = df_data, 
             pa = pa)
  return(swd)
}

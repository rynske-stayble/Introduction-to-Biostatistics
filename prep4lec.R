required.packages <- c("KernSmooth", "scatterplot3d", "maps", "mapdata", "plotly", "ape", "gplots", "fields", "cluster")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages, repos="http://cran.us.r-project.org")
}

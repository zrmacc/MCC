# Prevent R from writing Rplots.pdf when tests create plots (e.g. show(ggplot)).
# Any new default device opened during tests goes to null instead.
options(device = function(...) grDevices::pdf(file = NULL))

### Code externalisation
### https://yihui.name/knitr/demo/externalization/

# if(Sys.getenv("RSTUDIO") == "1"){
  # setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # This will be NULL if no document is open in Rstudio
  # logger::log_info("Working dir switched to codebase.")
# }else{
  # setwd(utils::getSrcDirectory(function(){})[1])
# }

## Adjust timeout for utils::download.file()
options(timeout = 60*30) # 30 minutes

## @knitr global_knitr_options
knitr::opts_chunk$set(
  fig.width = 12, 
  fig.height = 8, 
  fig.path = "markdown_figs/", 
  dev = "png", 
  eval = TRUE, 
  echo = FALSE, 
  warning = FALSE, 
  message = FALSE, 
  tidy = FALSE
)
## This switch allows for document-type dependent output (e.g. interactive graphs) https://trinkerrstuff.wordpress.com/2014/11/18/rmarkdown-alter-action-depending-on-document/
document_output_type <- knitr::opts_knit$get("rmarkdown.pandoc.to")
# print(document_output_type)
## ---- end


## @knitr package_management
# source("01a_install_packages.R")
## ---- end

## @knitr colorscheme
### Colours
library(RColorBrewer) ## Colour palette
HEATMAP_COL <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
tmp <- brewer.pal(9, "Set1")
mycol_2 <- rep(NA, length(tmp))
for (i in 1:length(tmp)) {
  mycol_2[i] <- rgb(red=col2rgb(tmp[i])[1], green=col2rgb(tmp[i])[2], blue=col2rgb(tmp[i])[3], alpha=0.5*255, maxColorValue = 255)
}
POINT_COL <- mycol_2
PAIRED <- brewer.pal(12, "Paired")[1:10]
tmp <- expand.grid( list(a = c("in", "out"), b = c("blue", "green", "red", "orange", "purple") ))
names(PAIRED) <- paste(tmp[, 2], tmp[, 1], sep="_")
## BTX palette
BTX_COL <-
	list(
		blue_d_in = c(189, 197, 255),
		blue_d_out = c(0, 62, 193),
		blue_l_in = c(179, 223, 255),
		blue_l_out = c(0, 151, 229),
		green_in = c(195, 249, 195),
		green_out = c(0, 169, 0),
		pink_in = c(254, 195, 252),
		pink_out = c(184, 56, 180),
		purple_in = c(212, 178, 255),
		purple_out = c(151, 73, 254),
		gray_in = c(214, 214, 214),
		gray_out = c(0, 0, 0)	
	)
tmp <- rep(NA, length(BTX_COL))	
for (i in 1:length(BTX_COL)) {
  tmp[i] <- 
	rgb(
		red = (BTX_COL[[i]][1]), 
		green = (BTX_COL[[i]][2]), 
		blue = (BTX_COL[[i]][3]),
		maxColorValue = 255
	)
}
names(tmp) <- names(BTX_COL)
BTX_COL <- as.list(tmp)
## ---- end

## @knitr logo

add_btx_logo <- function(p, l = "./IMG/2k B.png") {
	p <-
		stylehaven::add_logo(
			p, 
			image = l, 
			position = "right", 
			height = 0.025
		)
	p
	# Add horizontal line on top
	# It goes from x = 0 (left) to x = 1 (right) on the very top of the chart (y = 1)
	# You can think of 'gp' and 'gpar' as 'graphical parameters'.
	# There we indicate the line color and width
	# grid.lines(
	  # x = c(0, 1),
	  # y = 1,
	  # gp = gpar(col = BTX_COL$blue_d_out, lwd = 4)
	# )

	# # Add rectangle on top-left
	# # lwd = 0 means the rectangle does not have an outer line
	# # 'just' gives the horizontal and vertical justification
	# grid.rect(
	  # x = 0,
	  # y = 1,
	  # width = 0.05,
	  # height = 0.01,
	  # just = c("left", "top"),
	  # gp = gpar(fill = BTX_COL$blue_d_out, lwd = 0)
	# )
	
}
## ---- end


## @knitr fonts
## Fonts
## https://fonts.google.com/
library(sysfonts)
# sysfonts::font_families_google(db_cache = TRUE, handle = curl::new_handle()) # all options
sysfonts::font_add_google("Lato", "lato")
sysfonts::font_add_google("Open Sans", "opensans")
# sysfonts::font_add_google("Raleway", "raleway")
# sysfonts::font_add_google("Fira Sans", "firasans")
showtext::showtext_auto()

## Plotting themes
library(ggplot2)
old_theme <- 
  ggplot2::theme_set(ggplot2::theme_minimal()) +
  ggplot2::theme_update(
    legend.title = element_blank(),
    # legend.justification = c(0, 1), 
    # legend.position = c(.1, 1.075),
    legend.background = element_blank(),
    axis.title= element_text(family = "opensans", size = 10),
    plot.title = element_text(family = "lato", size = 20, margin = margin(b = 10)),
    plot.subtitle = element_text(family = "opensans", size = 10, color = "darkslategrey", margin = margin(b = 25)),
    plot.caption = element_text(family = "opensans", size = 8, margin = margin(t = 10), color = "grey70", hjust = 0)
  )
## ---- end





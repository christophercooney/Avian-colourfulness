
### Mapping colourfulness ###

rm(list=ls())

setwd("~/Dropbox/Projects/Current/Bird_colouration/Passerine_colourfulness/")

library(pavo)
library(sf)
library(parallel)
library(RColorBrewer)
library(letsR)
library(viridis)
library(dispRity)
library(spatstat)
library(maptools)
library(rgeos)
library(cetcolor)

source("Code/00_Custom_mapping_functions_v5.R")

# ------ #

## Custom mapping functions

sp_trait <- function(presab_dat, trait_dat, traits, type=c("average","q.richness","q.proportion"), logtrait=FALSE, geomean=FALSE, range_weights=FALSE, quantile=c(1,2,3,4), base_raster=NULL, plot_map=TRUE, clip_by_richness=0, cores=1) {
  
  # Libraries
  require(sf)
  require(fasterize)
  require(dplyr)
  require(raster)
  require(dispRity)
  require(ggplot2)
  
  # Create an empty raster
  base_raster2 <- base_raster
  base_raster2@data@inmemory <- TRUE
  
  myFun <- function(i) {
    
    tmp_val <- NA
    if (rowSums(presab_dat[i,,drop=F], na.rm=T) > 0) {
      sp_index <- !is.na(presab_dat[i,])
      tmp_nms <- colnames(presab_dat)[sp_index]
      idx <- which(trait_dat$binomial %in% tmp_nms) 
      if (length(idx)>0) {
        trait_vec <- as.matrix(trait_dat[idx, traits], ncol=length(traits))
        if (type == "average") {
          if (logtrait) {
            trait_vec <- log(trait_vec)
          }
          if (geomean) {
            if (range_weights) {
              range_vec <- as.matrix(trait_dat[idx, "ncells"], ncol=1)
              tmp_val <- weighted.geomean(trait_vec, 1/range_vec, na.rm=T)
            } else {
              tmp_val <- geomean(trait_vec, na.rm=TRUE)
            }
          } else {
            if (range_weights) {
              range_vec <- as.matrix(trait_dat[idx, "ncells"], ncol=1)
              tmp_val <- weighted.mean(trait_vec, 1/range_vec, na.rm=T)
            } else {
              tmp_val <- mean(trait_vec, na.rm=TRUE)
            }
          }
        }
        if (type == "q.richness") {
          tmp_val <- sum(trait_vec == quantile)
        }
        if (type == "q.proportion") {
          tmp_val <- sum(trait_vec == quantile) / nrow(trait_vec)
        }
      } else {
        tmp_val <- NA
      }
    }
    return(tmp_val)
  }
  
  output <- mclapply(1:dim(presab_dat)[1], myFun, mc.cores=cores) 
  
  tr_val <- unlist(output)
  
  base_raster2@data@values <- tr_val
  base_raster2@data@min <- min(tr_val, na.rm=TRUE)
  base_raster2@data@max <- max(tr_val, na.rm=TRUE)
  
  if (clip_by_richness > 0) {
    rich <- rowSums(presab_dat)
    values(base_raster2)[which(rich<clip_by_richness)] <- NA
  }
  
  if (plot_map==TRUE) {
    plot(base_raster2)
    #map(interior=FALSE, add=TRUE, lwd=0.5)
    plot(world_lines, add=T, col="black")
  }
  
  return(base_raster2)
  
}

# ----- #

color.legend2 <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
                           gradient = "x", border, ...) 
{
  oldcex <- par("cex")
  par(xpd = TRUE, cex = cex)
  gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
                gradient = gradient, border=border)
  if (gradient == "x") {
    xsqueeze <- (xr - xl)/(2 * length(rect.col))
    textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      texty <- yb - 0.2 * strheight("O")
      textadj <- c(0.5, 1)
    }
    else {
      texty <- yt + 0.2 * strheight("O")
      textadj <- c(0.5, 0)
    }
  }
  else {
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      textx <- xr + 0.2 * strwidth("O")
      textadj <- c(0, 0.5)
    }
    else {
      textx <- xl - 0.2 * strwidth("O")
      textadj <- c(1, 0.5)
    }
  }
  text(textx, texty, labels = legend, adj = textadj, ...)
  par(xpd = FALSE, cex = oldcex)
}

# ----- #

plot_bird_raster <- function(raster_data=NULL, scale="raw", nbins=20, breaks = NULL, col=viridis(20), cex=1.5, legend_pos=c(-8000000,-8000000, 8000000,-7500000), plot_world=TRUE, newproj=NULL){
  
  require(plotrix)
  require(ggplot2)
  require(viridis)
  
  if(scale=="eq_number") { 
    cuts <- cut_number(raster_data@data@values, n = nbins)
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="eq_width") {
    cuts <- cut_interval(raster_data@data@values, n = nbins)
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="raw") {
    cuts <- cut_interval(raster_data@data@values, n = length(na.omit(unique(raster_data@data@values))))
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="predefined") {
    cuts <- cut(raster_data@data@values, breaks = breaks, include.lowest = T)
    raster_data@data@values <- as.numeric(cuts)
  }
  
  if (scale=="predefined") {
    brks <- round(breaks[seq(1,length(breaks),length.out = 6)],2)
  } else {
    cuts_vec <- sort(na.omit(unique(cuts)))
    idx <- c(round(quantile(1:length(cuts_vec), c(0,0.2,0.4,0.6,0.8,1))))
    brks <- as.character(cuts_vec[idx])
    brks <- gsub("\\[|\\]|\\)|\\(", "", brks)
    brks <- as.numeric(unlist(strsplit(brks, ",")))[c(1,3,5,7,9,12)]
  }
  
  if (!is.null(newproj)) {
    fworld <- st_transform(world, crs = newproj)
    fraster_data <- projectRaster(raster_data, crs = newproj)
  } else {
    fworld <- world
    fraster_data <- raster_data
  }
  
  if (plot_world==TRUE) {
    plot(fworld$geometry, border="grey80", col="grey80")
    par(bty="n")
    plot(fraster_data, add=T, legend=F, col=col, axes=FALSE)
  }
  if (plot_world==FALSE) {par(bty="n")
    plot(fworld$geometry, border=NA, col=NA)
    par(bty="n")
    plot(fraster_data, add=T, legend=F, col=col, axes=FALSE)
  }
  
  color.legend2(legend_pos[1], legend_pos[2], legend_pos[3], legend_pos[4], brks, rect.col=col, bty="n", align="rb", cex=cex, border=NA)
  
}

# ------ #

## Get data
res <- read.csv("Outputs/02_Passerine_MF_colourfulness_scores_means.csv", strings=F)
pam <- readRDS("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs/PAM_birds_Behrman50km_Pres12_Orig12_Seas12_liberal_clipped.rds")

## Spatial parameters
null_rast <- pam[[2]]
pam <- pam[[1]]
colnames(pam) <- gsub(" ", "_", colnames(pam))
units = 100 # in km
newproj = as.character(crs(null_rast)) # Behrmann cylindrical equal area projection (standard parallels at 30 deg)

## Read in and project world (maptools: wrld_simpl)
data(wrld_simpl)
world = st_as_sf(wrld_simpl)
world = st_transform(world, crs = newproj)
world_lines = st_cast(world, to = "MULTILINESTRING")

# ------ #

### Identify top 25% colourful species based on sex-specific percentile assignments
res$volume.uv.quantile <- as.numeric(res$volume.uv.percentile > 75)
res$volume.vs.quantile <- as.numeric(res$volume.vs.percentile > 75)
res$loci.uv.quantile <- as.numeric(res$loci.uv.percentile > 75)
res$loci.vs.quantile <- as.numeric(res$loci.vs.percentile > 75)

# ------------------- #

### Global inputs
presab_dat <- pam
presab_dat <- presab_dat[,colnames(presab_dat) %in% res$binomial]
cores <- 12
trait_dat <- res
trait_dat <- trait_dat[trait_dat$binomial %in% colnames(presab_dat),]
base_raster <- null_rast
clip_by_richness <- 0
plot_map <- FALSE
logtrait <- FALSE
geomean <- FALSE

# ------ #

## Calculate cell averages using different approaches and datasets
res.cells <- data.frame(id=c(1:ncell(null_rast)), xyFromCell(null_rast, cell=c(1:ncell(null_rast))), nspp=rowSums(presab_dat, na.rm=T))
colnames(res.cells) <- c("id","long","lat","nspp")

# Mean - not range weighted
range_weights <- FALSE
res.cells$volume.uv.m.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.uv.f.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.m.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.f.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.m.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.f.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.m.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.f.mean <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values

# Mean - range weighted
range_weights <- TRUE
res.cells$volume.uv.m.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.uv.f.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.m.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.f.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.m.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.f.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.m.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.f.mean.wtd <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values

# Top 25% proportion
range_weights <- FALSE
res.cells$volume.uv.m.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.uv.f.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.m.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.uv.f.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.m.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "volume.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$volume.vs.f.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "volume.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.m.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="M",], "loci.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
res.cells$loci.vs.f.top25.prop <- sp_trait(presab_dat, trait_dat[trait_dat$sex=="F",], "loci.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values

# save
saveRDS(res.cells, "Outputs/04_Passerine_MF_colourfulness_scores_means_cell_assemblages.rds")

# ------------------------------------------- #

## Plot mean per-species colourfulness for grid cells

library(pavo)
library(sf)
library(parallel)
library(RColorBrewer)
library(letsR)
library(viridis)
library(dispRity)
library(spatstat)
library(maptools)
library(rgeos)
library(cetcolor)
library(inlmisc) # https://waterdata.usgs.gov/blog/tolcolors/
library(ggplot2)

# ------ #

producePlotRaster <- function (variable, clip_by_richness) {
  values <- res.cells[,c(variable)]
  values[res.cells$nspp < clip_by_richness] <- NA
  plotRaster <- null_rast
  values(plotRaster) <- values
  return(plotRaster)
}

# ------ #

## Plot rasters

pal <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00') # # colorRampPalette(GetColors(nbins, scheme = "ocean"))

## uv - mean
png(paste0("Plots/rasters/mean.uv.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfcol=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.uv.m.mean", min_cell_richness)
prastf <- producePlotRaster("volume.uv.f.mean", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.uv.m.mean", min_cell_richness)
prastf <- producePlotRaster("loci.uv.f.mean", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()

## uv - mean range weighted
png(paste0("Plots/rasters/wtd.mean.uv.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfrow=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.uv.m.mean.wtd", min_cell_richness)
prastf <- producePlotRaster("volume.uv.f.mean.wtd", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.uv.m.mean.wtd", min_cell_richness)
prastf <- producePlotRaster("loci.uv.f.mean.wtd", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()

## uv - top 25% proportion
png(paste0("Plots/rasters/top25.prop.uv.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfcol=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.uv.m.top25.prop", min_cell_richness)
prastf <- producePlotRaster("volume.uv.f.top25.prop", min_cell_richness)
nbins <- 8; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((breaks[i]),2), nsmall=2), "-", format(round((breaks[i+1]),2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((breaks[i]),2), nsmall=2), "-", format(round((breaks[i+1]),2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.uv.m.top25.prop", min_cell_richness)
prastf <- producePlotRaster("loci.uv.f.top25.prop", min_cell_richness)
nbins <- 8; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(breaks[i],2), nsmall=2), "-", format(round(breaks[i+1],2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(breaks[i],2), nsmall=2), "-", format(round(breaks[i+1],2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()

## vs - mean
png(paste0("Plots/rasters/mean.vs.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfrow=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.vs.m.mean", min_cell_richness)
prastf <- producePlotRaster("volume.vs.f.mean", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.vs.m.mean", min_cell_richness)
prastf <- producePlotRaster("loci.vs.f.mean", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()


## vs - mean range weighted
png(paste0("Plots/rasters/wtd.mean.vs.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfrow=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.vs.m.mean.wtd", min_cell_richness)
prastf <- producePlotRaster("volume.vs.f.mean.wtd", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((10^breaks[i])*1000,2), nsmall=2), "-", format(round((10^breaks[i+1])*1000,2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.vs.m.mean.wtd", min_cell_richness)
prastf <- producePlotRaster("loci.vs.f.mean.wtd", min_cell_richness)
nbins <- 10; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(10^breaks[i],0), nsmall=0), "-", format(round(10^breaks[i+1],0), nsmall=0))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()

## vs - top 25% proportion
png(paste0("Plots/rasters/top25.prop.vs.berhmann.png"), height=7.5, width=16, res=600, units="in")
par(mfcol=c(2,2), mar=c(1,1,1,1))
min_cell_richness <- 5
newproj <- NULL
lpos <- c(-6000,-18150, 6000,-17750) # i.e. not plotted; correct -> c(-6000,-8150, 6000,-7750)
lbpos <- c(-16000, 2500)
palfun <- colorRampPalette(pal)
# volume
prastm <- producePlotRaster("volume.vs.m.top25.prop", min_cell_richness)
prastf <- producePlotRaster("volume.vs.f.top25.prop", min_cell_richness)
nbins <- 8; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((breaks[i]),2), nsmall=2), "-", format(round((breaks[i+1]),2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round((breaks[i]),2), nsmall=2), "-", format(round((breaks[i+1]),2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
# loci
prastm <- producePlotRaster("loci.vs.m.top25.prop", min_cell_richness)
prastf <- producePlotRaster("loci.vs.f.top25.prop", min_cell_richness)
nbins <- 8; pal <- palfun(nbins)
breaks <- quantile(c(values(prastm),values(prastf)), seq(0,1,length.out=nbins+1), na.rm=T) # i.e. "eq_number"
plot_bird_raster(raster_data=prastm, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(breaks[i],2), nsmall=2), "-", format(round(breaks[i+1],2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
plot_bird_raster(raster_data=prastf, scale="predefined", breaks = breaks, col=pal, cex=0.75, newproj=newproj, legend_pos = lpos)
lb <- c(); for (i in 1:(length(breaks)-1)) { lb <- c(lb, paste0(format(round(breaks[i],2), nsmall=2), "-", format(round(breaks[i+1],2), nsmall=2))) }; legend(lbpos[1], lbpos[2], legend = lb, pch = 22, pt.bg = pal, cex = 1, pt.cex = 1.2)
dev.off()


# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #


### SAR analysis

rm(list=ls())

setwd("~/Dropbox/Projects/Current/Bird_colouration/Passerine_colourfulness/")

library(pavo)
library(sf)
library(parallel)
library(RColorBrewer)
library(letsR)
library(viridis)
library(dispRity)
library(spatstat)
library(maptools)
library(rgeos)
#library(rtoolz)
library(spdep)
library(spatialreg)
library(inlmisc) # https://waterdata.usgs.gov/blog/tolcolors/

# ------ #

## Get data
res <- read.csv("Outputs/02_Passerine_MF_colourfulness_scores_means.csv", strings=F)
res.cells <- readRDS("Outputs/04_Passerine_MF_colourfulness_scores_means_cell_assemblages.rds")
pam <- readRDS("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs/PAM_birds_Behrman50km_Pres12_Orig12_Seas12_liberal_clipped.rds")
pam <- data.frame(pam)

## Spatial parameters
null_rast <- pam[[2]]
pam <- pam[[1]]
colnames(pam) <- gsub(" ", "_", colnames(pam))
units = 100 # in km
newproj = as.character(crs(null_rast)) # Behrmann cylindrical equal area projection (standard parallels at 30 deg)

## Read in and project world (maptools: wrld_simpl)
data(wrld_simpl)
world = st_as_sf(wrld_simpl)
world = st_transform(world, crs = newproj)
world_lines = st_cast(world, to = "MULTILINESTRING")

# ------ #

## ECOREGIONS

# Summarise colourfulness by ecoregion

ecor.cells <- readRDS("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Outputs/50km/Behrmann_cea/WWF_ecoregions.rds")
ecor.poly <- read_sf("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Inputs/official_teow/wwf_terr_ecos.shp")
ecor.lookup <- read.csv("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Outputs/50km/Behrmann_cea/WWF_ecoregions_looukp.csv", strings=F)

ecor.poly <- ecor.poly[,c("ECO_ID")]
ecor.poly <- aggregate(ecor.poly, by=list(ecoregion=ecor.poly$ECO_ID), FUN=mean, do_union=T); ecor.poly$ECO_ID <- NULL
ecor.poly = st_transform(ecor.poly, crs = newproj)

pam$ecoregion <- ecor.cells$ecoregion_combined
res.cells$ecoregion <- ecor.cells$ecoregion_combined

fvars <- colnames(res.cells)[-(1:4)]

# -------- #

# Calculate ecoregion means

ecor.means <- aggregate(cbind(volume.uv.m.mean, volume.uv.f.mean, 
                              loci.uv.m.mean, loci.uv.f.mean,
                              volume.vs.m.mean, volume.vs.f.mean, 
                              loci.vs.m.mean, loci.vs.f.mean,
                              volume.uv.m.mean.wtd, volume.uv.f.mean.wtd, 
                              loci.uv.m.mean.wtd, loci.uv.f.mean.wtd,
                              volume.vs.m.mean.wtd, volume.vs.f.mean.wtd, 
                              loci.vs.m.mean.wtd, loci.vs.f.mean.wtd,
                              volume.uv.m.top25.prop, volume.uv.f.top25.prop,
                              loci.uv.m.top25.prop, loci.uv.f.top25.prop,  
                              volume.vs.m.top25.prop, volume.vs.f.top25.prop,
                              loci.vs.m.top25.prop, loci.vs.f.top25.prop) ~ ecoregion, res.cells, mean)

ecor.poly <- ecor.poly[ecor.poly$ecoregion %in% ecor.means$ecoregion,]

ecor.poly$lat <- st_coordinates(st_centroid(ecor.poly))[,"Y"] / units
ecor.poly$abs.lat <- abs(ecor.poly$lat)

ecor.poly$volume.uv.m.mean <- ecor.means$volume.uv.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.uv.f.mean <- ecor.means$volume.uv.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.m.mean <- ecor.means$loci.uv.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.f.mean <- ecor.means$loci.uv.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

ecor.poly$volume.vs.m.mean <- ecor.means$volume.vs.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.vs.f.mean <- ecor.means$volume.vs.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.m.mean <- ecor.means$loci.vs.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.f.mean <- ecor.means$loci.vs.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

ecor.poly$volume.uv.m.mean.wtd <- ecor.means$volume.uv.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.uv.f.mean.wtd <- ecor.means$volume.uv.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.m.mean.wtd <- ecor.means$loci.uv.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.f.mean.wtd <- ecor.means$loci.uv.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

ecor.poly$volume.vs.m.mean.wtd <- ecor.means$volume.vs.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.vs.f.mean.wtd <- ecor.means$volume.vs.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.m.mean.wtd <- ecor.means$loci.vs.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.f.mean.wtd <- ecor.means$loci.vs.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

ecor.poly$volume.uv.m.top25.prop <- ecor.means$volume.uv.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.uv.f.top25.prop <- ecor.means$volume.uv.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.m.top25.prop <- ecor.means$loci.uv.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.uv.f.top25.prop <- ecor.means$loci.uv.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

ecor.poly$volume.vs.m.top25.prop <- ecor.means$volume.vs.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$volume.vs.f.top25.prop <- ecor.means$volume.vs.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.m.top25.prop <- ecor.means$loci.vs.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
ecor.poly$loci.vs.f.top25.prop <- ecor.means$loci.vs.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]

psize <- plotrix:::rescale(log(as.numeric(st_area(ecor.poly))), newrange = c(0.2,1))
pcol <- "#000000"

saveRDS(ecor.poly, "Outputs/05_Ecoregion_Values.rds")

# ----- #

png("Plots/ecoregions/ecoregion_by_latitude_scatterplots_uv.png", height=6, width=7, res=600, units="in")
par(mfrow=c(2,2), mar=c(4,4,2,2))
yrange <- range(c(ecor.poly$volume.uv.m.mean, ecor.poly$volume.uv.f.mean))
plot(ecor.poly$volume.uv.m.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6)*1000,2), nsmall=2), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour volume (ecoregion mean)", line=2, cex=0.75)
plot(ecor.poly$volume.uv.f.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6)*1000,2), nsmall=2), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour volume (ecoregion mean)", line=2, cex=0.75)
yrange <- range(c(ecor.poly$loci.uv.m.mean, ecor.poly$loci.uv.f.mean))
plot(ecor.poly$loci.uv.m.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6),0), nsmall=0), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour loci (ecoregion mean)", line=2, cex=0.75)
plot(ecor.poly$loci.uv.f.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6),0), nsmall=0), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour loci (ecoregion mean)", line=2, cex=0.75)
dev.off()

png("Plots/ecoregions/ecoregion_by_latitude_scatterplots_vs.png", height=6, width=7, res=600, units="in")
par(mfrow=c(2,2), mar=c(4,4,2,2))
yrange <- range(c(ecor.poly$volume.vs.m.mean, ecor.poly$volume.vs.f.mean))
plot(ecor.poly$volume.vs.m.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6)*1000,2), nsmall=2), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour volume (ecoregion mean)", line=2, cex=0.75)
plot(ecor.poly$volume.vs.f.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6)*1000,2), nsmall=2), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour volume (ecoregion mean)", line=2, cex=0.75)
yrange <- range(c(ecor.poly$loci.vs.m.mean, ecor.poly$loci.vs.f.mean))
plot(ecor.poly$loci.vs.m.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6),0), nsmall=0), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour loci (ecoregion mean)", line=2, cex=0.75)
plot(ecor.poly$loci.vs.f.mean ~ ecor.poly$lat, log="", cex=psize, col=pcol, pch=21, bg=paste0(pcol,"25"), ylim=yrange, xlim=c(-75,75), xlab="", ylab="", axes=F)
axis(1, at = seq(-60,60,20), cex.axis=0.75)
axis(2, at = seq(yrange[1],yrange[2],length.out=6), labels = format(round(10^seq(yrange[1],yrange[2],length.out=6),0), nsmall=0), line=-1, las=1, cex.axis=0.75)
mtext(side=1, text="Latitude (ecoregion midpoint)", line=2.5, cex=0.75)
mtext(side=2, text="Colour loci (ecoregion mean)", line=2, cex=0.75)
dev.off()

# ---------- #

# Run SARs

nb <- poly2nb(ecor.poly, row.names = ecor.poly$ecoregion)

fit.vars <- c("volume.uv.m.mean", "volume.uv.f.mean", 
              "loci.uv.m.mean", "loci.uv.f.mean",
              "volume.vs.m.mean", "volume.vs.f.mean", 
              "loci.vs.m.mean", "loci.vs.f.mean",
              "volume.uv.m.mean.wtd", "volume.uv.f.mean.wtd", 
              "loci.uv.m.mean.wtd", "loci.uv.f.mean.wtd",
              "volume.vs.m.mean.wtd", "volume.vs.f.mean.wtd", 
              "loci.vs.m.mean.wtd", "loci.vs.f.mean.wtd",
              "volume.uv.m.top25.prop", "volume.uv.f.top25.prop", 
              "loci.uv.m.top25.prop", "loci.uv.f.top25.prop",
              "volume.vs.m.top25.prop", "volume.vs.f.top25.prop", 
              "loci.vs.m.top25.prop", "loci.vs.f.top25.prop")

out <- as.list(rep(NA, length(fit.vars)))
names(out) <- fit.vars

for (i in 1:length(fit.vars)) {
  print (paste(i, "of", length(fit.vars)))
  var <- fit.vars[i]
  # fit ols
  fit.ols <- lm(scale(ecor.poly[,var,drop=T]) ~ scale(ecor.poly[,"abs.lat",drop=T]))
  # build SAR models with different combinations of distance and neighbourhood style
  neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
  AICvec <- numeric(length = length(neiStyle))
  for (k in 1:length(neiStyle)) {
    nlw <- nb2listw(nb, style=neiStyle[k], zero.policy=TRUE)
    sar_e <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
    AICvec[k] <- AIC(sar_e)			
  }
  bestStyle <- neiStyle[which.min(AICvec)]
  nlw <- nb2listw(nb, style=bestStyle, zero.policy=TRUE)
  fit.sar <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  aic.ols <- AIC(fit.ols)
  aic.sar <- AIC(fit.sar)
  res <- list(fit.ols = fit.ols, fit.sar = fit.sar, aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
  out[[i]] <- res
}

# Summarise
x <- numeric(length(fit.vars))
dff <- data.frame(metric = fit.vars, AIC.ols=x, AIC.sar=x, dAIC.ols=x, dAIC.sar=x,
                  ols.slope=x, ols.pvalue=x,
                  sar.slope=x, sar.pvalue=x, slope.ratio=x,
                  moran.ols = x, moran.ols.pvalue = x, moran.sar = x, moran.sar.pvalue = x, stringsAsFactors=F)

for (i in 1:length(out)) {
  
  fres <- out[[i]]
  
  dff$AIC.ols[i] <- fres$aic.ols
  dff$AIC.sar[i] <- fres$aic.sar
  dff$dAIC.ols[i] <- fres$aic.ols - (min(fres$aic.ols, fres$aic.sar))
  dff$dAIC.sar[i] <- fres$aic.sar - (min(fres$aic.ols, fres$aic.sar))
  dff$ols.slope[i] <- coef(fres$fit.ols)[2]
  dff$ols.pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
  dff$sar.slope[i] <- coef(fres$fit.sar)[2]
  dff$sar.pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
  dff$slope.ratio[i] <- dff$sar.slope[i] / dff$ols.slope[i]
  dff$moran.ols[i] <- fres$moran.ols$estimate[1]
  dff$moran.ols.pvalue[i] <- fres$moran.ols$p.value
  dff$moran.sar[i] <- fres$moran.sar$estimate[1]
  dff$moran.sar.pvalue[i] <- fres$moran.sar$p.value
  
}


# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #


### Do SAR SES analysis

## RANDOMISE SPP TRAIT VALUES WRT LATITUDE

rm(list=ls())

setwd("~/Dropbox/Projects/Current/Bird_colouration/Passerine_colourfulness/")

library(pavo)
library(sf)
library(parallel)
library(RColorBrewer)
library(letsR)
library(viridis)
library(dispRity)
library(spatstat)
library(maptools)
library(rgeos)
library(cetcolor)

# ------ #

## Custom mapping functions

sp_trait <- function(presab_dat, trait_dat, traits, type=c("average","q.richness","q.proportion"), logtrait=FALSE, geomean=FALSE, range_weights=FALSE, quantile=c(1,2,3,4), base_raster=NULL, plot_map=TRUE, clip_by_richness=0, cores=1) {
  
  # Libraries
  require(sf)
  require(fasterize)
  require(dplyr)
  require(raster)
  require(dispRity)
  require(ggplot2)
  
  # Create an empty raster
  base_raster2 <- base_raster
  base_raster2@data@inmemory <- TRUE
  
  myFun <- function(i) {
    
    tmp_val <- NA
    if (rowSums(presab_dat[i,,drop=F], na.rm=T) > 0) {
      sp_index <- !is.na(presab_dat[i,])
      tmp_nms <- colnames(presab_dat)[sp_index]
      idx <- which(trait_dat$binomial %in% tmp_nms) 
      if (length(idx)>0) {
        trait_vec <- as.matrix(trait_dat[idx, traits], ncol=length(traits))
        if (type == "average") {
          if (logtrait) {
            trait_vec <- log(trait_vec)
          }
          if (geomean) {
            if (range_weights) {
              range_vec <- as.matrix(trait_dat[idx, "ncells"], ncol=1)
              tmp_val <- weighted.geomean(trait_vec, 1/range_vec, na.rm=T)
            } else {
              tmp_val <- geomean(trait_vec, na.rm=TRUE)
            }
          } else {
            if (range_weights) {
              range_vec <- as.matrix(trait_dat[idx, "ncells"], ncol=1)
              tmp_val <- weighted.mean(trait_vec, 1/range_vec, na.rm=T)
            } else {
              tmp_val <- mean(trait_vec, na.rm=TRUE)
            }
          }
        }
        if (type == "q.richness") {
          tmp_val <- sum(trait_vec == quantile)
        }
        if (type == "q.proportion") {
          tmp_val <- sum(trait_vec == quantile) / nrow(trait_vec)
        }
      } else {
        tmp_val <- NA
      }
    }
    return(tmp_val)
  }
  
  output <- mclapply(1:dim(presab_dat)[1], myFun, mc.cores=cores) 
  
  tr_val <- unlist(output)
  
  base_raster2@data@values <- tr_val
  base_raster2@data@min <- min(tr_val, na.rm=TRUE)
  base_raster2@data@max <- max(tr_val, na.rm=TRUE)
  
  if (clip_by_richness > 0) {
    rich <- rowSums(presab_dat)
    values(base_raster2)[which(rich<clip_by_richness)] <- NA
  }
  
  if (plot_map==TRUE) {
    plot(base_raster2)
    #map(interior=FALSE, add=TRUE, lwd=0.5)
    plot(world_lines, add=T, col="black")
  }
  
  return(base_raster2)
  
}

# ----- #

color.legend2 <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
                           gradient = "x", border, ...) 
{
  oldcex <- par("cex")
  par(xpd = TRUE, cex = cex)
  gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
                gradient = gradient, border=border)
  if (gradient == "x") {
    xsqueeze <- (xr - xl)/(2 * length(rect.col))
    textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      texty <- yb - 0.2 * strheight("O")
      textadj <- c(0.5, 1)
    }
    else {
      texty <- yt + 0.2 * strheight("O")
      textadj <- c(0.5, 0)
    }
  }
  else {
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      textx <- xr + 0.2 * strwidth("O")
      textadj <- c(0, 0.5)
    }
    else {
      textx <- xl - 0.2 * strwidth("O")
      textadj <- c(1, 0.5)
    }
  }
  text(textx, texty, labels = legend, adj = textadj, ...)
  par(xpd = FALSE, cex = oldcex)
}

# ----- #

plot_bird_raster <- function(raster_data=NULL, scale="raw", nbins=20, breaks = NULL, col=viridis(20), cex=1.5, legend_pos=c(-8000000,-8000000, 8000000,-7500000), plot_world=TRUE, newproj=NULL){
  
  require(plotrix)
  require(ggplot2)
  require(viridis)
  
  if(scale=="eq_number") { 
    cuts <- cut_number(raster_data@data@values, n = nbins)
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="eq_width") {
    cuts <- cut_interval(raster_data@data@values, n = nbins)
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="raw") {
    cuts <- cut_interval(raster_data@data@values, n = length(na.omit(unique(raster_data@data@values))))
    raster_data@data@values <- as.numeric(cuts)
  }
  if(scale=="predefined") {
    cuts <- cut(raster_data@data@values, breaks = breaks, include.lowest = T)
    raster_data@data@values <- as.numeric(cuts)
  }
  
  if (scale=="predefined") {
    brks <- round(breaks[seq(1,length(breaks),length.out = 6)],2)
  } else {
    cuts_vec <- sort(na.omit(unique(cuts)))
    idx <- c(round(quantile(1:length(cuts_vec), c(0,0.2,0.4,0.6,0.8,1))))
    brks <- as.character(cuts_vec[idx])
    brks <- gsub("\\[|\\]|\\)|\\(", "", brks)
    brks <- as.numeric(unlist(strsplit(brks, ",")))[c(1,3,5,7,9,12)]
  }
  
  if (!is.null(newproj)) {
    fworld <- st_transform(world, crs = newproj)
    fraster_data <- projectRaster(raster_data, crs = newproj)
  } else {
    fworld <- world
    fraster_data <- raster_data
  }
  
  if (plot_world==TRUE) {
    plot(fworld$geometry, border="grey80", col="grey80")
    par(bty="n")
    plot(fraster_data, add=T, legend=F, col=col, axes=FALSE)
  }
  if (plot_world==FALSE) {par(bty="n")
    plot(fworld$geometry, border=NA, col=NA)
    par(bty="n")
    plot(fraster_data, add=T, legend=F, col=col, axes=FALSE)
  }
  
  color.legend2(legend_pos[1], legend_pos[2], legend_pos[3], legend_pos[4], brks, rect.col=col, bty="n", align="rb", cex=cex, border=NA)
  
}

# ------ #

## Get data
res <- read.csv("Outputs/02_Passerine_MF_colourfulness_scores_means.csv", strings=F)
pam <- readRDS("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/PAMs/PAM_birds_Behrman50km_Pres12_Orig12_Seas12_liberal_clipped.rds")
pam <- data.frame(pam)

## Spatial parameters
null_rast <- pam[[2]]
pam <- pam[[1]]
colnames(pam) <- gsub(" ", "_", colnames(pam))
units = 100 # in km
newproj = as.character(crs(null_rast)) # Behrmann cylindrical equal area projection (standard parallels at 30 deg)

## Read in and project world (maptools: wrld_simpl)
data(wrld_simpl)
world = st_as_sf(wrld_simpl)
world = st_transform(world, crs = newproj)
world_lines = st_cast(world, to = "MULTILINESTRING")

# ------ #

### Identify top 25% colourful species based on sex-specific percentile assignments
res$volume.uv.quantile <- as.numeric(res$volume.uv.percentile > 75)
res$volume.vs.quantile <- as.numeric(res$volume.vs.percentile > 75)
res$loci.uv.quantile <- as.numeric(res$loci.uv.percentile > 75)
res$loci.vs.quantile <- as.numeric(res$loci.vs.percentile > 75)

# ------------------- #

### Global inputs
presab_dat <- pam
presab_dat <- presab_dat[,colnames(presab_dat) %in% res$binomial]
cores <- 20
trait_dat <- res
trait_dat <- trait_dat[trait_dat$binomial %in% colnames(presab_dat),]
base_raster <- null_rast
clip_by_richness <- 0
plot_map <- FALSE
logtrait <- FALSE
geomean <- FALSE

# ------ #

## Randomise the trait data (species colourfulness values) wrt species identity

random.res <- as.list(rep(NA,200))

for (i in 1:length(random.res)) {
  
  print (i)
  
  rand_spp_ids <- sample(trait_dat$binomial[trait_dat$sex=="M"], replace = F)
  rand_spp_ids <- c(rand_spp_ids, rand_spp_ids) # as Ms followed by Fs in the dataset
  rand_trait_dat <- trait_dat
  rand_trait_dat$binomial <- rand_spp_ids # overwrite true spp ids with randomised ids to result in randomised trait data associated with a given list of species in an assemblage
  
  ## Calculate cell averages using different approaches and datasets
  res.cells <- data.frame(id=c(1:ncell(null_rast)), xyFromCell(null_rast, cell=c(1:ncell(null_rast))), nspp=rowSums(presab_dat, na.rm=T))
  colnames(res.cells) <- c("id","long","lat","nspp")
  
  # Mean - not range weighted
  range_weights <- FALSE
  res.cells$volume.uv.m.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.uv.f.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.m.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.f.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.m.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.f.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.m.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.f.mean <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  
  # Mean - range weighted
  range_weights <- TRUE
  res.cells$volume.uv.m.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.uv.f.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.m.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.f.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.uv", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.m.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.f.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.m.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.f.mean.wtd <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.vs", "average", logtrait, geomean, range_weights, quantile=NA, base_raster, plot_map, clip_by_richness, cores)@data@values
  
  # Top 25% proportion
  range_weights <- FALSE
  res.cells$volume.uv.m.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.uv.f.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.m.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.uv.f.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.uv.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.m.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "volume.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$volume.vs.f.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "volume.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.m.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="M",], "loci.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  res.cells$loci.vs.f.top25.prop <- sp_trait(presab_dat, rand_trait_dat[rand_trait_dat$sex=="F",], "loci.vs.quantile", "q.proportion", logtrait, geomean, range_weights, quantile=1, base_raster, plot_map, clip_by_richness, cores)@data@values
  
  # Add to list
  random.res[[i]] <- res.cells
  
}

# Save (every iteration)
saveRDS(random.res, "Outputs/05.1_Randomised_cell_assemblage_colourfulness_values_200.rds")

# ------------------------------------------- #

## COMPILE SIMULATED SCORES FOR EACH ECOREGION

ecor.poly.list <- as.list(rep(NA, length(res.cells.list)))

for (xi in 1:length(res.cells.list)) {
  
  # Get individual randomised cells dataset
  
  print (xi)
  res.cells <- res.cells.list[[xi]]
  
  # Summarise colourfulness by ecoregion
  
  ecor.cells <- readRDS("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Outputs/50km/Behrmann_cea/WWF_ecoregions.rds")
  ecor.poly <- read_sf("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Inputs/official_teow/wwf_terr_ecos.shp")
  ecor.lookup <- read.csv("/Volumes/GoogleDrive/My Drive/CRCStorage/Datasets/SpatialData/BirdLife/BirdLife_Shapefiles_GHT/Layers/WWF_Ecoregions/Outputs/50km/Behrmann_cea/WWF_ecoregions_looukp.csv", strings=F)
  
  ecor.poly <- ecor.poly[,c("ECO_ID")]
  ecor.poly <- aggregate(ecor.poly, by=list(ecoregion=ecor.poly$ECO_ID), FUN=mean, do_union=T); ecor.poly$ECO_ID <- NULL
  ecor.poly = st_transform(ecor.poly, crs = newproj)
  
  pam$ecoregion <- ecor.cells$ecoregion_combined
  res.cells$ecoregion <- ecor.cells$ecoregion_combined
  
  fvars <- colnames(res.cells)[-(1:4)]
  
  # -------- #
  
  # Calculate ecoregion means
  
  ecor.means <- aggregate(cbind(volume.uv.m.mean, volume.uv.f.mean, 
                                loci.uv.m.mean, loci.uv.f.mean,
                                volume.vs.m.mean, volume.vs.f.mean, 
                                loci.vs.m.mean, loci.vs.f.mean,
                                volume.uv.m.mean.wtd, volume.uv.f.mean.wtd, 
                                loci.uv.m.mean.wtd, loci.uv.f.mean.wtd,
                                volume.vs.m.mean.wtd, volume.vs.f.mean.wtd, 
                                loci.vs.m.mean.wtd, loci.vs.f.mean.wtd,
                                volume.uv.m.top25.prop, volume.uv.f.top25.prop,
                                loci.uv.m.top25.prop, loci.uv.f.top25.prop,  
                                volume.vs.m.top25.prop, volume.vs.f.top25.prop,
                                loci.vs.m.top25.prop, loci.vs.f.top25.prop) ~ ecoregion, res.cells, mean)
  
  ecor.poly <- ecor.poly[ecor.poly$ecoregion %in% ecor.means$ecoregion,]
  
  ecor.poly$lat <- st_coordinates(st_centroid(ecor.poly))[,"Y"] / units
  ecor.poly$abs.lat <- abs(ecor.poly$lat)
  
  ecor.poly$volume.uv.m.mean <- ecor.means$volume.uv.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.uv.f.mean <- ecor.means$volume.uv.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.m.mean <- ecor.means$loci.uv.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.f.mean <- ecor.means$loci.uv.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$volume.vs.m.mean <- ecor.means$volume.vs.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.vs.f.mean <- ecor.means$volume.vs.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.m.mean <- ecor.means$loci.vs.m.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.f.mean <- ecor.means$loci.vs.f.mean[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$volume.uv.m.mean.wtd <- ecor.means$volume.uv.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.uv.f.mean.wtd <- ecor.means$volume.uv.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.m.mean.wtd <- ecor.means$loci.uv.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.f.mean.wtd <- ecor.means$loci.uv.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$volume.vs.m.mean.wtd <- ecor.means$volume.vs.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.vs.f.mean.wtd <- ecor.means$volume.vs.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.m.mean.wtd <- ecor.means$loci.vs.m.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.f.mean.wtd <- ecor.means$loci.vs.f.mean.wtd[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$volume.uv.m.top25.prop <- ecor.means$volume.uv.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.uv.f.top25.prop <- ecor.means$volume.uv.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.m.top25.prop <- ecor.means$loci.uv.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.uv.f.top25.prop <- ecor.means$loci.uv.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$volume.vs.m.top25.prop <- ecor.means$volume.vs.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$volume.vs.f.top25.prop <- ecor.means$volume.vs.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.m.top25.prop <- ecor.means$loci.vs.m.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  ecor.poly$loci.vs.f.top25.prop <- ecor.means$loci.vs.f.top25.prop[match(ecor.poly$ecoregion, ecor.means$ecoregion)]
  
  ecor.poly$geometry <- NULL
  ecor.poly.list[[xi]] <- ecor.poly
  
}

# ------- #

## CALCULATE PER-ECOREGION SES SCORES

ecor.poly <- readRDS("Outputs/05_Ecoregion_Values.rds")
ecor.poly.tmp <- ecor.poly
ecor.poly.tmp$geometry <- NULL

fit.vars <- c("volume.uv.m.mean", "volume.uv.f.mean", 
              "loci.uv.m.mean", "loci.uv.f.mean",
              "volume.vs.m.mean", "volume.vs.f.mean", 
              "loci.vs.m.mean", "loci.vs.f.mean",
              "volume.uv.m.mean.wtd", "volume.uv.f.mean.wtd", 
              "loci.uv.m.mean.wtd", "loci.uv.f.mean.wtd",
              "volume.vs.m.mean.wtd", "volume.vs.f.mean.wtd", 
              "loci.vs.m.mean.wtd", "loci.vs.f.mean.wtd",
              "volume.uv.m.top25.prop", "volume.uv.f.top25.prop", 
              "loci.uv.m.top25.prop", "loci.uv.f.top25.prop",
              "volume.vs.m.top25.prop", "volume.vs.f.top25.prop", 
              "loci.vs.m.top25.prop", "loci.vs.f.top25.prop")

for (j in 1:length(fit.vars)) {
  
  fvar <- fit.vars[j]
  fobs.dat <- ecor.poly.tmp[,fvar]
  fsim.dat <- c()
  for (xi in 1:length(res.cells.list)) {
    fsim.dat <- cbind(fsim.dat, ecor.poly.list[[xi]][,fvar])
  }
  fses <- (fobs.dat - apply(fsim.dat, 1, mean)) / apply(fsim.dat, 1, sd)
  ecor.poly[,fvar] <- fses
}

# ---------- #

# Run SARs on SES scores

nb <- poly2nb(ecor.poly, row.names = ecor.poly$ecoregion)

fit.vars <- c("volume.uv.m.mean", "volume.uv.f.mean", 
              "loci.uv.m.mean", "loci.uv.f.mean",
              "volume.vs.m.mean", "volume.vs.f.mean", 
              "loci.vs.m.mean", "loci.vs.f.mean",
              "volume.uv.m.mean.wtd", "volume.uv.f.mean.wtd", 
              "loci.uv.m.mean.wtd", "loci.uv.f.mean.wtd",
              "volume.vs.m.mean.wtd", "volume.vs.f.mean.wtd", 
              "loci.vs.m.mean.wtd", "loci.vs.f.mean.wtd",
              "volume.uv.m.top25.prop", "volume.uv.f.top25.prop", 
              "loci.uv.m.top25.prop", "loci.uv.f.top25.prop",
              "volume.vs.m.top25.prop", "volume.vs.f.top25.prop", 
              "loci.vs.m.top25.prop", "loci.vs.f.top25.prop")

out <- as.list(rep(NA, length(fit.vars)))
names(out) <- fit.vars

for (i in 1:length(fit.vars)) {
  
  print (paste(i, "of", length(fit.vars)))
  var <- fit.vars[i]
  
  # fit ols
  fit.ols <- lm(scale(ecor.poly[,var,drop=T]) ~ scale(ecor.poly[,"abs.lat",drop=T]))
  
  # build SAR models with different combinations of distance and neighbourhood style
  neiStyle <- c('W','B','S', 'C', 'U', 'minmax')
  AICvec <- numeric(length = length(neiStyle))
  for (k in 1:length(neiStyle)) {
    nlw <- nb2listw(nb, style=neiStyle[k], zero.policy=TRUE)
    sar_e <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
    AICvec[k] <- AIC(sar_e)			
  }
  bestStyle <- neiStyle[which.min(AICvec)]
  nlw <- nb2listw(nb, style=bestStyle, zero.policy=TRUE)
  fit.sar <- spautolm(fit.ols, listw = nlw, na.action = na.omit, zero.policy = TRUE, tol.solve = 1.0e-25)
  
  mtest1 <- moran.test(resid(fit.ols), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  mtest2 <- moran.test(resid(fit.sar), listw=nlw, na.action = na.omit, zero.policy = TRUE)
  
  aic.ols <- AIC(fit.ols)
  aic.sar <- AIC(fit.sar)
  
  res <- list(fit.ols = fit.ols, fit.sar = fit.sar, aic.ols = aic.ols, aic.sar = aic.sar, moran.ols = mtest1, moran.sar = mtest2)
  
  out[[i]] <- res
  
}

# ----- #

# Summarise

x <- numeric(length(fit.vars))
dff <- data.frame(metric = fit.vars, AIC.ols=x, AIC.sar=x, dAIC.ols=x, dAIC.sar=x,
                  ols.slope=x, ols.pvalue=x,
                  sar.slope=x, sar.pvalue=x, slope.ratio=x,
                  moran.ols = x, moran.ols.pvalue = x, moran.sar = x, moran.sar.pvalue = x, stringsAsFactors=F)

for (i in 1:length(out)) {
  
  fres <- out[[i]]
  
  dff$AIC.ols[i] <- fres$aic.ols
  dff$AIC.sar[i] <- fres$aic.sar
  dff$dAIC.ols[i] <- fres$aic.ols - (min(fres$aic.ols, fres$aic.sar))
  dff$dAIC.sar[i] <- fres$aic.sar - (min(fres$aic.ols, fres$aic.sar))
  dff$ols.slope[i] <- coef(fres$fit.ols)[2]
  dff$ols.pvalue[i] <- summary(fres$fit.ols)$coef[2,4]
  dff$sar.slope[i] <- coef(fres$fit.sar)[2]
  dff$sar.pvalue[i] <- summary(fres$fit.sar)$Coef[2,4]
  dff$slope.ratio[i] <- dff$sar.slope[i] / dff$ols.slope[i]
  dff$moran.ols[i] <- fres$moran.ols$estimate[1]
  dff$moran.ols.pvalue[i] <- fres$moran.ols$p.value
  dff$moran.sar[i] <- fres$moran.sar$estimate[1]
  dff$moran.sar.pvalue[i] <- fres$moran.sar$p.value
  
}


# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #


### Test latitude effect at the species_level


rm(list=ls())

setwd("~/Dropbox/Projects/Current/Bird_colouration/Passerine_colourfulness/")

library(ape)
library(phylolm)
library(MCMCglmm)
library(phytools)
library(geiger)
library(rr2)
library(inlmisc)
library(raster)
library(quantreg)
library(evobiR)
library(brms)
library(parallel)

# ------ #

## Phylogenies
phy2 <- ladderize(read.tree("Data/phylogeny/MCC_tree_stage2.tre"))
phy2 <- rescale(phy2, "depth", 1)

## Tree distributions
trees2 <- readRDS("Data/phylogeny/Tree_distribution_stage2.rds")
trees2 <- lapply(trees2, rescale, model="depth", depth=1)

# ------ #

## read in dataset
dat <- read.csv("Outputs/02_Passerine_MF_colourfulness_scores_means.csv", strings=F)

# ------ #

## separate sexes datasets
datm <- dat[dat$sex=="M",]
datf <- dat[dat$sex=="F",]
rownames(datm) <- datm$binomial
rownames(datf) <- datf$binomial
Ainv2 <- inverseA(phy2, nodes="ALL", scale=T)$Ainv
Ainv2.list <- as.list(rep(NA, length(trees2)))
for (i in 1:length(trees2)) {
  Ainv2.list[[i]] <- inverseA(trees2[[i]], nodes="ALL", scale=T)$Ainv
}

# ------ #

colDens <- function (x, y, colours=GetColors(256)) {
  df <- data.frame(x,y)
  xd <- densCols(x, y, colramp = colorRampPalette(c("black","white")))
  df$dens <- col2rgb(xd)[1,] + 1L
  cols <- colours
  df$col <- cols[df$dens]
  df <- df[order(df$dens),]
  return(df)
}

# ------ #

## Tests

# ------ #

myFun <- function(x, form, dat) {
  fit <- MCMCglmm(form, random=~phylo, ginverse=list(phylo=Ainv2.list[[x]]), data=dat, prior=prior, verbose=TRUE, nitt=110000, burnin=10000, thin=25)
  return(fit)
}

# ------ #

myFun2 <- function (fits) {
  sol <- c()
  vcv <- c()
  for (i in 1:length(fits)) {
    sol <- rbind(sol, fits[[i]]$Sol)
    vcv <- rbind(vcv, fits[[i]]$VCV)
  }
  full <- fits[[1]]
  full$Sol <- as.mcmc(sol)
  full$VCV <- as.mcmc(vcv)
  return(full)
}

# ------ #

myFun3 <- function(x, form, dat) {
  fit <- MCMCglmm(form, random=~phylo, ginverse=list(phylo=sr.Ainv2.list[[x]]), data=dat, prior=prior, verbose=TRUE, nitt=110000, burnin=10000, thin=25)
  return(fit)
}

# ------ #

## phylogenetic heritability
ntrees <- 100
prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

fit <- mclapply(1:ntrees, myFun, form=formula("lat.mean ~ 1"), dat=datm, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("volume.uv ~ 1"), dat=datm, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("volume.uv ~ 1"), dat=datf, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("loci.uv ~ 1"), dat=datm, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("loci.uv ~ 1"), dat=datf, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("volume.vs ~ 1"), dat=datm, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("volume.vs ~ 1"), dat=datf, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("loci.vs ~ 1"), dat=datm, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

fit <- mclapply(1:ntrees, myFun, form=formula("loci.vs ~ 1"), dat=datf, mc.cores=20)
full <- myFun2(fit); MCMCglmmExtras::varComp(full)

# ------ #

## phylogenetic regressions of colourfulness ~ absolute latitude using Bayesian Phylogenetic Mixed Models
ntrees <- 100
prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.uv) ~ scale(lat.abs)"), dat=datm, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.uv) ~ scale(lat.abs)"), dat=datf, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.uv) ~ scale(lat.abs)"), dat=datm, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.uv) ~ scale(lat.abs)"), dat=datf, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.vs) ~ scale(lat.abs)"), dat=datm, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.vs) ~ scale(lat.abs)"), dat=datf, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.vs) ~ scale(lat.abs)"), dat=datm, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.vs) ~ scale(lat.abs)"), dat=datf, mc.cores=20)
full <- myFun2(fit); summary(full)

# ------ #

## phylogenetic regressions of male and female colourfulness ~ absolute latitude * sex using Bayesian Phylogenetic Mixed Models
ntrees <- 100
prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

tmp <- dat
tmp$sex <- ifelse(tmp$sex == "M", 0, 1)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.uv) ~ scale(lat.abs) * scale(sex)"), dat=tmp, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.uv) ~ scale(lat.abs) * scale(sex)"), dat=tmp, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(volume.vs) ~ scale(lat.abs) * scale(sex)"), dat=tmp, mc.cores=20)
full <- myFun2(fit); summary(full)

fit <- mclapply(1:ntrees, myFun, form=formula("scale(loci.vs) ~ scale(lat.abs) * scale(sex)"), dat=tmp, mc.cores=20)
full <- myFun2(fit); summary(full)

# ------------------------ #

## Plots

## Phylogenetic distribution of colourfulness, with clades and latitude

pal <- colorRampPalette(c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00'))(20)
graypal <- gray.colors(20, start = 0.05, end = 0.95, rev = T)

clades <- table(datm$patch.clade)
clades <- clades[clades>20]

png("Plots/phylogenies/loci.uv.png", width = 3, height = 5, units = "in", res = 600)
  par(mar=c(1,1,1,1))
  plot.phylo(phy2, show.tip.label = F, x.lim = c(0,2.1), edge.width = 0.15) # phylogeny c(0,1), data c(1.05,2.05)
  for (i in 1:length(clades)) {
    tip.range <- range(which(phy2$tip.label %in% datm$binomial[datm$patch.clade==names(clades)[i]]))
    lats <- datm$lat.abs[datm$patch.clade==names(clades)[i]]
    ptrop <- sum(lats < 23.5) / length(lats)
    pcol <- graypal[as.numeric(cut(ptrop, breaks=seq(0,1,0.05)))]
    polygon(x=c(1.025,1.05,1.05,1.025), y=rep(tip.range, each=2), border = NA, col= pcol)
  }
  # loci.uv
  breaks <- quantile(dat$loci.uv, probs=seq(0,1,length.out=length(pal)+1))
  mcols <- pal[as.numeric(cut(datm$loci.uv, breaks, include.lowest=T))]
  fcols <- pal[as.numeric(cut(datf$loci.uv, breaks, include.lowest=T))]
  names(mcols) <- datm$binomial
  names(fcols) <- datf$binomial
  exts <- c(1.05,1.3,1.55)
  frescaled <- scales:::rescale(datf$loci.uv, from=range(dat$loci.uv), to=c(exts[2],exts[1]))
  mrescaled <- scales:::rescale(datm$loci.uv, from=range(dat$loci.uv), to=c(exts[2],exts[3]))
  names(frescaled) <- datf$binomial
  names(mrescaled) <- datm$binomial
  for (i in 1:Ntip(phy2)) {
    ftip <- phy2$tip.label[i]
    lines(c(exts[2], frescaled[ftip]), c(i,i), col=fcols[ftip], lwd=0.25)
    lines(c(exts[2], mrescaled[ftip]), c(i,i), col=mcols[ftip], lwd=0.25)
  }
  lines(c(exts[2],exts[2]), c(-10, Ntip(phy2)+10), lwd=1.5)
  # scales
  points(x=seq(0.1,0.4,length.out=length(pal)), y=rep(400,length(pal)), pch=15, col=pal, cex=0.5)
  text(x=mean(seq(0.1,0.4,length.out=5)), rep(290,5), labels="Colour loci", cex=0.3, font=2)
  text(x=seq(0.1,0.4,length.out=5), rep(350,5), labels=round(10^breaks[seq(1, length(breaks), length.out = 5)]), cex=0.25)
  points(x=seq(0.1,0.4,length.out=length(graypal)), y=rep(200,length(graypal)), pch=15, col=graypal, cex=0.5)
  text(x=mean(seq(0.1,0.4,length.out=5)), rep(90,5), labels="Proportion tropical", cex=0.3, font=2)
  text(x=seq(0.1,0.4,length.out=5), rep(150,5), labels=seq(0,1,length.out=5), cex=0.25)
dev.off()

# ------ #

## Latitude boxplots

png("Plots/species_level/Boxplot_absolute_latitude_uv.png", height=6, width=6, res=600, units="in")
  par(mfcol=c(2,2), mar=c(4,4,2,2))
  xcentre <- zoo:::rollmean(seq(0,70,5),2)
  xnudge <- 2
  bx <- boxplot(datm$volume.uv ~ cut(datm$lat.abs, breaks=seq(0,70,5)), outline=F, plot=F)$stats[c(2,3,4),] # IQR
  plot(0,0, xlab="", ylab="", las=1, bty="n", cex=0.5, xlim=c(0,70), ylim=range(bx), axes=F)
  axis(1, at = seq(0,70,10), cex.axis=0.75)
  axis(2, at = seq(range(bx)[1],range(bx)[2],length.out=6), labels = format(round(10^seq(range(bx)[1],range(bx)[2],length.out=6)*1000,2), nsmall=2), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Midpoint latitude (°, absolute)", line=2.5, cex=0.75)
  mtext(side=2, text="Colour volume", line=3, cex=0.75)
  for (i in 1:ncol(bx)) {
    polygon(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge, xcentre[i]+xnudge, xcentre[i]-xnudge), y=c(bx[1,i],bx[1,i],bx[3,i],bx[3,i]), col="grey90")
    lines(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge), y=c(bx[2,i],bx[2,i]))
  }
  bx <- boxplot(datf$volume.uv ~ cut(datf$lat.abs, breaks=seq(0,70,5)), outline=F, plot=F)$stats[c(2,3,4),] # IQR
  plot(0,0, xlab="", ylab="", las=1, bty="n", cex=0.5, xlim=c(0,70), ylim=range(bx), axes=F)
  axis(1, at = seq(0,70,10), cex.axis=0.75)
  axis(2, at = seq(range(bx)[1],range(bx)[2],length.out=6), labels = format(round(10^seq(range(bx)[1],range(bx)[2],length.out=6)*1000,2), nsmall=2), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Midpoint latitude (°, absolute)", line=2.5, cex=0.75)
  mtext(side=2, text="Colour volume", line=3, cex=0.75)
  for (i in 1:ncol(bx)) {
    polygon(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge, xcentre[i]+xnudge, xcentre[i]-xnudge), y=c(bx[1,i],bx[1,i],bx[3,i],bx[3,i]), col="grey90")
    lines(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge), y=c(bx[2,i],bx[2,i]))
  }
  bx <- boxplot(datm$loci.uv ~ cut(datm$lat.abs, breaks=seq(0,70,5)), outline=F, plot=F)$stats[c(2,3,4),] # IQR
  plot(0,0, xlab="", ylab="", las=1, bty="n", cex=0.5, xlim=c(0,70), ylim=range(bx), axes=F)
  axis(1, at = seq(0,70,10), cex.axis=0.75)
  axis(2, at = seq(range(bx)[1],range(bx)[2],length.out=6), labels = round(10^seq(range(bx)[1],range(bx)[2],length.out=6)), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Midpoint latitude (°, absolute)", line=2.5, cex=0.75)
  mtext(side=2, text="Colour loci", line=2.5, cex=0.75)
  for (i in 1:ncol(bx)) {
    polygon(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge, xcentre[i]+xnudge, xcentre[i]-xnudge), y=c(bx[1,i],bx[1,i],bx[3,i],bx[3,i]), col="grey90")
    lines(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge), y=c(bx[2,i],bx[2,i]))
  }
  bx <- boxplot(datf$loci.uv ~ cut(datf$lat.abs, breaks=seq(0,70,5)), outline=F, plot=F)$stats[c(2,3,4),] # IQR
  plot(0,0, xlab="", ylab="", las=1, bty="n", cex=0.5, xlim=c(0,70), ylim=range(bx), axes=F)
  axis(1, at = seq(0,70,10), cex.axis=0.75)
  axis(2, at = seq(range(bx)[1],range(bx)[2],length.out=6), labels = round(10^seq(range(bx)[1],range(bx)[2],length.out=6)), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Midpoint latitude (°, absolute)", line=2.5, cex=0.75)
  mtext(side=2, text="Colour loci", line=2.5, cex=0.75)
  for (i in 1:ncol(bx)) {
    polygon(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge, xcentre[i]+xnudge, xcentre[i]-xnudge), y=c(bx[1,i],bx[1,i],bx[3,i],bx[3,i]), col="grey90")
    lines(x=c(xcentre[i]-xnudge, xcentre[i]+xnudge), y=c(bx[2,i],bx[2,i]))
  }
dev.off()

# ------ #

## Correlation between the sexes

cols <- GetColors(256, scheme="smooth rainbow", stops=c(0.1,1), reverse=F)

png("Plots/species_level/Scatterplot_male_female_colourfulness.png", height=6, width=6, res=600, units="in")
  par(mfrow=c(2,2), mar=c(4,4,2,2))
  df <- colDens(x = datf$volume.uv, y = datm$volume.uv, colours = cols); rownames(df) <- datm$binomial
  plot(y ~ x, df, col=col, las=1, xlim=range(dat$volume.uv), ylim=range(dat$volume.uv), pch=16, cex=0.5, xlab="", ylab="", bty="n", axes=F); abline(0,1, lty=2)
  axis(1, at = seq(range(dat$volume.uv)[1],range(dat$volume.uv)[2],length.out=6), labels = format(round(10^seq(range(dat$volume.uv)[1],range(dat$volume.uv)[2],length.out=6)*1000,2), nsmall=2),  cex.axis=0.75)
  axis(2, at = seq(range(dat$volume.uv)[1],range(dat$volume.uv)[2],length.out=6), labels = format(round(10^seq(range(dat$volume.uv)[1],range(dat$volume.uv)[2],length.out=6)*1000,2), nsmall=2), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Female colour volume", line=2.5, cex=0.75)
  mtext(side=2, text="Male colour volume", line=3, cex=0.75)
  x <- df$x; names(x) <- rownames(df); y <- df$y; names(y) <- rownames(df); xy <- as.matrix(cbind(x,y))
  fit1 <- phyl.RMA(x, y, phy2, method="lambda"); fit1
  lines(fit1$RMA.beta[1] + (fit1$RMA.beta[2] * range(x)) ~ range(x), col="black", lwd=2)
  df <- colDens(x = datf$volume.vs, y = datm$volume.vs, colours = cols); rownames(df) <- datm$binomial
  plot(y ~ x, df, col=col, las=1, xlim=range(dat$volume.vs), ylim=range(dat$volume.vs), pch=16, cex=0.5, xlab="", ylab="", bty="n", axes=F); abline(0,1, lty=2)
  axis(1, at = seq(range(dat$volume.vs)[1],range(dat$volume.vs)[2],length.out=6), labels = format(round(10^seq(range(dat$volume.vs)[1],range(dat$volume.vs)[2],length.out=6)*1000,2), nsmall=2),  cex.axis=0.75)
  axis(2, at = seq(range(dat$volume.vs)[1],range(dat$volume.vs)[2],length.out=6), labels = format(round(10^seq(range(dat$volume.vs)[1],range(dat$volume.vs)[2],length.out=6)*1000,2), nsmall=2), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Female colour volume", line=2.5, cex=0.75)
  mtext(side=2, text="Male colour volume", line=3, cex=0.75)
  x <- df$x; names(x) <- rownames(df); y <- df$y; names(y) <- rownames(df); xy <- as.matrix(cbind(x,y))
  fit2 <- phyl.RMA(x, y, phy2, method="lambda"); fit2
  lines(fit2$RMA.beta[1] + (fit2$RMA.beta[2] * range(x)) ~ range(x), col="black", lwd=2)
  df <- colDens(x = datf$loci.uv, y = datm$loci.uv, colours = cols); rownames(df) <- datm$binomial
  plot(y ~ x, df, col=col, las=1, xlim=range(dat$loci.uv), ylim=range(dat$loci.uv), pch=16, cex=0.5, xlab="", ylab="", bty="n", axes=F); abline(0,1, lty=2)
  axis(1, at = seq(range(dat$loci.uv)[1],range(dat$loci.uv)[2],length.out=6), labels = round(10^seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6)),  cex.axis=0.75)
  axis(2, at = seq(range(dat$loci.uv)[1],range(dat$loci.uv)[2],length.out=6), labels = round(10^seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6)), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Female colour loci", line=2.5, cex=0.75)
  mtext(side=2, text="Male colour loci", line=3, cex=0.75)
  x <- df$x; names(x) <- rownames(df); y <- df$y; names(y) <- rownames(df); xy <- as.matrix(cbind(x,y))
  fit3 <- phyl.RMA(x, y, phy2, method="lambda"); fit3
  lines(fit3$RMA.beta[1] + (fit3$RMA.beta[2] * range(x)) ~ range(x), col="black", lwd=2)
  df <- colDens(x = datf$loci.vs, y = datm$loci.vs, colours = cols); rownames(df) <- datm$binomial
  plot(y ~ x, df, col=col, las=1, xlim=range(dat$loci.vs), ylim=range(dat$loci.vs), pch=16, cex=0.5, xlab="", ylab="", bty="n", axes=F); abline(0,1, lty=2)
  axis(1, at = seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6), labels = round(10^seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6)),  cex.axis=0.75)
  axis(2, at = seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6), labels = round(10^seq(range(dat$loci.vs)[1],range(dat$loci.vs)[2],length.out=6)), las=1, line = 0, cex.axis=0.75)
  mtext(side=1, text="Female colour loci", line=2.5, cex=0.75)
  mtext(side=2, text="Male colour loci", line=3, cex=0.75)
  x <- df$x; names(x) <- rownames(df); y <- df$y; names(y) <- rownames(df); xy <- as.matrix(cbind(x,y))
  fit4 <- phyl.RMA(x, y, phy2, method="lambda"); fit4
  lines(fit4$RMA.beta[1] + (fit4$RMA.beta[2] * range(x)) ~ range(x), col="black", lwd=2)
dev.off()


# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #


### Multipredictor models of colourfulness scores

rm(list=ls())

setwd("~/Dropbox/Projects/Current/Bird_colouration/Passerine_colourfulness/")

library(ape)
library(phylolm)
library(MCMCglmm)
library(phytools)
library(geiger)
library(rr2)
library(inlmisc)
library(raster)
library(quantreg)
library(evobiR)
library(brms)
library(car)
library(parallel)
library(MuMIn)

# ------ #

## Load data
load("Outputs/08_Predictor_model_dataset.RData") # trees, dataset (vars scaled and centred)
multipred <- c("BodyMass","FruitNect","Migratory","GroundNest","Territorial","Temp","Precip","SolarRad","UVB","NPP","ForestDep","CoexistSpp","Dichro")
resp <- c("volume.uv","loci.uv","volume.vs","loci.vs")

# ------ #

Ainv.list <- as.list(rep(NA, length(trees2)))
for (i in 1:length(trees2)) {
  Ainv.list[[i]] <- inverseA(trees2[[i]], nodes="ALL", scale=T)$Ainv
}

# ------ #

## Tree distribution
## MCMCGLMM - lambda

myFun <- function(x, form, dat) {
  fit <- MCMCglmm(form, random=~phylo, ginverse=list(phylo=Ainv.list[[x]]), data=dat, prior=prior, verbose=TRUE, nitt=110000, burnin=10000, thin=25)
  return(fit)
}

prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
ntrees <- 100

for (i in 1:length(resp)) {
  print(resp[i])
  fresp <- resp[i]
  fform <- formula(paste(fresp, "~", paste(multipred, collapse="+")))
  fit.list.m <- NULL
  fit.list.f <- NULL
  fit.list.m <- mclapply(1:ntrees, myFun, form=fform, dat=datm, mc.cores=20)
  saveRDS(fit.list.m, paste0("Outputs/09_Predictor_model_mcmcglmm/", fresp, ".m.mcmcglmm.rds"))
  fit.list.f <- mclapply(1:ntrees, myFun, form=fform, dat=datf, mc.cores=20)
  saveRDS(fit.list.f, paste0("Outputs/09_Predictor_model_mcmcglmm/", fresp, ".f.mcmcglmm.rds"))
}

# ------- #

myFun <- function (fits) { # phylo
  sol <- c()
  vcv <- c()
  for (i in 1:length(fits)) {
    sol <- rbind(sol, fits[[i]]$Sol)
    vcv <- rbind(vcv, fits[[i]]$VCV)
  }
  full <- fits[[1]]
  full$Sol <- as.mcmc(sol)
  full$VCV <- as.mcmc(vcv)
  full.out <- summary(full)$solutions
  full.out <- cbind(full.out, t(apply(full$Sol, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975))))
  full.out <- data.frame(full.out, print.est=paste(format(round(full.out[,"post.mean"],3), nsmall=3), " (", format(round(full.out[,"l-95% CI"],3), nsmall=3), ", ", format(round(full.out[,"u-95% CI"],3), nsmall=3), ")", sep=""))
  full.out <- data.frame(full.out, print.p=format(round(full.out[,"pMCMC"],3), nsmall=3))
  return(full.out)
}

# ------- #

porder <- c("(Intercept)",
            "Temp","Precip","NPP",
            "SolarRad","UVB",
            "BodyMass",
            "Dichro",
            "ForestDep","FruitNect",
            "Migratory","GroundNest","Territorial",
            "CoexistSpp")

# ------- #

## results tables

fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/volume.uv.m.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/volume.uv.m.mcmcglmm.coefs.csv")
fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/volume.uv.f.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/volume.uv.f.mcmcglmm.coefs.csv")

fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/loci.uv.m.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/loci.uv.m.mcmcglmm.coefs.csv")
fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/loci.uv.f.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/loci.uv.f.mcmcglmm.coefs.csv")

fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/volume.vs.m.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/volume.vs.m.mcmcglmm.coefs.csv")
fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/volume.vs.f.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/volume.vs.f.mcmcglmm.coefs.csv")

fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/loci.vs.m.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/loci.vs.m.mcmcglmm.coefs.csv")
fits <- readRDS("Outputs/09_Predictor_model_mcmcglmm/loci.vs.f.mcmcglmm.rds"); coefs <- myFun(fits)[porder,]; coefs
write.csv(coefs, "Outputs/09_Predictor_model_mcmcglmm/loci.vs.f.mcmcglmm.coefs.csv")

# ------- #

## plot

library(phylolm)
library(ggplot2)
library(inlmisc)

# ------ #

porder <- c("Temp","Precip","NPP",
            "SolarRad","UVB",
            "BodyMass",
            "Dichro",
            "ForestDep","FruitNect",
            "Migratory","GroundNest","Territorial",
            "CoexistSpp")

pcols <- rep(c("#fec44f","#fe9929","#ef3b2c","#4292c6","#a1d99b","#d4b9da","#df65b0"), c(3,2,1,1,2,3,1))

# ------ #

# loci.uv - mcmcglmm

mc <- read.csv("Outputs/09_Predictor_model_mcmcglmm/loci.uv.m.mcmcglmm.coefs.csv", row.names = 1)[-1,][porder,]
fc <- read.csv("Outputs/09_Predictor_model_mcmcglmm/loci.uv.f.mcmcglmm.coefs.csv", row.names = 1)[-1,][porder,]

pdf("Plots/predictor_models/mcmcglmm_loci_uv_coefs.pdf", height=4, width=5)
par(mfrow=c(1,2), mar=c(2,2,1,1))
dg <- 0.5
plot(0,0, xlim=c(-0.3,0.3), ylim=c(length(multipred)*2,1), type="n", bty="l", xlab="", ylab="", axes=F)
axis(1, at = round(seq(-0.3,0.3,0.1),2), labels = NA, cex.axis=0.7)
axis(2, at = c(1:length(multipred)*2), labels = NA)
mtext(round(seq(-0.3,0.3,0.1),2), side=1, at=round(seq(-0.3,0.3,0.1),2), line=0.5, cex=0.7)
abline(v=0, lty=2, col="grey80")
box(bty="l")
for (i in 1:length(porder)) {
  icol <- ifelse(mc$pMCMC[i] < 0.05, "black", "grey80")
  ilwd <- ifelse(mc$pMCMC[i] < 0.05, 1, 1)
  lines(mc[i,c("l.95..CI","u.95..CI")], c(i,i)*2, col=icol, lwd=ilwd)
  polygon(c(mc$X25.[i], mc$X75.[i], mc$X75.[i], mc$X25.[i]), c(i*2-dg, i*2-dg, i*2+dg, i*2+dg), col=pcols[i], border=icol, lwd=ilwd)
  lines(c(mc$X50.[i], mc$X50.[i]), c(i*2-dg,i*2+dg), lend="butt", col=icol, lwd=ilwd)
  if (mc$pMCMC[i] < 0.05 & mc$pMCMC[i] >= 0.01) {
    text(x=mc$X50.[i], y=i*2-dg*2, labels = "*", adj = 0.5)
  }
  if (mc$pMCMC[i] < 0.01 & mc$pMCMC[i] >= 0.001) {
    text(x=mc$X50.[i], y=i*2-dg*2, labels = "**", adj = 0.5)
  }
  if (mc$pMCMC[i] < 0.001) {
    text(x=mc$X50.[i], y=i*2-dg*2, labels = "***", adj = 0.5)
  }
}
plot(0,0, xlim=c(-0.3,0.3), ylim=c(length(multipred)*2,1), type="n", bty="l", xlab="", ylab="", axes=F)
axis(1, at = round(seq(-0.3,0.3,0.1),2), labels = NA, cex.axis=0.7)
axis(2, at = c(1:length(multipred)*2), labels = NA)
mtext(round(seq(-0.3,0.3,0.1),2), side=1, at=round(seq(-0.3,0.3,0.1),2), line=0.5, cex=0.7)
abline(v=0, lty=2, col="grey80")
box(bty="l")
for (i in 1:length(porder)) {
  icol <- ifelse(fc$pMCMC[i] < 0.05, "black", "grey80")
  ilwd <- ifelse(fc$pMCMC[i] < 0.05, 1, 1)
  lines(fc[i,c("l.95..CI","u.95..CI")], c(i,i)*2, col=icol, lwd=ilwd)
  polygon(c(fc$X25.[i], fc$X75.[i], fc$X75.[i], fc$X25.[i]), c(i*2-dg, i*2-dg, i*2+dg, i*2+dg), col=pcols[i], border=icol, lwd=ilwd)
  lines(c(fc$X50.[i], fc$X50.[i]), c(i*2-dg,i*2+dg), lend="butt", col=icol, lwd=ilwd)
  if (fc$pMCMC[i] < 0.05 & fc$pMCMC[i] >= 0.01) {
    text(x=fc$X50.[i], y=i*2-dg*2, labels = "*", adj = 0.5)
  }
  if (fc$pMCMC[i] < 0.01 & fc$pMCMC[i] >= 0.001) {
    text(x=fc$X50.[i], y=i*2-dg*2, labels = "**", adj = 0.5)
  }
  if (fc$pMCMC[i] < 0.001) {
    text(x=fc$X50.[i], y=i*2-dg*2, labels = "***", adj = 0.5)
  }
}
legend("right", legend=c("Climate","Irradiance","Body size","Sexual sel.","Ecology","Behaviour","Diversity"), bty="n", pch=15, col=unique(pcols), cex=0.4, pt.cex=0.7)
dev.off()


# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #



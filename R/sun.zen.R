sun.zen = function(ras,time, crs.code = CRS("+init=epsg:4326"),...){

	ras.grid = suppressWarnings(as(ras,"SpatialGridDataFrame"))
	gc()
	ras.trans = suppressWarnings(spTransform(ras.grid, CRS = crs.code))
	lat = mean(ras.trans@coords[,2])
	lon = mean(ras.trans@coords[,1])
	act.tm = solar(as.POSIXct(time))
	zen = suppressWarnings(zenith(act.tm, lat, lon))
	gc()
	rm(lat,lon,ras.grid)
	return(zen)

invisible(y)
}
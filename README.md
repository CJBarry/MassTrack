# MassTrack

Track Particle Mass along Calculated Pathlines in Groundwater

Designed as a post-processor to MODPATH, this package is a
model for tracking particle mass along a calculated set of pathlines by
analysing the flow balance of each cell through which the particles
pass.  The code takes advantage of quick subsetting with the data.table
package and NetCDF format for MODFLOW data provided by the Rflow
package.

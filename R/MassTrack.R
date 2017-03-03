# Christopher Barry, started on 11/01/2016 at University of Birmingham

# traces path lines through MODFLOW cells and tracks how much mass is lost on
#  the way, assuming complete mixing in each cell

# For each cell that the path passes through, the imbalance in flows entering
#  and exiting the cell faces and into or out of storage is calculated. If there
#  is a net loss, then the particle mass is reduced proportionately.  If there is
#  a net gain, the particle mass is unchanged. When the particle arrives at a
#  cell from which no water is exitting (i.e. a strong sink), its mass will be
#  reduced to zero (given enough time) and it will be considered captured.

# In reality, mass entering storage will be delayed.
# Without releasing more particles, one must either neglect the delay and retain
#  the mass, or consider the mass lost and ignore delayed release from storage.

# Particle mass set to 1 at release.

# MODPATH should be run with very lenient weak-sink options.

# constant porosity value only supported currently

# the headers that apply to a pathline data table from MODPATH-5
PTL.headers <- c("ptlno", "x", "y", "z_off", "z",
                 "t", "C", "R", "L", "timestep", "m")

#' MassTrack
#'
#' Track particle mass along a MODPATH-5 generated pathline
#'
#' @param ptl a data.frame, data.table or string giving the path to a data
#' file containing a pathline data set produced by MODPATH
#' @param gwnc a NetCDF of MODFLOW simulation data produced by GW.nc, or a
#' path to one
#' @param wtop a NetCDF of water height or path to one (MassTrack can
#' generate one if this hasn't been made)
#' @param phi porosity: single value, layer by layer (vector) or
#' distributed array
#' @param m0 initial mass of particles: single value, vector of values or
#' function of time only
#' @param retain.storage should mass flux to storage be retained on a
#' pathline?
#' @param truncate if TRUE, once a particle is completely drained, further
#' trajectory information is not returned with the result
#' @param confined not yet implemented
#' @param outmass return cell referenced information about where the mass
#' is abstracted?
#' @param loss return information about where particle mass is not released
#' @param outmass.array
#' logical \code{[1]};
#' if outmass = TRUE, then outmass.array = TRUE signals that the abstracted
#'  mass information should be returned in an array rather than a data.table
#' @param tlumpoutmass
#' logical \code{[1]};
#' should the outmass array only show the spatial cell location (TRUE), or
#'  show during which time step the abstractions occurred as well (in which
#'  case the outmass array will be four-dimensional and potentially very
#'  large)
#' @param linear.decay
#' double \code{[1]};
#' linear decay constant, defaults to 0 which is no decay
#' @param react.loss
#' logical \code{[1]};
#' similar to outmass, save locational reaction loss information?
#' @param end.t
#' double \code{[1]};
#' specify the end time of the simulation; if missing, it will be inferred
#'  from ptl, but the last particles to arrive at an abstraction may not be
#'  drained because the simulation will be assumed to have finished as soon
#'  as or just after they arrive
#'
#' @return
#' A list of results:
#' $trace, data.table: the pathline data table with new columns for mass
#'  and, if \code{outmass.array = FALSE}, mass loss and reactive mass
#'  loss\cr
#' $outmass, numeric array: if \code{outmass = TRUE} and
#'  \code{outmass.array = TRUE}, the mass removed from the system by
#'  column, row and layer and, if \code{tlumpoutmass = FALSE}, time step\cr
#' $react.loss, numeric array: as with outmass but for reactive loss\cr
#' $loss, numeric \code{[Np]}: mass which is lost otherwise (e.g. particle
#'  wasn't released)
#'
#' Only requested outputs are returned and if only traces are returned,
#'  then the data table is not put in a list.
#'
#' @import data.table
#' @import RNetCDF
#' @export
#'
#' @examples
#' # MODFLOW NetCDF dataset
#' library("RNetCDF")
#' mfdata <- open.nc(system.file("masstrack_mf_demo.nc",
#'                               package = "MassTrack"))
#'
#' # pathline file
#' library("data.table")
#' ptl <- fread(system.file("masstrack_mp_demo.ptl",
#'                          package = "MassTrack"), skip = 1L)
#'
#' # run MassTrack
#' # - this will generate a new NetCDF for the water top (see wtop
#' #    argument)
#' # - a simple example in which all particles start with mass 1, and there
#' #    is no degradationg
#' mt <- MassTrack(ptl, mfdata, phi = .1, end.t = 1500)
#'
#' # plot results
#' # - plot MODFLOW model boundaries
#' library("Rflow")
#' MFimage(is.na(var.get.nc(mfdata, "Head",
#'                          c(NA, NA, 1L, 1L),
#'                          c(NA, NA, 1L, 1L))),
#'         gccs(mfdata), grcs(mfdata), 0:1,
#'         c("transparent", "grey"), asp = 1)
#'
#' MFimage(var.get.nc(mfdata, "ConstantHead",
#'                    c(NA, NA, 1L, 1L),
#'                    c(NA, NA, 1L, 1L)) != 0,
#'         gccs(mfdata), grcs(mfdata), 0:1,
#'         c("transparent", "blue"), add = TRUE)
#'
#' MFimage(var.get.nc(mfdata, "RiverLeakage",
#'                    c(NA, NA, 1L, 1L),
#'                    c(NA, NA, 1L, 1L)) != 0,
#'         gccs(mfdata), grcs(mfdata), 0:1,
#'         c("transparent", "green"), add = TRUE)
#'
#' MFimage(var.get.nc(mfdata, "Wells",
#'                    c(NA, NA, 1L, 1L),
#'                    c(NA, NA, 1L, 1L)) != 0,
#'         gccs(mfdata), grcs(mfdata), 0:1,
#'         c("transparent", "red"), add = TRUE)
#'
#' # plot pathlines
#' # - note that data.table subsetting often doesn't work with lines
#' ptl <- fread(system.file("masstrack_mp_demo.ptl",
#'                          package = "MassTrack"), skip = 1L)
#' setnames(ptl, c("ptlno", "x", "y", "z_off", "z",
#'                 "t", "C", "R", "L", "timestep"))
#'
#' for(p in unique(ptl$ptlno)) ptl[ptlno == p, lines(x, y)]
#'
#' # plot particle masses
#' mt$traces[, points(x, y, pch = 16L, col = "darkgreen", cex = m)]
#'
#' # you can save the results to file using rlist::list.save, for example
#' \dontrun{
#'   rlist::list.save(mt, "MT_EXAMPLE.rds")
#' }
#'
MassTrack <- function(ptl = "pathline", gwnc, wtop = "wtop.nc",
                      phi, m0 = 1,
                      retain.storage = TRUE, truncate = TRUE,
                      confined = FALSE, outmass = TRUE, loss = TRUE,
                      outmass.array = FALSE, tlumpoutmass = FALSE,
                      linear.decay = FALSE, react.loss = TRUE, end.t){
  gwnc <- switch(class(gwnc),
                 NetCDF = gwnc,
                 character = open.nc(gwnc),
                 stop("gwnc must be a NetCDF dataset or a file path to one"))

  # initialise
  dtits <- c("NCOL", "NROW", "NLAY", "NTS")
  ds <- c(sapply(dtits, dim.inq.nc, ncfile = gwnc)["length",],
          recursive = TRUE)
  spds <- ds[1:3]
  nts <- ds[4L]

  if(outmass && outmass.array){
    out <- array(0, if(tlumpoutmass) spds else ds)
  }
  if(!linear.decay) react.loss <- FALSE

  if(react.loss && outmass.array){
    outr <- array(0, if(tlumpoutmass) spds else ds)
  }

  # get water top, write to new NetCDF if necessary
  if(is.character(wtop) && file.exists(wtop)){
    wtop <- open.nc(wtop)
    on.exit(close.nc(wtop), add = TRUE)
  }else if(is.character(wtop)){
    wtop <- create.nc(wtop, large = TRUE)
    on.exit(close.nc(wtop), add = TRUE)
    att.put.nc(wtop, "NC_GLOBAL", "title", "NC_CHAR",
               "height of saturated groundwater above datum")
    att.put.nc(wtop, "NC_GLOBAL", "explanation", "NC_CHAR",
               "Water top is the water table in unconfined regions and the top of the highest saturated cell in confined regions.")
    att.put.nc(wtop, "NC_GLOBAL", "history", "NC_CHAR",
               paste("Created on", date(), "by MassTrack"))

    Map(dim.def.nc, dimname = dtits, dimlength = ds,
        MoreArgs = list(ncfile = wtop))

    var.def.nc(wtop, "wtop", "NC_FLOAT", dtits)
    att.copy.nc(gwnc, "Head", "missing_value", wtop, "wtop")

    # lt is layer top
    lt <- c(var.get.nc(gwnc, "elev", count = spds))
    for(i in 1:nts){
      var.put.nc(wtop, "wtop", {
        # wt is water head
        wt <- c(var.get.nc(gwnc, "Head", c(1L, 1L, 1L, i), c(spds, 1L)))

        # layer top or water head, whichever is greater
        ifelse(lt > wt, wt, lt)
      }, c(1L, 1L, 1L, i), c(spds, 1L))
    }

    att.copy.nc(gwnc, "Head", "units", wtop, "wtop")
    att.copy.nc(gwnc, "NC_GLOBAL", "datum", wtop, "NC_GLOBAL")
  }
  stopifnot(identical(class(wtop), "NetCDF"))

  if(is.character(ptl)){
    ptl <- fread(ptl, skip = 1L)}
  else if(is.data.frame(ptl)){
    if(!is.data.table(ptl)) setDT(ptl)
  }else stop({
    "argument ptl.df must either be a 10-column data frame/ table or a file name to read as such"
  })
  # stopifnot(is.data.table(ptl))

  # delete unexpected extra columns
  # warnings are given if non-existent columns are set to NULL, but that is
  #  fine in this case
  suppressWarnings(set(ptl, NULL,
                       names(ptl)[!names(ptl) %in%
                                    c(PTL.headers[1:10],
                                      paste0("V", 1:10))], NULL))

  setnames(ptl, PTL.headers[1:10])
  setkey(ptl, ptlno)

  ## which flow parameters to check for imbalances
  # order is important for later (CRLts)
  afps <- c("FlowRightFace", "FlowFrontFace", "FlowLowerFace", "Storage")
  pnm <- var.get.nc(gwnc, "outvars")

  # links flow parameters to a dimension number
  fpdns <- structure(1:4, names = afps)[which({
    c("FlowRightFace", "FlowFrontFace", "FlowLowerFace",
      if(retain.storage) "Storage")
  } %chin% pnm)]

  CRLts <- as.matrix(ptl[, .(C, R, L, timestep)])

  ## exitting and entering flow
  # minimum and maximum C, R, L and ts indices
  minCRLts <- apply(CRLts, 2L, min)
  maxCRLts <- apply(CRLts, 2L, max)
  rgCRLts <- maxCRLts - minCRLts + 1L

  # this allows only the necessary data to be read into memory
  atedge <- minCRLts == 1L
  dstarts <- minCRLts - ifelse(atedge, 0L, 1L)
  dcounts <- rgCRLts + ifelse(atedge, 0L, 1L)

  # modify cell references
  CRLts2 <- vapply(1:4, function(i){
    CRLts[, i] - dstarts[i] + 1L
  }, integer(nrow(CRLts)))

  # all the following are out - in; that is, a positive value indicates
  #  that the cell is being drained
  # change in F{face}F
  # vapply used as faster than mapply (no need to call simplify2array)

  dJ <- vapply(1:length(fpdns), function(i){
    parname <- names(fpdns)[i]
    spino <- fpdns[i]

    # get the subarray from the data set
    F_Far <- var.get.nc(gwnc, parname, dstarts, dcounts, collapse = FALSE)

    # change in flow to the across the cell in the direction being analysed
    dF_F <-{
      # flow through faces to adjacent cells
      Jout <- F_Far[CRLts2]

      # flow through faces from adjacent cells
      if(!identical(fpdns, 4L)){
        imod <- array(0L, dim(CRLts2))

        # NA used to flag edge of model
        # -1 signifies value from adjacent cell flowing in
        imod[, spino] <- ifelse(CRLts2[, spino] == 1L, NA, -1L)

        # takes advantage of R's matrix subsetting (arr[imtx] where imtx
        #  has as many columns as arr has dimensions)
        Jin <- F_Far[CRLts2 + imod]

        # if at the edge of the model, then the flow coming from the
        #  adjacent cell is zero
        Jin[is.na(Jin)] <- 0

      }else{
        # case for analysing Storage
        # negative values signify flow to Storage
        Jout <- -Jout
        Jin <- 0
      }

      Jout - Jin
    }
  }, double(nrow(CRLts2)))

  # a cell with no flows to boundaries ought to be balanced, apart from
  #  machine imprecision
  # if water locked in storage is considered lost, then flows to Storage
  #  will also cause non-zero imbalance values
  # imbalance now considered negative if cell is being drained (net) and
  #  positive if it is being recharged (net)
  imbalance <- rowSums(dJ, na.rm = TRUE)

  ## cell water volumes
  # thickness
  thk <- {
    # takes advantage of R's matrix subsetting (arr[imtx] where imtx
    #  has as many columns as arr has dimensions)
    top <- {
      # water top in cells that pathlines pass through
      ar <- var.get.nc(wtop, "wtop", dstarts, dcounts, collapse = FALSE)
      ar[CRLts2]
    }
    bot <- {
      # layer base in cells that pathlines pass through
      ar <- var.get.nc(gwnc, "elev", dstarts[1:3] + c(0L, 0L, 1L),
                       dcounts[1:3], collapse = FALSE)
      ar[CRLts2[, -4L]]
    }

    top - bot
  }; rm(ar, top, bot)

  # horizontal dimensions
  dx <- diff(var.get.nc(gwnc, "gccs"))[CRLts[, 1L]]
  dy <- diff(var.get.nc(gwnc, "grcs"))[CRLts[, 2L]]

  # porosity
  por <- if(is.array(phi)) phi[CRLts[, 1:3, drop = FALSE]] else phi

  ## initial masses - may be single value, vector for each particle or
  ##  function of start time
  loss1 <- 0 # ensure loss1 exists
  m0p <- if(is.vector(m0)){
    if(identical(lnm0 <- length(m0), 1L)) rep(m0, max(ptl$ptlno)) else{
      if(identical(lnm0, length(unique(ptl$ptlno)))){
        m0tmp <- double(max(ptl[, ptlno], length(m0)))
        m0tmp[unique(ptl$ptlno)] <- m0
        m0tmp
      }else if(lnm0 >= max(ptl$ptlno)){
        # loss records mass not released due to particle being stranded, but
        #   cannot pick up particles that were not released because they
        #   were outside the model domain, because such particles are not
        #   recorded by MODPATH
        if(loss) loss1 <- ifelse(1:lnm0 %in% ptl$ptlno, 0, m0)
        m0
      }else stop("invalid m0: ",
                 "shorter than number of pathlines but not length 1\n")
    }
  }else if(is.function(m0)){
    vapply(vapply(ptl, `[`, double(1), 1L, "t"), m0, double(1))
  }else stop("invalid m0: ",
             "must be single value, vector or function of time\n")

  # add columns for flow imbalance and cell water volumes
  set(ptl, NULL, c("Qimb", "V"), list(imbalance, thk*dx*dy*por))

  # NA volumes signify dry or inactive cells, and any mass within them
  #  should be registered as loss
  if(any(is.na(ptl$V))){
    ptl[is.na(V),
        .SD[, loss1[ptlno] <<- loss1[ptlno] + sum(m0p[ptlno]), by = ptlno],
        .SD = c("ptlno", "V")]

    ptl <- ptl[!is.na(V)]
  }

  # total time range of simulation (used for when there is a single ptl entry for a particle)
  maxt <- if(missing(end.t)) max(ptl$t) else end.t
  totdur <- maxt - min(ptl$t)

  ## apply MassTrack calculations on each pathline
  # .N is the number of rows within the ptlno group
  ptl[, c("m", if(outmass) "ml", if(react.loss) "mrl") := if(identical(.N, 1L)){
    # case in which there is only one entry for the particle (probably instant capture)
    pVdiff <- Qimb*totdur/V

    m1 <- m0p[ptlno]*ifelse(pVdiff < 0, ifelse(pVdiff < -1, 0, 1 + pVdiff), 1)
    m1 <- ifelse(m1 < 0, 0, m1) # that's all folks (mass completely drained if m1 < 0, so reset to 0)

    # mass lost
    if(outmass){
      ml1nd <- -(m1 - m0p[ptlno])
      ml1 <- ml1nd*exp(-linear.decay*totdur/exp(1))
    }

    # linear decay - mass is removed after abstractions given to outmass
    if(linear.decay && m1){
      m1nd <- m1
      m1 <- m1nd*exp(-linear.decay*totdur)
    }else m1nd <- 0

    # mass loss due to reaction
    if(react.loss) mrl1 <- m1nd - m1 + if(outmass) ml1nd - ml1 else 0

    # return final mass and mass loss
    mget(c("m1", if(outmass) "ml1", if(react.loss) "mrl1"))
  }else{

    durs <- diff(c(t, maxt)) #amount of time spent on each step
    pVdiff <- Qimb*durs/V #proportional volume change occurring in the cells which the particle passes through

    m1 <- cumprod(ifelse(pVdiff < 0, 1 + pVdiff, 1))*m0p[ptlno] #evolving particle mass due to partial or complete abstraction
    if(any(m1 < 0)) m1[min(which(m1 < 0)):.N] <- 0 #all mass depleted when m1 < 0: set this and following to 0

    #update output mass
    if(outmass){
      ml1nd <- -diff(c(m0p[ptlno], m1))
      ml1 <- ml1nd*exp(-linear.decay*((c(t[-1L], maxt)*(exp(1) - 1) + t)/exp(1) - t[1L]))
    }

    #linear decay applied after abstracted mass given to outmass
    if(linear.decay){
      m1nd <- m1
      m1 <- m1nd*exp(-linear.decay*(t - t[1L]))
    }else m1nd <- 0

    if(react.loss) mrl1 <- diff(c(0, m1nd - m1 + if(outmass) ml1nd - ml1 else 0))

    #return mass track
    mget(c("m1", if(outmass) "ml1", if(react.loss) "mrl1"))
  }, by = ptlno]

  # apply mass loss to outmass array
  if(outmass && outmass.array){
    if(tlumpoutmass){
      ptl[, out[C, R, L] <<- sum(ml), by = list(C, R, L)]
    }else{
      ptl[, out[C, R, L, timestep] <<- sum(ml), by = list(C, R, L, timestep)]
    }
  }
  if(react.loss && outmass.array){
    if(tlumpoutmass){
      ptl[, outr[C, R, L] <<- sum(mrl), by = list(C, R, L)]
    }else{
      ptl[, outr[C, R, L, timestep] <<- sum(mrl), by = list(C, R, L, timestep)]
    }
  }

  # unneeded columns
  set(ptl, NULL, c("Qimb", "V", if(outmass && outmass.array) "ml",
                   if(react.loss && outmass.array) "mrl"), NULL)

  # package results into list
  res <- list(outmass = if(outmass && outmass.array) out,
              reactloss = if(react.loss && outmass.array) outr,
              loss = if(loss) loss1, traces = ptl)
  return(res[!sapply(res, is.null)])
}


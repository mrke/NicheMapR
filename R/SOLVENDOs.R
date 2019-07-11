#' SOLVENDOs
#'
#' R wrapper for Fortran binary of SOLVENDOs (endotherm model, multipart)
#' @param AMASS A
#' @export
SOLVENDOs <- function(SOLVENDOs.input, NPART, TAs, TGRDs, TSs, TSKYs, VELs, RHs, QSOLRs,
                     FLTYPEs, AMASSs, GMULTs, TCONDSBs, TBUSHs, SHADEs, MXWETs, AKMAXs,
                     NGEOMs, MAXPTVENs, PTCONDs, FURTHRMKs, DHAIRDs, DHAIRVs, LHAIRDs,
                     LHAIRVs, ZFURDs, ZFURVs, RHODs, RHOVs, REFLDs, REFLVs, EMISANs,
                     FATOBJs, FABUSHs, FGDREFs, FSKREFs, TCs, TCMAXs, AK1s, SKINWs, BAREVAPs,
                     PCTBAREVAPs, PCTDIFs, SUBQFATs, FATPCTs, PARTMULTs, ANDENSs){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/SOLVENDOs.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/SOLVENDOs.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/SOLVENDOs.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/SOLVENDOs.so'
  }
  if (!is.loaded('SOLVENDO')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("SOLVENDOs",
                as.double(SOLVENDOs.input),
                as.integer(NPART),
                as.double(TAs),
                as.double(TGRDs),
                as.double(TSs),
                as.double(TSKYs),
                as.double(VELs),
                as.double(RHs),
                as.double(QSOLRs),
                as.double(FLTYPEs),
                as.double(AMASSs),
                as.double(GMULTs),
                as.double(TCONDSBs),
                as.double(TBUSHs),
                as.double(SHADEs),
                as.double(MXWETs),
                as.double(AKMAXs),
                as.double(NGEOMs),
                as.double(MAXPTVENs),
                as.double(PTCONDs),
                as.double(FURTHRMKs),
                as.double(DHAIRDs),
                as.double(DHAIRVs),
                as.double(LHAIRDs),
                as.double(LHAIRVs),
                as.double(ZFURDs),
                as.double(ZFURVs),
                as.double(RHODs),
                as.double(RHOVs),
                as.double(REFLDs),
                as.double(REFLVs),
                as.double(EMISANs),
                as.double(FATOBJs),
                as.double(FABUSHs),
                as.double(FGDREFs),
                as.double(FSKREFs),
                as.double(TCs),
                as.double(TCMAXs),
                as.double(AK1s),
                as.double(SKINWs),
                as.double(BAREVAPs),
                as.double(PCTBAREVAPs),
                as.double(PCTDIFs),
                as.double(SUBQFATs),
                as.double(FATPCTs),
                as.double(PARTMULTs),
                as.double(ANDENSs),
                results=matrix(data = 0., nrow = 1, ncol = 30),
                SIMULSOLouts=matrix(data=0., nrow = 2*NPART, ncol = 17),
                PACKAGE = "SOLVENDOs")
  #dyn.unload("SOLVENDOs.dll")

  whole <- matrix(data = 0, nrow = 1, ncol = 30)
  parts <- matrix(data = 0, nrow = 2 * NPART, ncol = 17)

  storage.mode(whole)<-"double"
  storage.mode(parts)<-"double"
  parts <- a$results
  parts <- a$SIMULSOLouts

  return(list(whole, parts))
}

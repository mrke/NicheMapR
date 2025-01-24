#' plot_human - auxiliary function for HomoTherm
#'
#' This function plots the human body and insulation.
#' @encoding UTF-8
#' @param MASS = 70
#' @param HEIGHT = 170
#' @param MASSFRACs = c(0.07609801, 0.50069348, 0.04932963, 0.16227462)
#' @param DENSITYs = rep(1050, 4)
#' @param INSDEPDs = c(1e-02, rep(6e-03, 3))
#' @param INSDEPVs = c(1e-09, rep(6e-03, 3))
#' @param SUBQFATs = rep(1, 4)
#' @param FATPCTs = c(5, 36, 10, 23)
#' @param SHAPE_Bs = c(1.6, 1.9, 11, 7.0)
#' @param DHAIRDs = c(7.5e-5, rep(1E-06, 3))
#' @param DHAIRVs = c(7.5e-5, rep(1E-06, 3))
#' @param LHAIRDs = c(50e-3, 50e-3, 50e-3, 50e-3)
#' @param LHAIRVs = c(1e-9, 50e-3, 50e-3, 50e-3)
#' @param INSDENDs = rep(3e+08, 4)
#' @param INSDENVs = c(3e+05, rep(3e+08, 3))
#' @param PJOINs = c(0.02591961, 0.07879298, 0.01850589, 0.03333333)
#' @usage plot_person(MASS = 70, HEIGHT = 170, INSDEPDs = c(1e-02, rep(6e-03, 3)), INSDEPVs = c(1e-09, rep(6e-03, 3)), FATPCTs = c(5, 36, 10, 23), SHAPE_Bs = c(1.6, 1.9, 11, 7.0))
#' @export
plot_human <- function(
    MASS = 70,
    HEIGHT = 170,
    MASSFRACs = c(0.07609801, 0.50069348, 0.04932963, 0.16227462),
    SHAPEs = c(4, 1, 1, 1),
    DENSITYs = rep(1050, 4),
    INSDEPDs = c(1e-02, rep(6e-03, 3)),
    INSDEPVs = c(1e-09, rep(6e-03, 3)),
    SUBQFATs = rep(1, 4),
    FATPCTs = c(5, 36, 10, 23),
    SHAPE_Bs = c(1.6, 1.9, 11, 7.0),
    DHAIRDs = c(7.5e-5, rep(1E-06, 3)),
    DHAIRVs = c(7.5e-5, rep(1E-06, 3)),
    LHAIRDs = c(50e-3, 50e-3, 50e-3, 50e-3),
    LHAIRVs = c(1e-9, 50e-3, 50e-3, 50e-3),
    INSDENDs = rep(3e+08, 4),
    INSDENVs = c(3e+05, rep(3e+08, 3)),
    PJOINs = c(0.02591961, 0.07879298, 0.01850589, 0.03333333)
    ){

  if (!require("plotrix", quietly = TRUE)) {
    stop("package 'plotrix' is needed. Please install it.",
         call. = FALSE)
  }

  SHAPEs <- c(4, 1, 1, 1)
  ORIENTs <- rep(0, 4)
  MASSs <- MASS * MASSFRACs

  SAMODE <- 0
  GEOM.head <- GEOM_ENDO(MASSs[1], DENSITYs[1], FATPCTs[1], SHAPEs[1], INSDEPDs[1], SUBQFATs[1], SHAPE_Bs[1], SHAPE_Bs[1], DHAIRDs[1], INSDENDs[1], PJOINs[1], SAMODE, ORIENTs[1], ZEN = 90)
  GEOM.trunk <- GEOM_ENDO(MASSs[2], DENSITYs[2], FATPCTs[2], SHAPEs[2], INSDEPDs[2], SUBQFATs[2], SHAPE_Bs[2], SHAPE_Bs[2], DHAIRDs[2], INSDENDs[2], PJOINs[2], SAMODE, ORIENTs[2], ZEN = 90)
  GEOM.arm <- GEOM_ENDO(MASSs[3], DENSITYs[3], FATPCTs[3], SHAPEs[3], INSDEPDs[3], SUBQFATs[3], SHAPE_Bs[3], SHAPE_Bs[3], DHAIRDs[3], INSDENDs[3], PJOINs[3], SAMODE, ORIENTs[3], ZEN = 90)
  GEOM.leg <- GEOM_ENDO(MASSs[4], DENSITYs[4], FATPCTs[4], SHAPEs[4], INSDEPDs[4], SUBQFATs[4], SHAPE_Bs[4], SHAPE_Bs[4], DHAIRDs[4], INSDENDs[4], PJOINs[4], SAMODE, ORIENTs[4], ZEN = 90)

  GEOM.heads <- as.data.frame(GEOM.head)
  GEOM.trunks <- as.data.frame(GEOM.trunk)
  GEOM.arms <- as.data.frame(GEOM.arm)
  GEOM.legs <- as.data.frame(GEOM.leg)
  GEOM.lab <- c("VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH",
                "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP", "GMASS", "AREASKIN",
                "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK",
                "CONVAR", "R1", "R2")
  colnames(GEOM.heads) <- GEOM.lab
  colnames(GEOM.trunks) <- GEOM.lab
  colnames(GEOM.arms) <- GEOM.lab
  colnames(GEOM.legs) <- GEOM.lab
  AREA <- max(GEOM.heads$ATOT + GEOM.trunks$ATOT + GEOM.arms$ATOT * 2 + GEOM.legs$ATOT * 2)

  head.ASEMAJ <- GEOM.heads$ASEMAJ
  head.BSEMIN <- GEOM.heads$BSEMIN
  head.CSEMIN <- GEOM.heads$CSEMIN
  head.INSDEPD <- INSDEPDs[1]
  head.INSDEPV <- INSDEPVs[1]
  head.FATTHK <- GEOM.heads$FATTHK

  trunk.ALENGTH <- GEOM.trunks$ALENTH
  trunk.AWIDTH <- GEOM.trunks$AWIDTH
  trunk.AHEIGHT <- GEOM.trunks$AHEIT
  trunk.INSDEPD <- INSDEPDs[2]
  trunk.INSDEPV <- INSDEPVs[2]
  trunk.FATTHK <- GEOM.trunks$FATTHK

  arm.ALENGTH <- GEOM.arms$ALENTH
  arm.AWIDTH <- GEOM.arms$AWIDTH
  arm.AHEIGHT <- GEOM.arms$AHEIT
  arm.INSDEPD <- INSDEPDs[3]
  arm.INSDEPV <- INSDEPVs[3]
  arm.FATTHK <- GEOM.arms$FATTHK

  leg.ALENGTH <- GEOM.legs$ALENTH
  leg.AWIDTH <- GEOM.legs$AWIDTH
  leg.AHEIGHT <- GEOM.legs$AHEIT
  leg.INSDEPD <- INSDEPDs[4]
  leg.INSDEPV <- INSDEPVs[4]
  leg.FATTHK <- GEOM.legs$FATTHK

  # coronal

  #plot.height <- (head.ASEMAJ * 2 + trunk.ALENGTH + leg.ALENGTH) * 1.05
  plot.height <- 2
  plot.width <- c(-1.5, 1.5)
  trunk.left <- 0 - trunk.AWIDTH / 2

  plot(plot.width, c(0, plot.height), type="n", main="saggital section / coronal section / transverse section", ylab = 'metres', xlab = 'metres', asp=1, xaxs = 'i', yaxs = 'i')
  #abline(v = 0, lty = 2)

  # trunk fur
  rect(xleft = trunk.left - INSDEPDs[2],
       ybottom = leg.ALENGTH - INSDEPDs[2],
       xright = trunk.left + trunk.AWIDTH + INSDEPDs[2],
       ytop = leg.ALENGTH + trunk.ALENGTH + INSDEPDs[2], col = 'lightblue')

  # head fur
  top.of.head <- head.ASEMAJ + leg.ALENGTH + trunk.ALENGTH
  draw.ellipse(c(0, 0), c(top.of.head, top.of.head),
               b = head.ASEMAJ + INSDEPDs[1],
               a = head.BSEMIN + INSDEPDs[1], col = "brown",
               border = NA, segment = c(0, 180), deg = TRUE)
  #head naked
  draw.ellipse(c(0, 0),
               c(top.of.head,
                 top.of.head),
               col = "grey", b = head.ASEMAJ, a = head.BSEMIN, lty = 2)
  draw.ellipse(c(0, 0),
               c(top.of.head,
                 top.of.head),
               col = "pink", b = head.ASEMAJ - head.FATTHK, a = head.BSEMIN - head.FATTHK, lty = 0)

  # left leg
  xleft.leg.fur <- trunk.left - trunk.AWIDTH / 10 - INSDEPDs[4]
  xleft.leg.naked <- trunk.left - trunk.AWIDTH / 10
  xright.leg.fur <- trunk.left - trunk.AWIDTH / 10 + leg.AWIDTH + INSDEPDs[4]
  xright.leg.naked <- trunk.left - trunk.AWIDTH / 10 + leg.AWIDTH
  ytop.leg.fur <- leg.ALENGTH + INSDEPDs[4]
  ytop.leg.naked <- leg.ALENGTH

  rect(xleft = xleft.leg.fur,
       ybottom = 0 - INSDEPDs[4],
       xright = xright.leg.fur,
       ytop = ytop.leg.fur, col = 'khaki')
  rect(xleft = xleft.leg.naked + INSDEPDs[4],
       ybottom = 0,
       xright = xright.leg.naked - INSDEPDs[4],
       ytop = ytop.leg.naked, col = 'lightgrey', lty = 2)
  rect(xleft = xleft.leg.naked + INSDEPDs[4] + leg.FATTHK,
       ybottom = 0 + leg.FATTHK,
       xright = xright.leg.naked - INSDEPDs[4] - leg.FATTHK,
       ytop = ytop.leg.naked - leg.FATTHK, col = 'pink', lty = 0)

  # right leg
  rect(xleft = xright.leg.fur * -1,
       ybottom = 0 - INSDEPDs[4],
       xright = xleft.leg.fur * -1,
       ytop = ytop.leg.fur, col = 'khaki')
  rect(xleft = xright.leg.naked * -1 + INSDEPDs[4],
       ybottom = 0,
       xright = xleft.leg.naked * -1 - INSDEPDs[4],
       ytop = ytop.leg.naked, col = 'lightgrey', lty = 2)
  rect(xleft = xright.leg.naked * -1 + INSDEPDs[4] + leg.FATTHK,
       ybottom = 0 + leg.FATTHK,
       xright = xleft.leg.naked * -1 - INSDEPDs[4] - leg.FATTHK,
       ytop = ytop.leg.naked - leg.FATTHK, col = 'pink', lty = 0)

  # trunk naked
  rect(xleft = trunk.left + INSDEPDs[2],
       ybottom = leg.ALENGTH,
       xright = trunk.left + trunk.AWIDTH - INSDEPDs[2],
       ytop = leg.ALENGTH + trunk.ALENGTH,
       col = 'lightgrey', lty = 2)
  rect(xleft = trunk.left + trunk.FATTHK + INSDEPDs[2],
       ybottom = leg.ALENGTH + trunk.FATTHK,
       xright = trunk.left + trunk.AWIDTH - INSDEPDs[2] - trunk.FATTHK,
       ytop = leg.ALENGTH + trunk.ALENGTH - trunk.FATTHK,
       col = 'pink', lty = 0)

  # left arm
  xleft.arm.fur <- trunk.left - arm.AWIDTH - arm.AWIDTH / 4 - INSDEPDs[3]
  xleft.arm.naked <- trunk.left - arm.AWIDTH - arm.AWIDTH / 4
  xright.arm.fur <- trunk.left - arm.AWIDTH / 4 - INSDEPDs[2]
  xright.arm.naked <- trunk.left - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3]
  ytop.arm.fur <- leg.ALENGTH + trunk.ALENGTH + INSDEPDs[3]
  ytop.arm.naked <- leg.ALENGTH + trunk.ALENGTH
  ybottom.arm.fur <- head.ASEMAJ + leg.ALENGTH + INSDEPDs[4] + trunk.ALENGTH - arm.ALENGTH - INSDEPDs[3]
  ybottom.arm.naked <- head.ASEMAJ + leg.ALENGTH + INSDEPDs[4] + trunk.ALENGTH - arm.ALENGTH
  rect(xleft = xleft.arm.fur,
       ybottom = ybottom.arm.fur,
       xright = xright.arm.fur,
       ytop = ytop.arm.fur, col = 'lightblue')
  rect(xleft = xleft.arm.naked,
       ybottom = ybottom.arm.naked,
       xright = xright.arm.naked,
       ytop = ytop.arm.naked, col = 'grey', lty = 2)
  rect(xleft = xleft.arm.naked + arm.FATTHK,
       ybottom = ybottom.arm.naked + arm.FATTHK,
       xright = xright.arm.naked - arm.FATTHK,
       ytop = ytop.arm.naked - arm.FATTHK, col = 'pink', lty = 0)

  # right arm
  rect(xleft = xright.arm.fur * -1,
       ybottom = ybottom.arm.fur,
       xright = xleft.arm.fur * -1,
       ytop = ytop.arm.fur, col = 'lightblue')
  rect(xleft = xright.arm.naked * -1,
       ybottom = ybottom.arm.naked,
       xright = xleft.arm.naked * -1,
       ytop = ytop.arm.naked, col = 'grey', lty = 2)
  rect(xleft = xright.arm.naked * -1 + arm.FATTHK,
       ybottom = ybottom.arm.naked + arm.FATTHK,
       xright = xleft.arm.naked * -1 - arm.FATTHK,
       ytop = ytop.arm.naked - arm.FATTHK, col = 'pink', lty = 0)

  # sagittal

  # trunk fur
  rect(xleft = trunk.left - INSDEPDs[2] - 1,
       ybottom = leg.ALENGTH - INSDEPDs[2] + INSDEPDs[4],
       xright = trunk.left + trunk.AWIDTH + INSDEPDs[2] - 1,
       ytop = leg.ALENGTH + trunk.ALENGTH + INSDEPDs[2] + INSDEPDs[4],
       col = 'lightblue')

  #head fur
  draw.ellipse(c(0, 0) - 1, c(top.of.head, top.of.head),
               b = head.ASEMAJ + INSDEPDs[1],
               a = head.CSEMIN + INSDEPDs[1], col = "brown",
               border = NA, segment = c(0, 180), deg = TRUE)
  #head naked
  draw.ellipse(c(0, 0) - 1,
               c(top.of.head,
                 top.of.head),
               col = "grey", b = head.ASEMAJ, a = head.CSEMIN, lty = 2)
  draw.ellipse(c(0, 0) - 1,
               c(top.of.head,
                 top.of.head),
               col = "pink", b = head.ASEMAJ - head.FATTHK, a = head.CSEMIN - head.FATTHK, lty = 0)


  # left leg
  xleft.leg.fur <- -leg.AWIDTH / 2 - INSDEPDs[4] - 1
  xleft.leg.naked <- -leg.AWIDTH / 2 - 1
  xright.leg.fur <- -leg.AWIDTH / 2 - 1 + leg.AWIDTH + INSDEPDs[4]
  xright.leg.naked <- -leg.AWIDTH / 2 - 1 + leg.AWIDTH

  # left leg
  rect(xleft = xleft.leg.fur,
       ybottom = 0,
       xright = xright.leg.fur,
       ytop = ytop.leg.fur, col = 'khaki')
  rect(xleft = xleft.leg.naked,
       ybottom = 0 + INSDEPDs[4],
       xright = xright.leg.naked,
       ytop = ytop.leg.naked, col = 'grey', lty = 2)
  rect(xleft = xleft.leg.naked + leg.FATTHK,
       ybottom = 0 + INSDEPDs[4] + leg.FATTHK,
       xright = xright.leg.naked - leg.FATTHK,
       ytop = ytop.leg.naked - leg.FATTHK, col = 'pink', lty = 0)

  # trunk naked
  rect(xleft = trunk.left - 1,
       ybottom = leg.ALENGTH + INSDEPDs[4],
       xright = trunk.left + trunk.AWIDTH - 1,
       ytop = leg.ALENGTH + trunk.ALENGTH + INSDEPDs[4], col = 'grey', lty = 2)
  rect(xleft = trunk.left - 1 + trunk.FATTHK,
       ybottom = leg.ALENGTH + INSDEPDs[4] + trunk.FATTHK,
       xright = trunk.left + trunk.AWIDTH - 1 - trunk.FATTHK,
       ytop = leg.ALENGTH + trunk.ALENGTH + INSDEPDs[4] - trunk.FATTHK, col = 'pink', lty = 0)

  # left arm
  xleft.arm.fur <- -arm.AWIDTH / 2 - INSDEPDs[3] - 1
  xleft.arm.naked <- -arm.AWIDTH / 2 - 1
  xright.arm.fur <- -arm.AWIDTH / 2 - 1 + arm.AWIDTH + INSDEPDs[3]
  xright.arm.naked <- -arm.AWIDTH / 2 - 1 + arm.AWIDTH
  rect(xleft = xleft.arm.fur,
       ybottom = ybottom.arm.fur,
       xright = xright.arm.fur,
       ytop = ytop.arm.fur, col = 'lightblue', lty = 2)
  rect(xleft = xleft.arm.naked,
       ybottom = ybottom.arm.naked,
       xright = xright.arm.naked,
       ytop = ytop.arm.naked, col = 'grey', lty = 2)
  rect(xleft = xleft.arm.naked + arm.FATTHK,
       ybottom = ybottom.arm.naked + arm.FATTHK,
       xright = xright.arm.naked - arm.FATTHK,
       ytop = ytop.arm.naked - arm.FATTHK, col = 'pink', lty = 0)


  # bird's eye

  mid.height <- (leg.ALENGTH + trunk.ALENGTH + head.ASEMAJ * 2) / 2
  top.height <- mid.height * 1.5
  bottom.height <- mid.height / 2
  abline(h = mid.height, lty = 2, col = 'grey')
  #abline(h = top.height, lty = 2, col = 'grey')
  abline(h = bottom.height, lty = 2, col = 'grey')

  # legs
  draw.ellipse(c(-trunk.AWIDTH / 20 - leg.AWIDTH / 2 - INSDEPDs[4],
                 trunk.AWIDTH / 20 + leg.AWIDTH / 2 + INSDEPDs[4]) + 1,
               c(0, 0) + top.height, col = "khaki",
               b = leg.AWIDTH / 2 + INSDEPDs[4],
               a = leg.AWIDTH / 2 + INSDEPDs[4])
  if(INSDEPDs[4] < 0.0001){
    draw.ellipse(c(-trunk.AWIDTH / 20 - leg.AWIDTH / 2 - INSDEPDs[4],
                   trunk.AWIDTH / 20 + leg.AWIDTH / 2 + INSDEPDs[4]) + 1,
                 c(0, 0) + top.height, col = "pink",
                 b = leg.AWIDTH / 2,
                 a = leg.AWIDTH / 2, lty = 2)
  }
  # trunk
  draw.ellipse(c(0, 0) + 1,
               c(0, 0) + top.height, col = "lightblue",
               b = trunk.AWIDTH / 2 + INSDEPDs[2],
               a = trunk.AWIDTH / 2 + INSDEPDs[2])
  if(INSDEPDs[2] < 0.0001){
    draw.ellipse(c(0, 0) + 1,
                 c(0, 0) + top.height, col = "pink",
                 b = trunk.AWIDTH / 2,
                 a = trunk.AWIDTH / 2, lty = 2)
  }
  # head
  if(INSDEPDs[2] < 0.0001){
    draw.ellipse(c(0, 0) + 1,
                 c(0, 0) + top.height,
                 col = "pink",
                 b = head.ASEMAJ + INSDEPDs[1],
                 a = head.BSEMIN + INSDEPDs[1])
  }else{
    draw.ellipse(c(0, 0) + 1,
                 c(0, 0) + top.height,
                 col = "brown",
                 b = head.ASEMAJ + INSDEPDs[1],
                 a = head.BSEMIN + INSDEPDs[1])
  }
  # arms
  draw.ellipse(c(-trunk.AWIDTH / 2 - arm.AWIDTH / 2 - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3],
                 trunk.AWIDTH / 2 + arm.AWIDTH / 2 + arm.AWIDTH / 4 + INSDEPDs[2] + INSDEPDs[3]) + 1,
               c(0, 0) + top.height, col = "lightblue",
               b = arm.AWIDTH / 2 + INSDEPDs[3],
               a = arm.AWIDTH / 2 + INSDEPDs[3])
  if(INSDEPDs[4] < 0.0001){
    draw.ellipse(c(-trunk.AWIDTH / 2 - arm.AWIDTH / 2 - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3],
                   trunk.AWIDTH / 2 + arm.AWIDTH / 2 + arm.AWIDTH / 4 + INSDEPDs[2] + INSDEPDs[3]) + 1,
                 c(0, 0) + top.height, col = "pink",
                 b = arm.AWIDTH / 2,
                 a = arm.AWIDTH / 2, lty = 2)
  }
  # transverse 1

  # trunk
  draw.ellipse(c(0, 0) + 1,
               c(0, 0) + mid.height, col = "lightblue",
               b = trunk.AWIDTH / 2 + INSDEPDs[2],
               a = trunk.AWIDTH / 2 + INSDEPDs[2])
  draw.ellipse(c(0, 0) + 1,
               c(0, 0) + mid.height, col = "grey",
               b = trunk.AWIDTH / 2,
               a = trunk.AWIDTH / 2, lty = 2)
  draw.ellipse(c(0, 0) + 1,
               c(0, 0) + mid.height, col = "pink",
               b = trunk.AWIDTH / 2 - trunk.FATTHK,
               a = trunk.AWIDTH / 2 - trunk.FATTHK, lty = 0)
  # arms
  draw.ellipse(c(-trunk.AWIDTH / 2 - arm.AWIDTH / 2 - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3],
                 trunk.AWIDTH / 2 + arm.AWIDTH / 2 + arm.AWIDTH / 4 + INSDEPDs[2] + INSDEPDs[3]) + 1,
               c(0, 0) + mid.height, col = "lightblue",
               b = arm.AWIDTH / 2 + INSDEPDs[3],
               a = arm.AWIDTH / 2 + INSDEPDs[3])
  draw.ellipse(c(-trunk.AWIDTH / 2 - arm.AWIDTH / 2 - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3],
                 trunk.AWIDTH / 2 + arm.AWIDTH / 2 + arm.AWIDTH / 4 + INSDEPDs[2] + INSDEPDs[3]) + 1,
               c(0, 0) + mid.height, col = "grey",
               b = arm.AWIDTH / 2,
               a = arm.AWIDTH / 2, lty = 2)
  draw.ellipse(c(-trunk.AWIDTH / 2 - arm.AWIDTH / 2 - arm.AWIDTH / 4 - INSDEPDs[2] - INSDEPDs[3],
                 trunk.AWIDTH / 2 + arm.AWIDTH / 2 + arm.AWIDTH / 4 + INSDEPDs[2] + INSDEPDs[3]) + 1,
               c(0, 0) + mid.height, col = "pink",
               b = arm.AWIDTH / 2 - arm.FATTHK,
               a = arm.AWIDTH / 2 - arm.FATTHK, lty = 2)


  # transverse 2

  # legs
  draw.ellipse(c(-trunk.AWIDTH / 20 - leg.AWIDTH / 2 - INSDEPDs[4],
                 trunk.AWIDTH / 20 + leg.AWIDTH / 2 + INSDEPDs[4]) + 1,
               c(0, 0) + bottom.height, col = "khaki",
               b = leg.AWIDTH / 2 + INSDEPDs[4],
               a = leg.AWIDTH / 2 + INSDEPDs[4])
  draw.ellipse(c(-trunk.AWIDTH / 20 - leg.AWIDTH / 2 - INSDEPDs[4],
                 trunk.AWIDTH / 20 + leg.AWIDTH / 2 + INSDEPDs[4]) + 1,
               c(0, 0) + bottom.height, col = "grey",
               b = leg.AWIDTH / 2,
               a = leg.AWIDTH / 2, lty = 2)
  draw.ellipse(c(-trunk.AWIDTH / 20 - leg.AWIDTH / 2 - INSDEPDs[4],
                 trunk.AWIDTH / 20 + leg.AWIDTH / 2 + INSDEPDs[4]) + 1,
               c(0, 0) + bottom.height, col = "pink",
               b = leg.AWIDTH / 2 - leg.FATTHK,
               a = leg.AWIDTH / 2 - leg.FATTHK, lty = 0)
  height <- GEOM.heads$ALENTH + GEOM.trunks$ALENTH + GEOM.legs$ALENTH # m
  abline(h = height + INSDEPDs[4], lty = 2)
  abline(h = HEIGHT / 100 + INSDEPDs[4], col = 'red', lty = 2)
  return(height)
}

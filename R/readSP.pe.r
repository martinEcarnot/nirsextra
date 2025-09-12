#' readSP.pe
#'
#' Read one PerkinElmer file
#'
#' @export
#'

readSP.pe <- function(filename){

#% Reads in spectra from PerkinElmer block structured files.
# This version supports 'Spectrum' SP files.
# Note that earlier 'Data Manager' formats are not supported.
#
#   This function gives as oputpout a List with xData, yData
#   and other miscellanous information in name,value pairs

# Copyright for Matlab Version spload.m (C)2007 PerkinElmer Life and Analytical Sciences
# Stephen Westlake, Seer Green
# spload.m was transformed to readSP.pe.R by Karl Molt, Universit?t Duisburg-Essen, 2009
#
# History
# 2007-04-24 SW     Initial Matlab version
# 2009-10-04 KM     Initial R version
# Block IDs

DSet2DC1DIBlock               =  120
HistoryRecordBlock            =  121
InstrHdrHistoryRecordBlock    =  122
InstrumentHeaderBlock         =  123
IRInstrumentHeaderBlock       =  124
UVInstrumentHeaderBlock       =  125
FLInstrumentHeaderBlock       =  126
# Data member IDs
DataSetDataTypeMember              =  -29839
DataSetAbscissaRangeMember         =  -29838
DataSetOrdinateRangeMember         =  -29837
DataSetIntervalMember              =  -29836
DataSetNumPointsMember             =  -29835
DataSetSamplingMethodMember        =  -29834
DataSetXAxisLabelMember            =  -29833
DataSetYAxisLabelMember            =  -29832
DataSetXAxisUnitTypeMember         =  -29831
DataSetYAxisUnitTypeMember         =  -29830
DataSetFileTypeMember              =  -29829
DataSetDataMember                  =  -29828
DataSetNameMember                  =  -29827
DataSetChecksumMember              =  -29826
DataSetHistoryRecordMember         =  -29825
DataSetInvalidRegionMember         =  -29824
DataSetAliasMember                 =  -29823
DataSetVXIRAccyHdrMember           =  -29822
DataSetVXIRQualHdrMember           =  -29821
DataSetEventMarkersMember          =  -29820
# Type code IDs
ShortType               = 29999
UShortType              = 29998
IntType                 = 29997
UIntType                = 29996
LongType                = 29995
BoolType                = 29988
CharType                = 29987
CvCoOrdPointType        = 29986
StdFontType             = 29985
CvCoOrdDimensionType    = 29984
CvCoOrdRectangleType    = 29983
RGBColorType            = 29982
CvCoOrdRangeType        = 29981
DoubleType              = 29980
CvCoOrdType             = 29979
ULongType               = 29978
PeakType                = 29977
CoOrdType               = 29976
RangeType               = 29975
CvCoOrdArrayType        = 29974
EnumType                = 29973
LogFontType             = 29972

nbyt=0
fid <- file(filename, "rb")
Signature <- readChar(fid,4)

if (Signature !=  "PEPE"){
stop("This is not a PerkinElmer block structured file!")
}
description <- readChar(fid,40)
xLen <- 0

# The rest of the file is a list of blocks
spdata.gelesen = FALSE
Ausgabe.Liste <- NULL

# while (T) {
  # while (!isEof.connection(fid) & !spdata.gelesen){
  while (!spdata.gelesen){

    blockID <- readBin(fid,"integer",1,size=2)
    blockSize <- readBin(fid,"integer",1,size=4)


# if (blockID==-29829) {browser()}
# if (isEof.connection(fid)){
# stop("Reading of spectral data was not possible! End of file was reached!")
# break
# }


if (blockID==DSet2DC1DIBlock){
        # Wrapper block.  Read nothing.
        }
        else{     #else 1
        if (blockID==DataSetAbscissaRangeMember){
            innerCode = readBin(fid,"integer",1,size=2)
            #%_ASSERTE(CvCoOrdRangeType == nInnerCode)
            x0 = readBin(fid, "double",1)
            xEnd = readBin(fid, "double",1)
            Ausgabe.Liste$x0 <- x0
            Ausgabe.Liste$xEnd <- xEnd
            }
            else{  #else 2

        if (blockID== DataSetIntervalMember){
            innerCode = readBin(fid,"integer",1,size=2)
            xDelta = readBin(fid, "double",1)
            Ausgabe.Liste$xDelta <- xDelta
            xData <- seq(Ausgabe.Liste$x0,Ausgabe.Liste$xEnd,Ausgabe.Liste$xDelta)
            Ausgabe.Liste$xData <- xData
            }
             else{        # else 3
        if (blockID== DataSetNumPointsMember){
            innerCode = readBin(fid,"integer",1,size=2)
            xLen = readBin(fid, "integer", 1,size=4)
            Ausgabe.Liste$xLen <- xLen
            if (xLen != length(xData)){
            stop("Length of xData is not compatible with xLen!")
            }
            }
            else{     # else 4

        if (blockID== DataSetXAxisLabelMember){
            innerCode = readBin(fid,"integer",1,size=2)
            len = readBin(fid,"integer", 1,size=2)
            xLabel = readChar(fid, len)
            Ausgabe.Liste$xLabel <- xLabel
            }

            else{ #else 5
         if (blockID == DataSetYAxisLabelMember){
            innerCode = readBin(fid,"integer",1,size=2)
            len = readBin(fid,"integer", 1,size=2)
            yLabel = readChar(fid, len)
            Ausgabe.Liste$yLabel <- yLabel
            }

         else{ # else 6

        if (blockID == DataSetAliasMember){
            innerCode = readBin(fid,"integer",1,size=2)
            len = readBin(fid,"integer", 1,size=2)
            Alias = readChar(fid, len)
            Ausgabe.Liste$Alias <- Alias
            }
            else{   # else 7

        if (blockID == DataSetNameMember){
            innerCode = readBin(fid,"integer",1,size=2)
            len = readBin(fid,"integer", 1,size=2)
            originalName = readChar(fid, len)
            Ausgabe.Liste$originalName <- originalName
            }
            else{ # else 8

        if (blockID == DataSetDataMember){
            innerCode = readBin(fid,"integer",1,size=2)
            len = readBin(fid,"integer", 1,size=4);
            # innerCode should be CvCoOrdArrayType
            # len should be xLen * 8
            if (xLen == 0){xLen = len / 8}
            spdata = readBin(fid, "double",xLen)
            # R Platform uses endian="little"
            Ausgabe.Liste$yData <- spdata
            spdata.gelesen = TRUE
            }
            else{  #else 9
                       # unknown block, just seek past it
           oldpos <- seek(fid,blockSize,"current")

}}}}}}}}}   # Close all 9 else

} #end while
close(fid)
#
# if (xLen == 0){
#     stop("The file does not contain spectral data")
#     }
#
    return(Ausgabe.Liste)
#     close(fid)
#
} # end function


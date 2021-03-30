read_brimrose=function(f) {

# f="C:\\Users\\seedmeister\\Desktop\\brimrose coffee 04012016\\spectres & resultats\\transmission\\sedmeister coffee 1nm trans_gain4.dat"

f='D:\\Brimrose\\test.dat'
fi <- file(f, "rb")
# hdr=readChar(fi,2)

BFF_FILE_LENGTH=260
BFF_USERID_LENGTH=30
BFF_DESCR_LENGTH=50
MAX_PROCS=20
INGS_MAX=20
INGS_LEN=15
SENS_MAX=4
SENS_LEN=15

seek(fi,0)
id=readBin(fi, "int",1,4)
userid=readChar(fi,BFF_USERID_LENGTH,useBytes=TRUE)
descr=readChar(fi,BFF_DESCR_LENGTH,useBytes=TRUE)
subdes=readChar(fi,BFF_DESCR_LENGTH,useBytes=TRUE)

seek(fi)
subcnt=readBin(fi, integer(),1,4)
start=readBin(fi, integer(),1,2)
stop=readBin(fi, integer(),1,2)
increment=readBin(fi, integer(),1,2)
numpoints=readBin(fi, integer(),1,2)
numscans=readBin(fi, integer(),1,2)
entries=readBin(fi, integer(),1,4)
gain=readBin(fi, integer(),1,2)
channel=readBin(fi, integer(),1,2)
scantype=readBin(fi, integer(),1,2)
trigger=readBin(fi, integer(),1,2)
background=readChar(fi,BFF_FILE_LENGTH,useBytes=TRUE)
node=readBin(fi, integer(),1,1)
preproc=readBin(fi, integer(),1,1)
prelist=readChar(fi,1,useBytes=TRUE)
prelist_const=readBin(fi, "numeric",MAX_PROCS,4)
ings=readBin(fi, integer(),1,1)
ingname=readChar(fi,INGS_MAX*INGS_LEN,useBytes=TRUE)
sensors=readBin(fi, integer(),1,1)
sensorname=readBin(fi, integer(),1,2)
scanid=readBin(fi, "int",1,4)                 ## to calculate date and time (number of seconds since 1.1.1970)
descr=readChar(fi,BFF_DESCR_LENGTH,useBytes=TRUE) ## subhead's description
ing_vals=readBin(fi, "numeric",INGS_MAX,4)[] ## ingredient values (one for each defined in MAINHEAD)
sensor_vals=readBin(fi, "numeric",SENS_MAX,4)[] ##sensor data (one for each defined in MAINHEAD)
reserved=readBin(fi, integer(),1,2)


close(fi)
#
# data type     # of bits    # of bytes
# char              8            1
# unsigned char      8           1
# short             16           2
# unsigned short    16           2
# long              32           4
# unsigned long     32           4
# float             32           4
#
#
# typedef struct {
#   long id;                             // FILE ID
#   char userid[BFF_USERID_LENGTH];
#   char descr[BFF_DESCR_LENGTH];        // file decription
#   char subdes[BFF_DESCR_LENGTH];       // subhead's description ( used for appending )
#   unsigned long  subcnt;               // current number in sub group (used for appending )
#   unsigned short start;                // start wavelength [nm]
#   unsigned short stop;                 // end wavelength [nm]
#   unsigned short increment;            // wavelength increment [nm]
#
#   unsigned short numpoints;            // (stop - start)/increment + 1
#   unsigned short numscans;             // number of spectra taken per average
#   unsigned long  entries;              // total number of spectra (subheaders + data)
#   unsigned short gain;                 // spectrometer's gain (1,2,4,8)
#
#   unsigned short channel;              // SAMPLESCAN - sample channel scan,
#   // REFSCAN - internal reference channel scan,
#   // POLYSCAN - internal standard channel scan
#
#   unsigned short scantype;             // RATIO - ratio scan
#   // UNITS
#   // BACKCORR
#   unsigned short trigger;              // = 0 - no trigger (othewise trigger = delay)
#
#   char background[BFF_FILE_LENGTH];    // background file name
#
#   unsigned char node;                  // spectrometer number (1-255)
#
#   unsigned char preproc;               // number of times files was processed
#   unsigned short prelist[MAX_PROCS];   // list of processing (see #defines above)
#                                                               float prelist_const[MAX_PROCS];      // to keep constants for the processing
#
#                                                               unsigned char ings;                  // number of ingredients (constituents) defined
#                                                               char ingname[INGS_MAX][INGS_LEN];    // ingredient names
#
#                                                               unsigned char sensors;               // number of external sensors used (max = 4)
#                                                               char sensorname[SENS_MAX][SENS_LEN]; // sensor names
#
#                                                               unsigned short reserved;             // unused
# } MAINHEAD_4;
# #pragma pack()
#
#
# #pragma pack(1)
# typedef struct {
#   long scanid;                         // to calculate date and time (number of seconds since 1.1.1970)
#   char descr[BFF_DESCR_LENGTH];        // subhead's description
#   float ing_vals[INGS_MAX];            // ingredient values (one for each defined in MAINHEAD)
#   float sensor_vals[SENS_MAX];         // sensor data (one for each defined in MAINHEAD)
#   unsigned short reserved;             // unused
# } SUBHEAD;
#   #pragma pack()
#


return(sp)

}

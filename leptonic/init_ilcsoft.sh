export ILCSOFT=/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16


#--------------------------------------------------------------------------------
#     MySQL
#--------------------------------------------------------------------------------
export MYSQL_HOME="/afs/cern.ch/sw/lcg/external/mysql/5.1.45/x86_64-slc5-gcc43-opt"
export MYSQL_LIBDIR="$MYSQL_HOME/lib64/mysql"
export MYSQL_PATH="$MYSQL_HOME"
export MYSQL="$MYSQL_HOME"
export PATH="$MYSQL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MYSQL_HOME/lib64/mysql:$MYSQL_HOME/lib64:$MYSQL_HOME/lib/mysql:$MYSQL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Java
#--------------------------------------------------------------------------------
export JAVA_HOME="/afs/cern.ch/sw/lcg/external/Java/JDK/1.6.0/amd64"
export JDK_HOME="$JAVA_HOME"
export PATH="$JDK_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$JDK_HOME/jre/lib/i386:$JDK_HOME/jre/lib/i386/client:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CERNLIB
#--------------------------------------------------------------------------------
export CERN_ROOT="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/CERNLIBS/2006-gfortran"
export CERN="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/CERNLIBS"
export CERN_LEVEL="2006-gfortran"
export CVSCOSRC="$CERN_ROOT/src"
export PATH="$CERN_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$CERN_ROOT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-14/geant4/9.5.p01"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
export G4SYSTEM="Linux-g++"
export LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p1/x86_64-slc5-gcc43-opt/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     QT
#--------------------------------------------------------------------------------
export QTDIR="/afs/cern.ch/sw/lcg/external/qt/4.7.4/x86_64-slc5-gcc43-opt"
export QMAKESPEC="$QTDIR/mkspecs/linux-g++"
export PATH="$QTDIR/bin:$PATH"
export LD_LIBRARY_PATH="$QTDIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/afs/cern.ch/sw/lcg/external/clhep/2.1.1.0/x86_64-slc5-gcc41-opt"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ilcutil/v01-00"
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/lcio/v02-03-01"
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     ROOT
#--------------------------------------------------------------------------------
export ROOTSYS="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ROOT/v5-28-00f"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/gear/v01-02-02"
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/KalTest/v01-05-01"
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     KalDet
#--------------------------------------------------------------------------------
export KALDET="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/KalDet/v01-11"
export LD_LIBRARY_PATH="$KALDET/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Marlin/v01-04"
export PATH="$MARLIN/bin:$PATH"
# export MARLIN_DLL="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinReco/v01-05/lib/libMarlinReco.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraAnalysis/v00-04/lib/libPandoraAnalysis.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinPandora/v00-09-02/lib/libMarlinPandora.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01/lib/libLCFIVertex.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/lib/libCEDViewer.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Overlay/v00-13/lib/libOverlay.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJetClustering/v00-02/lib/libFastJetClustering.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinFastJet/v00-01/lib/libMarlinFastJet.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCTuple/v01-01/lib/libLCTuple.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinKinfit/v00-01-02/lib/libMarlinKinfit.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09/lib/libMarlinTrkProcessors.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Clupatra/v00-09-01/lib/libClupatra.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Eutelescope/v00-08-00-rc1/lib/libEutelescope.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10/lib/libMarlinTPC.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Garlic/v2.10.1/lib/libGarlic.so:$MARLIN_DLL"
#export MARLIN_DLL="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinReco/v01-05/lib/libMarlinReco.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraAnalysis/v00-04/lib/libPandoraAnalysis.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinPandora/v00-09-02/lib/libMarlinPandora.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01/lib/libLCFIVertex.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/lib/libCEDViewer.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Overlay/v00-13/lib/libOverlay.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJetClustering/v00-02/lib/libFastJetClustering.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinFastJet/v00-01/lib/libMarlinFastJet.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCTuple/v01-01/lib/libLCTuple.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinKinfit/v00-01-02/lib/libMarlinKinfit.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09/lib/libMarlinTrkProcessors.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Clupatra/v00-09-01/lib/libClupatra.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10/lib/libMarlinTPC.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Garlic/v2.10.1/lib/libGarlic.so:$MARLIN_DLL"

#export MARLIN_DLL="/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/MarlinReco/v01-05/lib/libMarlinReco.so.1.5.0:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraAnalysis/v00-04/lib/libPandoraAnalysis.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinPandora/v00-09-02/lib/libMarlinPandora.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01/lib/libLCFIVertex.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/lib/libCEDViewer.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Overlay/v00-13/lib/libOverlay.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJetClustering/v00-02/lib/libFastJetClustering.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinFastJet/v00-01/lib/libMarlinFastJet.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCTuple/v01-01/lib/libLCTuple.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinKinfit/v00-01-02/lib/libMarlinKinfit.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09/lib/libMarlinTrkProcessors.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Clupatra/v00-09-01/lib/libClupatra.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10/lib/libMarlinTPC.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Garlic/v2.10.1/lib/libGarlic.so:$MARLIN_DLL"

#export MARLIN_DLL="/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/MarlinReco/v01-05/lib/libMarlinReco.so.1.5.0:/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/MarlinFastJet/lib/libMarlinFastJet.so.0.1.0:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraAnalysis/v00-04/lib/libPandoraAnalysis.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinPandora/v00-09-02/lib/libMarlinPandora.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01/lib/libLCFIVertex.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/lib/libCEDViewer.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Overlay/v00-13/lib/libOverlay.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJetClustering/v00-02/lib/libFastJetClustering.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCTuple/v01-01/lib/libLCTuple.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinKinfit/v00-01-02/lib/libMarlinKinfit.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09/lib/libMarlinTrkProcessors.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Clupatra/v00-09-01/lib/libClupatra.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10/lib/libMarlinTPC.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Garlic/v2.10.1/lib/libGarlic.so:$MARLIN_DLL"

#export MARLIN_DLL="/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/TauFinder/lib/libTauFinder.so.0.1.0:/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/MarlinReco/v01-05/lib/libMarlinReco.so.1.5.0:/afs/cern.ch/user/s/sredford/Documents/ttH/myMarlin/MarlinFastJet/lib/libMarlinFastJet.so.0.1.0:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraAnalysis/v00-04/lib/libPandoraAnalysis.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinPandora/v00-09-02/lib/libMarlinPandora.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01/lib/libLCFIVertex.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/lib/libCEDViewer.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Overlay/v00-13/lib/libOverlay.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJetClustering/v00-02/lib/libFastJetClustering.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCTuple/v01-01/lib/libLCTuple.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinKinfit/v00-01-02/lib/libMarlinKinfit.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09/lib/libMarlinTrkProcessors.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Clupatra/v00-09-01/lib/libClupatra.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10/lib/libMarlinTPC.so:/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Garlic/v2.10.1/lib/libGarlic.so:$MARLIN_DLL"

export LD_LIBRARY_PATH="$MARLIN/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/lccd/v01-02"
export LD_LIBRARY_PATH="$LCCD/lib:/afs/cern.ch/sw/lcg/external/mysql/5.1.45/x86_64-slc5-gcc43-opt/lib/mysql:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export COND_DB_DEBUGLOG="/dev/stdout"
export CondDBMySQL="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CondDBMySQL/CondDBMySQL_ILC-0-9-5"
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/RAIDA/v01-06-02"
export PATH="$RAIDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$RAIDA_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
export MARLINUTIL="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinUtil/v01-05-03"
export LD_LIBRARY_PATH="$MARLINUTIL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/GSL/1.14"
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
export CED="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CED/v01-07"
export PATH="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CED/v01-07/bin:$PATH"
export LD_LIBRARY_PATH="$CED/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Mokka
#--------------------------------------------------------------------------------
export MOKKA="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Mokka/mokka-08-00-03"
export PATH="$MOKKA/bin:$PATH"
export LD_LIBRARY_PATH="$MOKKA/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinReco
#--------------------------------------------------------------------------------
export MARLINRECO="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinReco/v01-05"
export LD_LIBRARY_PATH="$MARLINRECO/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/PandoraPFANew/v00-09"
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinPandora
#--------------------------------------------------------------------------------
export PANDORASETTINGS="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/Config/v01-14-p00/StandardConfig/current/PandoraSettings.xml"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
export LCFIVertex="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIVertex/v00-06-01"
export LD_LIBRARY_PATH="$LCFIVertex/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
export PATH="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/CEDViewer/v01-06-01/bin:$PATH"


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/FastJet/2.4.2"
export PATH="$FastJet_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$FastJet_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinTrk
#--------------------------------------------------------------------------------
export MARLINTRK="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrk/v01-10-01"
export LD_LIBRARY_PATH="$MARLINTRK/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     KiTrack
#--------------------------------------------------------------------------------
export KITRACK="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/KiTrack/v01-04"
export LD_LIBRARY_PATH="$KITRACK/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     KiTrackMarlin
#--------------------------------------------------------------------------------
export KITRACKMARLIN="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/KiTrackMarlin/v01-04"
export LD_LIBRARY_PATH="$KITRACKMARLIN/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinTrkProcessors
#--------------------------------------------------------------------------------
export MARLINTRKPROCESSORS="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTrkProcessors/v01-09"
export LD_LIBRARY_PATH="$MARLINTRKPROCESSORS/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIPlus
#--------------------------------------------------------------------------------
export LCFIPlus="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/LCFIPlus/v00-05-02"
export LD_LIBRARY_PATH="$LCFIPlus/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Eutelescope
#--------------------------------------------------------------------------------
export EUDAQ="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Eutelescope/v00-08-00-rc1/eudaq/trunk"
export MILLEPEDEII="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Eutelescope/v00-08-00-rc1/millepede2/trunk"
export MILLEPEDEII_VERSION="trunk"
export EUDAQ_VERSION="trunk"
export EUTELESCOPE="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Eutelescope/v00-08-00-rc1"
export PATH="$MILLEPEDEII:$EUTELESCOPE/bin:$PATH"
export LD_LIBRARY_PATH="$EUTELESCOPE/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     PathFinder
#--------------------------------------------------------------------------------
export PATHFINDER="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/pathfinder/v00-02"
export LD_LIBRARY_PATH="$PATHFINDER/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MarlinTPC
#--------------------------------------------------------------------------------
export MARLINTPC="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/MarlinTPC/v00-10"
export PATH="$MARLINTPC/bin:$PATH"
export LD_LIBRARY_PATH="$MARLINTPC/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Druid
#--------------------------------------------------------------------------------
export DRUID="/afs/cern.ch/eng/clic/software/x86_64-slc5-gcc41/ILCSOFT/v01-16/Druid/1.8"
export PATH="$DRUID/bin:$PATH"


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
export PATH="/afs/cern.ch/sw/lcg/external/CMake/2.8.6/x86_64-slc5-gcc43-opt/bin:$PATH"

# --- source GEANT4 INIT script ---
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }

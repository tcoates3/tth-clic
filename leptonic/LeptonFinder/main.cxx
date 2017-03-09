#include "tth_analysis.h"

#include <TString.h>

#include <iostream>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>

using namespace std;

int main()
{
  Int_t tpe = 0;

  DIR *dp;

  std::string dir = "/eos/user/y/yixuan/clic/mymarlin/outputs/"; // has 

  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }
  while ((dirp = readdir(dp)) != NULL) {
    std::string file(dirp->d_name);
    if(file[0]==std::string(".")) continue;

    int lastindex = file.find_last_of(".");
    string rawname = file.substr(0, lastindex);

    tpe += leptonInvest(dir+file, "root/"+rawname+".root");

  }
  closedir(dp);
  /*
  */

  //leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2447_97.slcio","tau_investC.root"); // ttZ
  //leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/2441/000/PreMarlinDSTs/PreMarlinDST_2441_99.slcio","leptonInvest_99.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_97.slcio","tau_investC_97.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_96.slcio","tau_investC_96.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_95.slcio","tau_investC_95.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_93.slcio","tau_investC_93.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_89.slcio","tau_investC_89.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_87.slcio","tau_investC_87.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_86.slcio","tau_investC_86.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_85.slcio","tau_investC_85.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_84.slcio","tau_investC_84.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_83.slcio","tau_investC_83.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_82.slcio","tau_investC_82.root"); // ttH
  // leptonInvest("/afs/cern.ch/work/s/sredford/tth_data/PreMarlinDST_2441_81.slcio","tau_investC_81.root"); // ttH

  //postLeptonFinder("/afs/cern.ch/work/s/sredford/tth_data/mytaus0p25_oldold.slcio","postLeptonFinder.root");

}

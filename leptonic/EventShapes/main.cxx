#include "tth_analysis.h"

#include <TString.h>

#include <iostream>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>

using namespace std;

int main()
{

  // reprocessing with all 2441 sample
  DIR *dp;
  std::string dir = "/eos/user/y/yixuan/clic/mymarlin/outputs/";

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

    eventShapes(dir+file,"root/EventShapes_Selected"+rawname+".root","Isolep_Selected","kt_6jets","SelectedPandoraPFOCollection","SelectedPandoraPFOsInJets");
    //eventShapes(dir+file,"2441/000/EventShapes_Tight"+rawname+".root","Isolep_TightSelected","kt_6jets_T_R1p0","TightSelectedPandoraPFOCollection","SelectedPandoraPFOsInJets_T_R1p0");
    //eventShapes(dir+file,"2441/000/EventShapes_Loose"+rawname+".root","Isolep_LooseSelected","kt_6jets_L_R1p0","LooseSelectedPandoraPFOCollection","SelectedPandoraPFOsInJets_L_R1p0");

  }
  closedir(dp);


}

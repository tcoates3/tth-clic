#include "tth_analysis.h"

#include <TString.h>

#include <iostream>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>

#include <ctime>

using namespace std;

int main()
{ 


  Int_t tpe = 0;

  DIR *dp;

  std::string dir = "/eos/user/y/yixuan/clic/mymarlin/";

  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2417/000/FlavourTaggedDSTs/"; // 236493(sl) 235767(had) 236493 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2417/001/FlavourTaggedDSTs/"; // 34040(sl) 33140(had) 34340 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2420/000/FlavourTaggedDSTs/"; // 45620(sl) 46429(had) 46629 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2423/000/FlavourTaggedDSTs/"; // 48062(sl) 48462(had) 48462 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2426/000/FlavourTaggedDSTs/"; // 47432(sl) 47432(had) 47432 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2429/000/FlavourTaggedDSTs/"; // 48031(sl) 48202(had) 48202 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2432/000/FlavourTaggedDSTs/"; // 48770(sl) 48970(had) 48970 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2435/000/FlavourTaggedDSTs/"; // 46853(sl) 46453(had) 46453 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2438/000/FlavourTaggedDSTs/"; // 49569(sl) 49369(had) 49369 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2441/000/FlavourTaggedDSTs/"; // 45761(sl) 45561(had) 45761 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2444/000/FlavourTaggedDSTs/"; // 47423(sl) 47423(had) 47423 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2447/000/FlavourTaggedDSTs/"; // 44359(sl) 44463(had) 47603 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2450/000/FlavourTaggedDSTs/"; // 46965(sl) 47700(had) 49765 (new had)
  //std::string dir = "/afs/cern.ch/work/s/sredford/tth_data/2453/000/FlavourTaggedDSTs/"; // 46594(sl) 47018(had) 47358 (new had)

  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
  }
  
  while ((dirp = readdir(dp)) != NULL) {
    std::string file(dirp->d_name);
    if(file[0]==std::string(".")) continue;

    //if (file.find("_had")!= std::string::npos) continue; // skips had
    //if (file.find("_sl")!= std::string::npos) continue; // skips sl

    int lastindex = file.find_last_of(".");
    string rawname = file.substr(0, lastindex);

    //tpe += treeMaker(dir+file, "2417/000new/"+rawname+".root");    
    //tpe += treeMaker(dir+file, "2417/001new/"+rawname+".root");    
    //tpe += treeMaker(dir+file, "2420/000new/"+rawname+".root");    
    //tpe += treeMaker(dir+file, "2423/000new/"+rawname+".root");    
    //tpe += treeMaker(dir+file, "2426/000new/"+rawname+".root");  
    //tpe += treeMaker(dir+file, "2429/000new/"+rawname+".root");  
    //tpe += treeMaker(dir+file, "2432/000new/"+rawname+".root");  
    //tpe += treeMaker(dir+file, "2435/000new/"+rawname+".root"); 
    //tpe += treeMaker(dir+file, "2438/000new/"+rawname+".root"); 
    //tpe += treeMaker(dir+file, "2441/000new/"+rawname+".root");    
    //tpe += treeMaker(dir+file, "2444/000new/"+rawname+".root"); 
    //tpe += treeMaker(dir+file, "2447/000new/"+rawname+".root");
    //tpe += treeMaker(dir+file, "2450/000new/"+rawname+".root");
    //tpe += treeMaker(dir+file, "2453/000new/"+rawname+".root");

    //tpe += treeMaker(dir+file, "root/"+rawname+".root"); 

    //tpe += treeMaker_hadronic(dir+file, "2417/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2417/001new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2420/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2423/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2426/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2429/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2432/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2435/000new/"+rawname+".root");   
    //tpe += treeMaker_hadronic(dir+file, "2438/000new/"+rawname+".root");   
    //tpe += treeMaker_hadronic(dir+file, "2441/000new/"+rawname+".root");    
    //tpe += treeMaker_hadronic(dir+file, "2444/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2447/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2450/000new/"+rawname+".root");
    //tpe += treeMaker_hadronic(dir+file, "2453/000new/"+rawname+".root");

    }
  //closedir(dp);
  //cout << tpe << endl;

  //treeMaker("/afs/cern.ch/work/s/sredford/tth_data/tth-ln4q-hbb_rec_2441_99_flavourtagged_sl.slcio","test_tree_sl99.root");
  //treeMaker("/afs/cern.ch/work/s/sredford/tth_data/tth-ln4q-hbb_rec_2441_97_flavourtagged_sl.slcio","test_tree_sl97.root");

  //treeMaker_hadronic("/afs/cern.ch/work/s/sredford/tth_data/tth-6q-hbb_rec_2435_99_flavourtagged_had.slcio","test_tree_had99.root");
  //treeMaker_hadronic("/afs/cern.ch/work/s/sredford/tth_data/tt_rec_2417_608_flavourtagged_had.slcio","test_tree_had608.root");
  //treeMaker(dir+"tth_100.slcio","root/tth_100.root");
  treeMaker(dir+"testVertexed.slcio","root/bkg/tt/tt.root");
  /*
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_103.slcio_vertexed.slcio","tth_103.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_105.slcio_vertexed.slcio","tth_105.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_106.slcio_vertexed.slcio","tth_106.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_107.slcio_vertexed.slcio","tth_107.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_108.slcio_vertexed.slcio","tth_108.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_109.slcio_vertexed.slcio","tth_109.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_10.slcio_vertexed.slcio","tth_10.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_111.slcio_vertexed.slcio","tth_111.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_114.slcio_vertexed.slcio","tth_114.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_115.slcio_vertexed.slcio","tth_115.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_116.slcio_vertexed.slcio","tth_116.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_117.slcio_vertexed.slcio","tth_117.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_118.slcio_vertexed.slcio","tth_118.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_121.slcio_vertexed.slcio","tth_121.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_122.slcio_vertexed.slcio","tth_122.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_123.slcio_vertexed.slcio","tth_123.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_125.slcio_vertexed.slcio","tth_125.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_126.slcio_vertexed.slcio","tth_126.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_127.slcio_vertexed.slcio","tth_127.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_128.slcio_vertexed.slcio","tth_128.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_129.slcio_vertexed.slcio","tth_129.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_130.slcio_vertexed.slcio","tth_130.root");
  treeMaker(dir+"tth-ln4q-hbb_rec_2441_131.slcio_vertexed.slcio","tth_131.root");*/
}

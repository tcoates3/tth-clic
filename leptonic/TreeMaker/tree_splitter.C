// code to take a merged tree and split it in two for two independent TMVA samples
// root [0] .L tree_splitter.C
// root [1] tree_splitter()


void tree_splitter() {
  
  // sl
  /*
  TFile* full_file_2417a_sl = new TFile("2417/000new/merged_tree_2417_000_sl.root","READ");
  splitter(full_file_2417a_sl);

  TFile* full_file_2417b_sl = new TFile("2417/001new/merged_tree_2417_001_sl.root","READ");
  splitter(full_file_2417b_sl);

  TFile* full_file_2420_sl = new TFile("2420/000new/merged_tree_2420_000_sl.root","READ");
  splitter(full_file_2420_sl);

  TFile* full_file_2423_sl = new TFile("2423/000new/merged_tree_2423_000_sl.root","READ");
  splitter(full_file_2423_sl);

  TFile* full_file_2426_sl = new TFile("2426/000new/merged_tree_2426_000_sl.root","READ");
  splitter(full_file_2426_sl);

  TFile* full_file_2429_sl = new TFile("2429/000new/merged_tree_2429_000_sl.root","READ");
  splitter(full_file_2429_sl);

  TFile* full_file_2432_sl = new TFile("2432/000new/merged_tree_2432_000_sl.root","READ");
  splitter(full_file_2432_sl);

  TFile* full_file_2435_sl = new TFile("2435/000new/merged_tree_2435_000_sl.root","READ");
  splitter(full_file_2435_sl);

  TFile* full_file_2438_sl = new TFile("2438/000new/merged_tree_2438_000_sl.root","READ");
  splitter(full_file_2438_sl);

  TFile* full_file_2441_sl = new TFile("2441/000new/merged_tree_2441_000_sl.root","READ");
  splitter(full_file_2441_sl);

  TFile* full_file_2444_sl = new TFile("2444/000new/merged_tree_2444_000_sl.root","READ");
  splitter(full_file_2444_sl);

  TFile* full_file_2447_sl = new TFile("2447/000new/merged_tree_2447_000_sl.root","READ");
  splitter(full_file_2447_sl);

  TFile* full_file_2450_sl = new TFile("2450/000new/merged_tree_2450_000_sl.root","READ");
  splitter(full_file_2450_sl);

  TFile* full_file_2453_sl = new TFile("2453/000new/merged_tree_2453_000_sl.root","READ");
  splitter(full_file_2453_sl);
  */

  // had
  TFile* full_file_2417a_had = new TFile("2417/000new/merged_tree_2417_000_had.root","READ");
  splitter(full_file_2417a_had);

  TFile* full_file_2417b_had = new TFile("2417/001new/merged_tree_2417_001_had.root","READ");
  splitter(full_file_2417b_had);

  TFile* full_file_2420_had = new TFile("2420/000new/merged_tree_2420_000_had.root","READ");
  splitter(full_file_2420_had);

  TFile* full_file_2423_had = new TFile("2423/000new/merged_tree_2423_000_had.root","READ");
  splitter(full_file_2423_had);

  TFile* full_file_2426_had = new TFile("2426/000new/merged_tree_2426_000_had.root","READ");
  splitter(full_file_2426_had);

  TFile* full_file_2429_had = new TFile("2429/000new/merged_tree_2429_000_had.root","READ");
  splitter(full_file_2429_had);

  TFile* full_file_2432_had = new TFile("2432/000new/merged_tree_2432_000_had.root","READ");
  splitter(full_file_2432_had);

  TFile* full_file_2435_had = new TFile("2435/000new/merged_tree_2435_000_had.root","READ");
  splitter(full_file_2435_had);

  TFile* full_file_2438_had = new TFile("2438/000new/merged_tree_2438_000_had.root","READ");
  splitter(full_file_2438_had);

  TFile* full_file_2441_had = new TFile("2441/000new/merged_tree_2441_000_had.root","READ");
  splitter(full_file_2441_had);

  TFile* full_file_2444_had = new TFile("2444/000new/merged_tree_2444_000_had.root","READ");
  splitter(full_file_2444_had);

  TFile* full_file_2447_had = new TFile("2447/000new/merged_tree_2447_000_had.root","READ");
  splitter(full_file_2447_had);

  TFile* full_file_2450_had = new TFile("2450/000new/merged_tree_2450_000_had.root","READ");
  splitter(full_file_2450_had);

  TFile* full_file_2453_had = new TFile("2453/000new/merged_tree_2453_000_had.root","READ");
  splitter(full_file_2453_had);

}

void splitter(TFile *full_file){

  TTree* full_tree = full_file->Get("tth_tree");

  string full_file_name = full_file->GetName();
  int lastindex = full_file_name.find_last_of(".");
  TString rawname = full_file_name.substr(0, lastindex);

  half1_file_name = rawname+"_half1.root";
  half2_file_name = rawname+"_half2.root";

  TFile* half1_file = new TFile(half1_file_name, "RECREATE");
  TTree* half1_tree = full_tree->CloneTree(0);

  TFile* half2_file = new TFile(half2_file_name, "RECREATE");
  TTree* half2_tree = full_tree->CloneTree(0);

  TRandom *r3 = new TRandom3();

  for (i=0; i<full_tree->GetEntries(); i++){
    full_tree->GetEntry(i);
    Double_t rndm = r3->Rndm();

    if (rndm < 0.5) {
      half1_tree->Fill();
    } 
    else {
      half2_tree->Fill();
    }
  }
  half1_file->Write();
  half2_file->Write();

  cout << full_file_name << " entries " << full_tree->GetEntries() << endl;
  cout << half1_file_name << " entries " << half1_tree->GetEntries() << endl;
  cout << half2_file_name << " entries " << half2_tree->GetEntries() << endl;
  cout << "1+2 entries: " << half1_tree->GetEntries() + half2_tree->GetEntries() << endl;
  cout << " " << endl;
}


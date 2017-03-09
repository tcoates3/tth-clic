# file to plot nLeptons seperately for each sample, to see how the isolated lepton finder is performing

from util import *

# preliminary
'''
tup_2417 = TChain("tth_tree")
tup_2417.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/000/merged_tree_2417_000.root")
tup_2417.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/001/merged_tree_2417_001.root")

file_2420_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2420/000/merged_tree_2420_000.root") 
file_2423_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2423/000/merged_tree_2423_000.root") 
file_2426_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2426/000/merged_tree_2426_000.root") 
file_2429_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2429/000/merged_tree_2429_000.root") 
file_2432_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2432/000/merged_tree_2432_000.root") 
file_2435_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2435/000/merged_tree_2435_000.root") 
file_2438_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2438/000/merged_tree_2438_000.root") 
file_2441_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2441/000/merged_tree_2441_000.root") 
file_2444_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2444/000/merged_tree_2444_000.root") 
file_2447_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2447/000/merged_tree_2447_000.root") 
file_2450_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2450/000/merged_tree_2450_000.root") 
file_2453_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2453/000/merged_tree_2453_000.root") 
'''

# sl
tup_2417 = TChain("tth_tree")
tup_2417.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/000new/merged_tree_2417_000_sl.root")
tup_2417.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/001new/merged_tree_2417_001_sl.root")

file_2420_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2420/000new/merged_tree_2420_000_sl.root") 
file_2423_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2423/000new/merged_tree_2423_000_sl.root") 
file_2426_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2426/000new/merged_tree_2426_000_sl.root") 
file_2429_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2429/000new/merged_tree_2429_000_sl.root") 
file_2432_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2432/000new/merged_tree_2432_000_sl.root") 
file_2435_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2435/000new/merged_tree_2435_000_sl.root") 
file_2438_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2438/000new/merged_tree_2438_000_sl.root") 
file_2441_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2441/000new/merged_tree_2441_000_sl.root") 
file_2444_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2444/000new/merged_tree_2444_000_sl.root") 
file_2447_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2447/000new/merged_tree_2447_000_sl.root") 
file_2450_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2450/000new/merged_tree_2450_000_sl.root") 
file_2453_000 = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2453/000new/merged_tree_2453_000_sl.root") 

tup_2420 = file_2420_000.Get("tth_tree")
tup_2423 = file_2423_000.Get("tth_tree")
tup_2426 = file_2426_000.Get("tth_tree")
tup_2429 = file_2429_000.Get("tth_tree")
tup_2432 = file_2432_000.Get("tth_tree")
tup_2435 = file_2435_000.Get("tth_tree")
tup_2438 = file_2438_000.Get("tth_tree")
tup_2441 = file_2441_000.Get("tth_tree")
tup_2444 = file_2444_000.Get("tth_tree")
tup_2447 = file_2447_000.Get("tth_tree")
tup_2450 = file_2450_000.Get("tth_tree")
tup_2453 = file_2453_000.Get("tth_tree")

c = TCanvas("c","c",0,0,400,400)

for tup,name in zip([tup_2417,tup_2420,tup_2423,tup_2426,tup_2429,tup_2432,tup_2435,tup_2438,tup_2441,tup_2444,tup_2447,tup_2450,tup_2453],["2417","2420","2423","2426","2429","2432","2435","2438","2441","2444","2447","2450","2453"]):
    print "Doing", tup.GetName(), name

    hist = TH1F("hist","",6,-0.5,5.5)

    tup.Draw("nLeptons >> hist")

    hist.Scale(1./hist.GetEntries())

    hist.SetLabelSize(0.05)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitleOffset(0.9)
    hist.GetYaxis().SetLabelSize(0.05)
    hist.GetXaxis().SetTitle("nLeptons")
    hist.SetNdivisions(505)
    hist.SetStats(0)
    hist.SetLineWidth(3)
    hist.Draw()
    
    if name == "2417": hist.SetTitle("tt")
    if name == "2420": hist.SetTitle("ttbb 2l2nbb all")
    if name == "2423": hist.SetTitle("ttbb 6q all")
    if name == "2426": hist.SetTitle("ttbb ln4q all")
    if name == "2429": hist.SetTitle("ttH 2l2nbb Hbb")
    if name == "2432": hist.SetTitle("ttH 2l2nbb Hnonbb")
    if name == "2435": hist.SetTitle("ttH 6q Hbb")
    if name == "2438": hist.SetTitle("ttH 6q Hnonbb")
    if name == "2441": hist.SetTitle("ttH ln4q Hbb")
    if name == "2444": hist.SetTitle("ttH ln4q Hnonbb")
    if name == "2447": hist.SetTitle("ttZ 2l2nbb all")
    if name == "2450": hist.SetTitle("ttZ 6q all")
    if name == "2453": hist.SetTitle("ttZ ln4q all")

    c.Update()
    c.SaveAs("plots/nLeptons/nLeptons_%s.png" %name)
    c.SaveAs("plots/nLeptons/nLeptons_%s.pdf" %name)

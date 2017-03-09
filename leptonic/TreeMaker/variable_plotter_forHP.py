# to plot the variables produced in my ntuple making
# for the Higgs Paper

# usng the CLICdp ROOT style
from ROOT import *
gROOT.ProcessLine(".L /afs/cern.ch/user/s/sredford/Documents/RootStyle/SampleStyle/rootstyle/CLICStyle.C")
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

# sl
tup_2417_sl = TChain("tth_tree")
tup_2417_sl.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/000new/merged_tree_2417_000_sl.root")
tup_2417_sl.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/001new/merged_tree_2417_001_sl.root")

file_2420_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2420/000new/merged_tree_2420_000_sl.root") 
file_2423_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2423/000new/merged_tree_2423_000_sl.root") 
file_2426_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2426/000new/merged_tree_2426_000_sl.root") 
file_2429_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2429/000new/merged_tree_2429_000_sl.root") 
file_2432_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2432/000new/merged_tree_2432_000_sl.root") 
file_2435_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2435/000new/merged_tree_2435_000_sl.root") 
file_2438_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2438/000new/merged_tree_2438_000_sl.root") 
file_2441_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2441/000new/merged_tree_2441_000_sl.root") 
file_2444_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2444/000new/merged_tree_2444_000_sl.root") 
file_2447_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2447/000new/merged_tree_2447_000_sl.root") 
file_2450_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2450/000new/merged_tree_2450_000_sl.root") 
file_2453_000_sl = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2453/000new/merged_tree_2453_000_sl.root") 

tup_2420_sl = file_2420_000_sl.Get("tth_tree")
tup_2423_sl = file_2423_000_sl.Get("tth_tree")
tup_2426_sl = file_2426_000_sl.Get("tth_tree")
tup_2429_sl = file_2429_000_sl.Get("tth_tree")
tup_2432_sl = file_2432_000_sl.Get("tth_tree")
tup_2435_sl = file_2435_000_sl.Get("tth_tree")
tup_2438_sl = file_2438_000_sl.Get("tth_tree")
tup_2441_sl = file_2441_000_sl.Get("tth_tree")
tup_2444_sl = file_2444_000_sl.Get("tth_tree")
tup_2447_sl = file_2447_000_sl.Get("tth_tree")
tup_2450_sl = file_2450_000_sl.Get("tth_tree")
tup_2453_sl = file_2453_000_sl.Get("tth_tree")

# had
tup_2417_had = TChain("tth_tree")
tup_2417_had.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/000new/merged_tree_2417_000_had.root")
tup_2417_had.Add("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2417/001new/merged_tree_2417_001_had.root")

file_2420_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2420/000new/merged_tree_2420_000_had.root") 
file_2423_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2423/000new/merged_tree_2423_000_had.root") 
file_2426_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2426/000new/merged_tree_2426_000_had.root") 
file_2429_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2429/000new/merged_tree_2429_000_had.root") 
file_2432_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2432/000new/merged_tree_2432_000_had.root") 
file_2435_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2435/000new/merged_tree_2435_000_had.root") 
file_2438_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2438/000new/merged_tree_2438_000_had.root") 
file_2441_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2441/000new/merged_tree_2441_000_had.root") 
file_2444_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2444/000new/merged_tree_2444_000_had.root") 
file_2447_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2447/000new/merged_tree_2447_000_had.root") 
file_2450_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2450/000new/merged_tree_2450_000_had.root") 
file_2453_000_had = TFile("/afs/cern.ch/user/s/sredford/Documents/ttH/TreeMaker/2453/000new/merged_tree_2453_000_had.root") 

tup_2420_had = file_2420_000_had.Get("tth_tree")
tup_2423_had = file_2423_000_had.Get("tth_tree")
tup_2426_had = file_2426_000_had.Get("tth_tree")
tup_2429_had = file_2429_000_had.Get("tth_tree")
tup_2432_had = file_2432_000_had.Get("tth_tree")
tup_2435_had = file_2435_000_had.Get("tth_tree")
tup_2438_had = file_2438_000_had.Get("tth_tree")
tup_2441_had = file_2441_000_had.Get("tth_tree")
tup_2444_had = file_2444_000_had.Get("tth_tree")
tup_2447_had = file_2447_000_had.Get("tth_tree")
tup_2450_had = file_2450_000_had.Get("tth_tree")
tup_2453_had = file_2453_000_had.Get("tth_tree")

def main():
    global c, sh, bh
    gROOT.SetBatch(True)

    PlotForHiggsPaper("had","h_bb_mass",50,0,200, "Higgs candidate mass (GeV/c^{2})","nLeptons == 0 && nJets == 8","h_bb_mass")
    PlotForHiggsPaper("had","btag_2",50,0,1, "Third highest b-tag","nLeptons == 0 && nJets == 8","btag_2")
    PlotDecayAnglesForHiggsPaper("had","dec_ang_min","dec_ang_max",40,0,3.15,"Higgs decay angle (rad)","nLeptons == 0 && nJets == 8","dec_ang")


def PlotDecayAnglesForHiggsPaper(channel,var1,var2,nbins,lbin,ubin,axisname,cut,savename):
    # only valid for hadronic channel
    print var1, var2

    c = TCanvas("c","c",0,0,400,400)

    # by channel
    h_tt = TH1F("h_tt","",nbins,lbin,ubin) # tt
    h_ttbb = TH1F("h_ttbb","",nbins,lbin,ubin) # ttbb
    h_tth = TH1F("h_tth","",nbins,lbin,ubin) # other tth
    h_sig = TH1F("h_sig","",nbins,lbin,ubin) # signal
    h_ttz = TH1F("h_ttz","",nbins,lbin,ubin) # ttz

    tup_2417_had.Draw("%s >> h_tt" %var1, "%s" %cut)
    tup_2420_had.Draw("%s >> h_ttbb" %var1, "%s" %cut)
    tup_2423_had.Draw("%s >> +h_ttbb" %var1, "%s" %cut)
    tup_2426_had.Draw("%s >> +h_ttbb" %var1, "%s" %cut)
    tup_2429_had.Draw("%s >> h_tth" %var1, "%s" %cut)
    tup_2432_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2435_had.Draw("%s >> h_sig" %var1, "%s" %cut)
    tup_2438_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2441_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2444_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2447_had.Draw("%s >> h_ttz" %var1, "%s" %cut)
    tup_2450_had.Draw("%s >> +h_ttz" %var1, "%s" %cut)
    tup_2453_had.Draw("%s >> +h_ttz" %var1, "%s" %cut)
    
    tup_2417_had.Draw("%s >> +h_tt" %var2, "%s" %cut)
    tup_2420_had.Draw("%s >> +h_ttbb" %var2, "%s" %cut)
    tup_2423_had.Draw("%s >> +h_ttbb" %var2, "%s" %cut)
    tup_2426_had.Draw("%s >> +h_ttbb" %var2, "%s" %cut)
    tup_2429_had.Draw("%s >> +h_tth" %var2, "%s" %cut)
    tup_2432_had.Draw("%s >> +h_tth" %var2, "%s" %cut)
    tup_2435_had.Draw("%s >> +h_sig" %var2, "%s" %cut)
    tup_2438_had.Draw("%s >> +h_tth" %var2, "%s" %cut)
    tup_2441_had.Draw("%s >> +h_tth" %var2, "%s" %cut)
    tup_2444_had.Draw("%s >> +h_tth" %var2, "%s" %cut)
    tup_2447_had.Draw("%s >> +h_ttz" %var2, "%s" %cut)
    tup_2450_had.Draw("%s >> +h_ttz" %var2, "%s" %cut)
    tup_2453_had.Draw("%s >> +h_ttz" %var2, "%s" %cut)

    for hist in [h_tt, h_ttbb, h_tth, h_sig, h_ttz]:
        if (hist.Integral() > 0):
            hist.Scale(1/hist.Integral())

    h_tt.SetLineColor(kRed)
    h_ttbb.SetLineColor(80)
    h_tth.SetLineColor(61)
    h_ttz.SetLineColor(91)
    h_sig.SetLineColor(kBlack)
    h_sig.SetFillStyle(3004)
    h_sig.SetFillColor(kBlack)

    max_height = 0
    for hist in [h_tt, h_ttbb, h_tth, h_sig, h_ttz]:
        if hist.GetMaximum() > max_height:
            max_height = hist.GetMaximum()

    h_sig.GetYaxis().SetRangeUser(0,1.4*max_height)

    h_sig.SetTitle("")
    h_sig.GetXaxis().SetTitleOffset(0.9)
    h_sig.GetXaxis().SetTitle("%s" %axisname)
    h_sig.Draw()
    h_tt.Draw("same")
    h_ttbb.Draw("same")
    h_tth.Draw("same")
    h_ttz.Draw("same")

    leg = TLegend(0.12,0.68,0.6,0.89)
    leg.SetFillStyle(0)
    leg.AddEntry(h_sig,"t#bar{t}H, fully hadronic, H #rightarrow b#bar{b}","f")
    leg.AddEntry(h_tth,"Other t#bar{t}H","l")
    leg.AddEntry(h_ttbb,"t#bar{t}b#bar{b}","l")
    leg.SetTextSize(0.05)
    leg.Draw()

    leg2 = TLegend(0.5,0.68,0.89,0.83)
    leg2.SetFillStyle(0)
    leg2.AddEntry(h_ttz,"t#bar{t}Z","l")
    leg2.AddEntry(h_tt,"t#bar{t}","l")
    leg2.SetTextSize(0.05)
    leg2.Draw()

    c.Update()
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.png" %(savename,channel))
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.pdf" %(savename,channel))
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.C" %(savename,channel)) 

def PlotForHiggsPaper(channel,var1,nbins,lbin,ubin,axisname,cut,savename):
    # only valid for hadronic channel
    print var1

    c = TCanvas("c","c",0,0,400,400)

    # by channel
    h_tt = TH1F("h_tt","",nbins,lbin,ubin) # tt
    h_ttbb = TH1F("h_ttbb","",nbins,lbin,ubin) # ttbb
    h_tth = TH1F("h_tth","",nbins,lbin,ubin) # other tth
    h_sig = TH1F("h_sig","",nbins,lbin,ubin) # signal
    h_ttz = TH1F("h_ttz","",nbins,lbin,ubin) # ttz

    tup_2417_had.Draw("%s >> h_tt" %var1, "%s" %cut)
    tup_2420_had.Draw("%s >> h_ttbb" %var1, "%s" %cut)
    tup_2423_had.Draw("%s >> +h_ttbb" %var1, "%s" %cut)
    tup_2426_had.Draw("%s >> +h_ttbb" %var1, "%s" %cut)
    tup_2429_had.Draw("%s >> h_tth" %var1, "%s" %cut)
    tup_2432_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2435_had.Draw("%s >> h_sig" %var1, "%s" %cut)
    tup_2438_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2441_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2444_had.Draw("%s >> +h_tth" %var1, "%s" %cut)
    tup_2447_had.Draw("%s >> h_ttz" %var1, "%s" %cut)
    tup_2450_had.Draw("%s >> +h_ttz" %var1, "%s" %cut)
    tup_2453_had.Draw("%s >> +h_ttz" %var1, "%s" %cut)

    for hist in [h_tt, h_ttbb, h_tth, h_sig, h_ttz]:
        if (hist.Integral() > 0):
            hist.Scale(1/hist.Integral())

    h_tt.SetLineColor(kRed)
    h_ttbb.SetLineColor(80)
    h_tth.SetLineColor(61)
    h_sig.SetLineColor(kBlack)
    h_sig.SetFillStyle(3004)
    h_sig.SetFillColor(kBlack)
    h_ttz.SetLineColor(91)

    max_height = 0
    for hist in [h_tt, h_ttbb, h_tth, h_sig, h_ttz]:
        if hist.GetMaximum() > max_height:
            max_height = hist.GetMaximum()

    h_sig.GetYaxis().SetRangeUser(0,1.2*max_height)

    h_sig.SetTitle("")
    h_sig.GetXaxis().SetTitleOffset(0.9)
    h_sig.GetXaxis().SetTitle("%s" %axisname)
    h_sig.Draw()
    h_ttbb.Draw("same")
    h_tth.Draw("same")
    h_ttz.Draw("same")
    h_tt.Draw("same")

    if var1 == "btag_2":
        leg = TLegend(0.22,0.54,0.89,0.89)
    else:
        leg = TLegend(0.12,0.54,0.6,0.89)
    leg.SetFillStyle(0)
    leg.AddEntry(h_sig,"t#bar{t}H, fully hadronic, H #rightarrow b#bar{b}","f")
    leg.AddEntry(h_tth,"Other t#bar{t}H","l")
    leg.AddEntry(h_ttbb,"t#bar{t}b#bar{b}","l")
    leg.AddEntry(h_ttz,"t#bar{t}Z","l")
    leg.AddEntry(h_tt,"t#bar{t}","l")
    leg.SetTextSize(0.05)
    leg.Draw()

    c.Update()
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.png" %(savename,channel))
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.pdf" %(savename,channel))
    c.SaveAs("plots/forHiggsPaper/%s_%s_forHP.C" %(savename,channel)) 

if __name__ == "__main__":
    main()

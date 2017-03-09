# to plot event shape variables

from util import *
gStyle.SetOptStat(0)
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(11)
gStyle.SetTitleSize(0.05,"xyz")
gStyle.SetTitleSize(0.05,"t")
gStyle.SetLabelSize(0.05,"xyz")
gStyle.SetHistLineWidth(3)
gStyle.SetFuncColor(kBlue)
gStyle.SetFuncWidth(2)
gStyle.SetLegendFont(22)
gStyle.SetLineWidth(2)

thrust_hists = []
sphericity_hists = []
aplanarity_hists = []
oblateness_hists = []

c = TCanvas("c","",0,0,400,400)

selectedfile = TFile("2441/000/merged_tree_EventShape_SelectedVertexed_2441.root")
selectedtup = selectedfile.Get("tth_tree")

looseselectedfile = TFile("2441/000/merged_tree_EventShape_LooseVertexed_2441.root")
looseselectedtup = looseselectedfile.Get("tth_tree")

tightselectedfile = TFile("2441/000/merged_tree_EventShape_TightVertexed_2441.root")
tightselectedtup = tightselectedfile.Get("tth_tree")

#thrust
selected_thrust_hist = TH1F("selected_thrust_hist","",20,0.5,1)
selectedtup.Draw("pfo_thrust >> selected_thrust_hist")
thrust_hists.append(selected_thrust_hist)

looseselected_thrust_hist = TH1F("looseselected_thrust_hist","",20,0.5,1)
looseselectedtup.Draw("pfo_thrust >> looseselected_thrust_hist")
thrust_hists.append(looseselected_thrust_hist)

tightselected_thrust_hist = TH1F("tightselected_thrust_hist","",20,0.5,1)
tightselectedtup.Draw("pfo_thrust >> tightselected_thrust_hist")
thrust_hists.append(tightselected_thrust_hist)

thrust_hists[0].Scale(1./thrust_hists[0].GetEntries())
thrust_hists[1].Scale(1./thrust_hists[1].GetEntries())
thrust_hists[2].Scale(1./thrust_hists[2].GetEntries())

thrust_hists[0].SetTitle("")
thrust_hists[0].GetXaxis().SetTitle("Thrust")
thrust_hists[0].Draw()
thrust_hists[1].SetLineColor(kGreen+2)
thrust_hists[1].Draw("same")
thrust_hists[2].SetLineColor(kRed)
thrust_hists[2].Draw("same")
thrust_hists[0].Draw("same")

leg = TLegend(0.15,0.6,0.45,0.89)
leg.AddEntry(thrust_hists[1],"#color[418]{Loose}","l")
leg.AddEntry(thrust_hists[0],"Selected","l")
leg.AddEntry(thrust_hists[2],"#color[2]{Tight}","l")
leg.SetLineColor(kWhite)
leg.SetFillColor(kWhite)
leg.SetTextSize(0.05)
leg.Draw()
c.SaveAs("plots/thrust_allPFOs.png")
c.SaveAs("plots/thrust_allPFOs.pdf")

#sphericity
selected_sphericity_hist = TH1F("selected_sphericity_hist","",20,0,1)
selectedtup.Draw("pfo_sphericity >> selected_sphericity_hist")
sphericity_hists.append(selected_sphericity_hist)

looseselected_sphericity_hist = TH1F("looseselected_sphericity_hist","",20,0,1)
looseselectedtup.Draw("pfo_sphericity >> looseselected_sphericity_hist")
sphericity_hists.append(looseselected_sphericity_hist)

tightselected_sphericity_hist = TH1F("tightselected_sphericity_hist","",20,0,1)
tightselectedtup.Draw("pfo_sphericity >> tightselected_sphericity_hist")
sphericity_hists.append(tightselected_sphericity_hist)

sphericity_hists[0].Scale(1./sphericity_hists[0].GetEntries())
sphericity_hists[1].Scale(1./sphericity_hists[1].GetEntries())
sphericity_hists[2].Scale(1./sphericity_hists[2].GetEntries())

sphericity_hists[0].SetTitle("")
sphericity_hists[0].GetXaxis().SetTitle("Sphericity")
sphericity_hists[0].Draw()
sphericity_hists[1].SetLineColor(kGreen+2)
sphericity_hists[1].Draw("same")
sphericity_hists[2].SetLineColor(kRed)
sphericity_hists[2].Draw("same")
sphericity_hists[0].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/sphericity_allPFOs.png")
c.SaveAs("plots/sphericity_allPFOs.pdf")

#aplanarity
selected_aplanarity_hist = TH1F("selected_aplanarity_hist","",20,0,0.5)
selectedtup.Draw("pfo_aplanarity >> selected_aplanarity_hist")
aplanarity_hists.append(selected_aplanarity_hist)

looseselected_aplanarity_hist = TH1F("looseselected_aplanarity_hist","",20,0,0.5)
looseselectedtup.Draw("pfo_aplanarity >> looseselected_aplanarity_hist")
aplanarity_hists.append(looseselected_aplanarity_hist)

tightselected_aplanarity_hist = TH1F("tightselected_aplanarity_hist","",20,0,0.5)
tightselectedtup.Draw("pfo_aplanarity >> tightselected_aplanarity_hist")
aplanarity_hists.append(tightselected_aplanarity_hist)

aplanarity_hists[0].Scale(1./aplanarity_hists[0].GetEntries())
aplanarity_hists[1].Scale(1./aplanarity_hists[1].GetEntries())
aplanarity_hists[2].Scale(1./aplanarity_hists[2].GetEntries())

aplanarity_hists[0].SetTitle("")
aplanarity_hists[0].GetXaxis().SetTitle("Aplanarity")
aplanarity_hists[0].Draw()
aplanarity_hists[1].SetLineColor(kGreen+2)
aplanarity_hists[1].Draw("same")
aplanarity_hists[2].SetLineColor(kRed)
aplanarity_hists[2].Draw("same")
aplanarity_hists[0].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/aplanarity_allPFOs.png")
c.SaveAs("plots/aplanarity_allPFOs.pdf")

#oblateness
selected_oblateness_hist = TH1F("selected_oblateness_hist","",20,-0.5,1)
selectedtup.Draw("pfo_oblateness >> selected_oblateness_hist")
oblateness_hists.append(selected_oblateness_hist)

looseselected_oblateness_hist = TH1F("looseselected_oblateness_hist","",20,-0.5,1)
looseselectedtup.Draw("pfo_oblateness >> looseselected_oblateness_hist")
oblateness_hists.append(looseselected_oblateness_hist)

tightselected_oblateness_hist = TH1F("tightselected_oblateness_hist","",20,-0.5,1)
tightselectedtup.Draw("pfo_oblateness >> tightselected_oblateness_hist")
oblateness_hists.append(tightselected_oblateness_hist)

oblateness_hists[0].Scale(1./oblateness_hists[0].GetEntries())
oblateness_hists[1].Scale(1./oblateness_hists[1].GetEntries())
oblateness_hists[2].Scale(1./oblateness_hists[2].GetEntries())

oblateness_hists[0].SetTitle("")
oblateness_hists[0].GetXaxis().SetTitle("Oblateness")
oblateness_hists[0].Draw()
oblateness_hists[1].SetLineColor(kGreen+2)
oblateness_hists[1].Draw("same")
oblateness_hists[2].SetLineColor(kRed)
oblateness_hists[2].Draw("same")
oblateness_hists[0].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/oblateness_allPFOs.png")
c.SaveAs("plots/oblateness_allPFOs.pdf")

## PFOS in Jets
#thrust
selected_thrust_hist = TH1F("selected_thrust_hist","",20,0.5,1)
selectedtup.Draw("pfos_injets_thrust >> selected_thrust_hist")
thrust_hists.append(selected_thrust_hist)

looseselected_thrust_hist = TH1F("looseselected_thrust_hist","",20,0.5,1)
looseselectedtup.Draw("pfos_injets_thrust >> looseselected_thrust_hist")
thrust_hists.append(looseselected_thrust_hist)

tightselected_thrust_hist = TH1F("tightselected_thrust_hist","",20,0.5,1)
tightselectedtup.Draw("pfos_injets_thrust >> tightselected_thrust_hist")
thrust_hists.append(tightselected_thrust_hist)

thrust_hists[3].Scale(1./thrust_hists[3].GetEntries())
thrust_hists[4].Scale(1./thrust_hists[4].GetEntries())
thrust_hists[5].Scale(1./thrust_hists[5].GetEntries())

thrust_hists[3].SetTitle("")
thrust_hists[3].GetXaxis().SetTitle("Thrust")
thrust_hists[3].Draw()
thrust_hists[4].SetLineColor(kGreen+2)
thrust_hists[4].Draw("same")
thrust_hists[5].SetLineColor(kRed)
thrust_hists[5].Draw("same")
thrust_hists[3].Draw("same")

leg.SetX1NDC(0.15)
leg.SetX2NDC(0.45)
leg.Draw()
c.SaveAs("plots/thrust_PFOsinJets.png")
c.SaveAs("plots/thrust_PFOsinJets.pdf")

#sphericity
selected_sphericity_hist = TH1F("selected_sphericity_hist","",20,0,1)
selectedtup.Draw("pfos_injets_sphericity >> selected_sphericity_hist")
sphericity_hists.append(selected_sphericity_hist)

looseselected_sphericity_hist = TH1F("looseselected_sphericity_hist","",20,0,1)
looseselectedtup.Draw("pfos_injets_sphericity >> looseselected_sphericity_hist")
sphericity_hists.append(looseselected_sphericity_hist)

tightselected_sphericity_hist = TH1F("tightselected_sphericity_hist","",20,0,1)
tightselectedtup.Draw("pfos_injets_sphericity >> tightselected_sphericity_hist")
sphericity_hists.append(tightselected_sphericity_hist)

sphericity_hists[3].Scale(1./sphericity_hists[3].GetEntries())
sphericity_hists[4].Scale(1./sphericity_hists[4].GetEntries())
sphericity_hists[5].Scale(1./sphericity_hists[5].GetEntries())

sphericity_hists[3].SetTitle("")
sphericity_hists[3].GetXaxis().SetTitle("Sphericity")
sphericity_hists[3].Draw()
sphericity_hists[4].SetLineColor(kGreen+2)
sphericity_hists[4].Draw("same")
sphericity_hists[5].SetLineColor(kRed)
sphericity_hists[5].Draw("same")
sphericity_hists[3].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/sphericity_PFOsinJets.png")
c.SaveAs("plots/sphericity_PFOsinJets.pdf")

#aplanarity
selected_aplanarity_hist = TH1F("selected_aplanarity_hist","",20,0,0.5)
selectedtup.Draw("pfos_injets_aplanarity >> selected_aplanarity_hist")
aplanarity_hists.append(selected_aplanarity_hist)

looseselected_aplanarity_hist = TH1F("looseselected_aplanarity_hist","",20,0,0.5)
looseselectedtup.Draw("pfos_injets_aplanarity >> looseselected_aplanarity_hist")
aplanarity_hists.append(looseselected_aplanarity_hist)

tightselected_aplanarity_hist = TH1F("tightselected_aplanarity_hist","",20,0,0.5)
tightselectedtup.Draw("pfos_injets_aplanarity >> tightselected_aplanarity_hist")
aplanarity_hists.append(tightselected_aplanarity_hist)

aplanarity_hists[3].Scale(1./aplanarity_hists[3].GetEntries())
aplanarity_hists[4].Scale(1./aplanarity_hists[4].GetEntries())
aplanarity_hists[5].Scale(1./aplanarity_hists[5].GetEntries())

aplanarity_hists[3].SetTitle("")
aplanarity_hists[3].GetXaxis().SetTitle("Aplanarity")
aplanarity_hists[3].Draw()
aplanarity_hists[4].SetLineColor(kGreen+2)
aplanarity_hists[4].Draw("same")
aplanarity_hists[5].SetLineColor(kRed)
aplanarity_hists[5].Draw("same")
aplanarity_hists[3].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/aplanarity_PFOsinJets.png")
c.SaveAs("plots/aplanarity_PFOsinJets.pdf")

#oblateness
selected_oblateness_hist = TH1F("selected_oblateness_hist","",20,-0.5,1)
selectedtup.Draw("pfos_injets_oblateness >> selected_oblateness_hist")
oblateness_hists.append(selected_oblateness_hist)

looseselected_oblateness_hist = TH1F("looseselected_oblateness_hist","",20,-0.5,1)
looseselectedtup.Draw("pfos_injets_oblateness >> looseselected_oblateness_hist")
oblateness_hists.append(looseselected_oblateness_hist)

tightselected_oblateness_hist = TH1F("tightselected_oblateness_hist","",20,-0.5,1)
tightselectedtup.Draw("pfos_injets_oblateness >> tightselected_oblateness_hist")
oblateness_hists.append(tightselected_oblateness_hist)

oblateness_hists[3].Scale(1./oblateness_hists[3].GetEntries())
oblateness_hists[4].Scale(1./oblateness_hists[4].GetEntries())
oblateness_hists[5].Scale(1./oblateness_hists[5].GetEntries())

oblateness_hists[3].SetTitle("")
oblateness_hists[3].GetXaxis().SetTitle("Oblateness")
oblateness_hists[3].Draw()
oblateness_hists[4].SetLineColor(kGreen+2)
oblateness_hists[4].Draw("same")
oblateness_hists[5].SetLineColor(kRed)
oblateness_hists[5].Draw("same")
oblateness_hists[3].Draw("same")

leg.SetX1NDC(0.59)
leg.SetX2NDC(0.89)
leg.Draw()
c.SaveAs("plots/oblateness_PFOsinJets.png")
c.SaveAs("plots/oblateness_PFOsinJets.pdf")


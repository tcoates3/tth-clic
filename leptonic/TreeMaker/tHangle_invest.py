# to investigate how the tH angle looks in signal and background

from util import *
from ROOT import gSystem
gSystem.Load("${LCIO}/lib/liblcio.so")
gSystem.Load("${LCIO}/lib/liblcioDict.so")
from ROOT import IOIMPL
from ROOT import TLorentzVector
import math

gStyle.SetOptStat(0)

fac_ins = IOIMPL.LCFactory.getInstance()

higgs_files = [100,102,103,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120]
higgs_hist = TH1F("higgs_hist","",20,0,3.2)
higgs_coshist = TH1F("higgs_coshist","",20,-1,1)
for higgs_file in higgs_files:

    reader = fac_ins.createLCReader()
    reader.open('/afs/cern.ch/work/s/sredford/tth_data/2441/000/FlavourTaggedDSTs/tth-ln4q-hbb_rec_2441_%i_flavourtagged.slcio' %higgs_file)

    event = reader.readNextEvent()
    while (event):
        mcp = event.getCollection('MCParticle')
    
        for i in xrange(20):
            mcpi = mcp.getElementAt(i)
            if mcpi.getPDG() == 25:
                higgs = mcpi
            elif mcpi.getPDG() == 21:
                higgs = mcpi
            elif mcpi.getPDG() == 23:
                higgs = mcpi

            if mcpi.getPDG() == 6:
                top = mcpi
            if mcpi.getPDG() == -6:
                topbar = mcpi

        higgs_4v = TLorentzVector(higgs.getMomentum())
        top_4v = TLorentzVector(top.getMomentum())
        topbar_4v = TLorentzVector(topbar.getMomentum())

        tHang1 = higgs_4v.Angle(top_4v.Vect())
        tHang2 = higgs_4v.Angle(topbar_4v.Vect())
        
        if (tHang1 < tHang2):
            higgs_hist.Fill(tHang1)
            higgs_coshist.Fill(math.cos(tHang1))
        else:
            higgs_hist.Fill(tHang2)
            higgs_coshist.Fill(math.cos(tHang2))

        event = reader.readNextEvent()

gluon_files = [100,102,103,104,105,106,109,110,111,112,113,114,115,116,118,119,121,122,123]
gluon_hist = TH1F("gluon_hist","",20,0,3.2)
gluon_coshist = TH1F("gluon_coshist","",20,-1,1)
for gluon_file in gluon_files:

    reader = fac_ins.createLCReader()
    reader.open('/afs/cern.ch/work/s/sredford/tth_data/ttbb-ln4q-all_rec_2426_%i_flavourtagged.slcio' %gluon_file)

    event = reader.readNextEvent()
    while (event):
        mcp = event.getCollection('MCParticle')
        
        for i in xrange(20):
            mcpi = mcp.getElementAt(i)
            if mcpi.getPDG() == 25:
                higgs = mcpi
            elif mcpi.getPDG() == 21:
                higgs = mcpi
            elif mcpi.getPDG() == 23:
                higgs = mcpi

            if mcpi.getPDG() == 6:
                top = mcpi
            if mcpi.getPDG() == -6:
                topbar = mcpi

        higgs_4v = TLorentzVector(higgs.getMomentum())
        top_4v = TLorentzVector(top.getMomentum())
        topbar_4v = TLorentzVector(topbar.getMomentum())

        tHang1 = higgs_4v.Angle(top_4v.Vect())
        tHang2 = higgs_4v.Angle(topbar_4v.Vect())

        if (tHang1 < tHang2):
            gluon_hist.Fill(tHang1)
            gluon_coshist.Fill(math.cos(tHang1))
        else:
            gluon_hist.Fill(tHang2)
            gluon_coshist.Fill(math.cos(tHang2))

        event = reader.readNextEvent()

z_files = [100,101,102,103,104,105,106,107,108,110,111,112,113,114,115,116,117,118,119]
z_hist = TH1F("z_hist","",20,0,3.2)
z_coshist = TH1F("z_coshist","",20,-1,1)
for z_file in z_files:
    reader = fac_ins.createLCReader()
    reader.open('/afs/cern.ch/work/s/sredford/tth_data/ttz-2l2nbb-all_rec_2447_%i_flavourtagged.slcio' %z_file)
    event = reader.readNextEvent()

    while (event):
        mcp = event.getCollection('MCParticle')
    
        for i in xrange(20):
            mcpi = mcp.getElementAt(i)
            if mcpi.getPDG() == 25:
                higgs = mcpi
            elif mcpi.getPDG() == 21:
                higgs = mcpi
            elif mcpi.getPDG() == 23:
                higgs = mcpi

            if mcpi.getPDG() == 6:
                top = mcpi
            if mcpi.getPDG() == -6:
                topbar = mcpi

        higgs_4v = TLorentzVector(higgs.getMomentum())
        top_4v = TLorentzVector(top.getMomentum())
        topbar_4v = TLorentzVector(topbar.getMomentum())

        tHang1 = higgs_4v.Angle(top_4v.Vect())
        tHang2 = higgs_4v.Angle(topbar_4v.Vect())
        
        if (tHang1 < tHang2):
            z_hist.Fill(tHang1)
            z_coshist.Fill(math.cos(tHang1))
        else:
            z_hist.Fill(tHang2)
            z_coshist.Fill(math.cos(tHang2))
            
        event = reader.readNextEvent()

gluon_hist.SetMinimum(0)
gluon_coshist.SetMinimum(0)

for hist in [higgs_hist,gluon_hist,z_hist,higgs_coshist,gluon_coshist,z_coshist]:
    hist.SetLineWidth(2)

higgs_hist.SetLineColor(kBlue)
gluon_hist.SetLineColor(kRed)
z_hist.SetLineColor(kOrange)
higgs_coshist.SetLineColor(kBlue)
gluon_coshist.SetLineColor(kRed)
z_coshist.SetLineColor(kOrange)

higgs_hist.Scale(1./higgs_hist.GetEntries())
gluon_hist.Scale(1./gluon_hist.GetEntries())
z_hist.Scale(1./z_hist.GetEntries())
higgs_coshist.Scale(1./higgs_hist.GetEntries())
gluon_coshist.Scale(1./gluon_hist.GetEntries())
z_coshist.Scale(1./z_hist.GetEntries())

gluon_hist.GetXaxis().SetTitle("Min top-boson angle (rad)")
gluon_hist.GetXaxis().SetTitleSize(0.05)
gluon_hist.GetXaxis().SetTitleOffset(0.9)
gluon_hist.GetXaxis().SetLabelSize(0.05)
gluon_hist.GetYaxis().SetLabelSize(0.05)

gluon_coshist.GetXaxis().SetTitle("cos( min top-boson angle )")
gluon_coshist.GetXaxis().SetTitleSize(0.05)
gluon_coshist.GetXaxis().SetTitleOffset(0.9)
gluon_coshist.GetXaxis().SetLabelSize(0.05)
gluon_coshist.GetYaxis().SetLabelSize(0.05)
gluon_coshist.GetXaxis().SetNdivisions(505)

leg = TLegend(0.15,0.7,0.6,0.89)
leg.SetFillColor(kWhite)
leg.SetTextFont(22)
leg.SetLineColor(kWhite)
leg.AddEntry(higgs_hist,"Higgs",'l')
leg.AddEntry(gluon_hist,"Gluon",'l')
leg.AddEntry(z_hist,"Z^{0}",'l')

c = TCanvas("c","",0,0,800,400)
c.Divide(2,1)
c.cd(1)
gluon_hist.Draw()
higgs_hist.Draw("same")
z_hist.Draw("same")
c.cd(2)
gluon_coshist.Draw()
higgs_coshist.Draw("same")
z_coshist.Draw("same")
leg.Draw()
c.Update()

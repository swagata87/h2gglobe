#!/bin/env python
from ROOT import TFile, TF1, TCanvas, kBlack, kBlue, kRed, kGreen, kCyan, kMagenta, kWhite, TGraph, TF2, TF12, kFullCircle, gROOT, TLegend, TPaveText
from sys import argv

def eff(fin,name,var):
    num=fin.Get("my_nvtx_rv_cat0_tot")
    den=fin.Get("my_nvtx_cat0_tot")
    eff = num.Clone(name)
    eff.Divide(num, den, 1., 1., "B")
    return eff

gROOT.LoadMacro("../Macros/rootglobestyle.C")
from ROOT import setTDRStyle, gStyle
setTDRStyle()
gStyle.SetOptTitle(0)
file=argv[1]
fin=TFile.Open(file)
eff_vtx = eff(fin,"eff2","nvtx")
print eff_vtx
eff_vtx.SetMarkerStyle(21);
eff_vtx.SetLineColor(kBlack);
pt = TPaveText(51.37856,0.9114967,190.8217,0.9923778,"br")
pt.SetFillColor(0)
pt.SetLineColor(0)
pt.SetTextAlign(13)
pt.SetTextFont(42)
txt=argv[2]
pt.AddText(txt)
pt.Draw("same")

c0 = TCanvas("vtxEffPt","vtxEffPt")
c3 = TCanvas("vtxEffNvtx","vtxEffNvtx")
eff3 = eff_vtx.Clone()
eff3.SetTitle("True vertex eff.;N_{vtx};Fraction of events")
eff3.Draw("e1")
leg3=TLegend(0.434564,0.784669,0.711409,0.850871)
leg3.SetShadowColor(kWhite), leg3.SetLineColor(kWhite), leg3.SetFillColor(kWhite), leg3.SetTextFont(60)
leg3.AddEntry(eff3,"| z_{reco} - z_{true} | < 10 mm","lpf")
#leg3.Draw("same")
pt.Draw("same")
for c in c0, c3:
    for fmt in "png", "pdf":
        c.SaveAs( "%s.%s" % (c.GetName(), fmt) )
        

#!/usr/bin/env python
#Finn Rebassoo, LLNL 12-01-2016
from os import listdir
from os.path import isfile, join
from ROOT import *
#import ROOT
import math as m
import sys

def calculate_t(p1,p3):
    p_comb=p1-p3
    test=p_comb*p_comb
    return test

#mypath="Ntupler"
#ListOfFiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
chain = TChain('demo/SlimmedNtuple')

#Old files:
#chain.Add("Ntupler"+"/SlimmedNtuple.root")
chain.Add("SlimmedNtuple.root")
i=0

#for file in ListOfFiles:
#    i=i+1
#    print file
#    chain.Add("Ntupler"+"/SlimmedNtuple.root")


files = [f for f in listdir('.') if isfile(f)]

#h_proton_rapidity_vs_mass=TH2F("h_proton_rapidity_vs_mass","",300,0,3000,50,-1.5,1.5)

h_xi_45=TH1F("h_xi_45",";#xi_{45};",100,0,1)
h_xi_56=TH1F("h_xi_56",";#xi_{56};",100,0,1)

h_xi_45_reco_arm45=TH1F("h_xi_45_reco_arm45",";#xi_{45};",100,0,1)
h_xi_56_reco_arm56=TH1F("h_xi_56_reco_arm56",";#xi_{56};",100,0,1)
h_xi_45_reco_pot_3=TH1F("h_xi_45_reco_pot_3",";#xi_{45};",100,0,1)
h_xi_45_reco_pot_23=TH1F("h_xi_45_reco_pot_23",";#xi_{45};",100,0,1)
h_xi_56_reco_pot_103=TH1F("h_xi_56_reco_pot_103",";#xi_{56};",100,0,1)
h_xi_56_reco_pot_123=TH1F("h_xi_56_reco_pot_123",";#xi_{56};",100,0,1)

h_xi_45_reco_arm45_eff=TH1F("h_xi_45_reco_arm45_eff",";#xi_{45};",100,0,1)
h_xi_56_reco_arm56_eff=TH1F("h_xi_56_reco_arm56_eff",";#xi_{56};",100,0,1)
h_xi_45_reco_pot_3_eff=TH1F("h_xi_45_reco_pot_3_eff",";#xi_{45};",100,0,1)
h_xi_45_reco_pot_23_eff=TH1F("h_xi_45_reco_pot_23_eff",";#xi_{45};",100,0,1)
h_xi_56_reco_pot_103_eff=TH1F("h_xi_56_reco_pot_103_eff",";#xi_{56};",100,0,1)
h_xi_56_reco_pot_123_eff=TH1F("h_xi_56_reco_pot_123_eff",";#xi_{56};",100,0,1)

h_num_pixels_45=TH1F("h_num_pixels_45",";Number pixels;",10,-0.5,9.5)
h_num_pixels_56=TH1F("h_num_pixels_56",";Number pixels;",10,-0.5,9.5)

h_mass_efficiency=TH1F("h_mass_efficiency",";Mass (GeV);",24,0,2400)
h_mass_efficiency_doubleArm=TH1F("h_mass_efficiency_doubleArm",";Mass (GeV);",24,0,2400)
h_mass_efficiency_doubleArm_pix=TH1F("h_mass_efficiency_doubleArm_pix",";Mass (GeV);",24,0,2400)
h_mass_efficiency_doubleArm_eff=TH1F("h_mass_efficiency_doubleArm_eff",";Mass (GeV);",24,0,2400)
h_mass_efficiency_doubleArm_pix_eff=TH1F("h_mass_efficiency_doubleArm_pix_eff",";Mass (GeV);",24,0,2400)

h_mass_reco_vs_sim_doubleArm=TH2F("h_mass_reco_vs_sim_doubleArm",";Gen Mass (GeV);Reco Mass (GeV)",48,0,2400,48,0,2400)
h_mass_reco_vs_sim_doubleArm_pix=TH2F("h_mass_reco_vs_sim_doubleArm_pix",";Gen Mass (GeV);Reco Mass (GeV)",48,0,2400,48,0,2400)

h_xi_diff_vs_xi_gen_45=TH2F("h_xi_diff_vs_xi_gen_45",";#xi_{gen};#xi_{reco}-#xi_{gen}",200,0,1,100,-0.2,0.2)
h_xi_diff_vs_xi_gen_56=TH2F("h_xi_diff_vs_xi_gen_56",";#xi_{gen};#xi_{reco}-#xi_{gen}",200,0,1,100,-0.2,0.2)

h_xi_diff_vs_xi_gen_45_p=TH2F("h_xi_diff_vs_xi_gen_45_p",";#xi_{gen};#xi_{reco}-#xi_{gen}",200,0,1,100,-0.2,0.2)
h_xi_diff_vs_xi_gen_56_p=TH2F("h_xi_diff_vs_xi_gen_56_p",";#xi_{gen};#xi_{reco}-#xi_{gen}",200,0,1,100,-0.2,0.2)

h_xi_diff_divided_xi_gen_45=TH1F("h_xi_diff_divided_xi_gen_45","Multi-RP S45;#xi_{reco}-#xi_{gen}/#xi_{gen};",100,-1,1)
h_xi_diff_divided_xi_gen_56=TH1F("h_xi_diff_divided_xi_gen_56","Multi-RP S56;#xi_{reco}-#xi_{gen}/#xi_{gen};",100,-1,1)
h_xi_diff_divided_xi_gen_45_p=TH1F("h_xi_diff_divided_xi_gen_45_p","Pixel S45;#xi_{reco}-#xi_{gen}/#xi_{gen};",100,-1,1)
h_xi_diff_divided_xi_gen_56_p=TH1F("h_xi_diff_divided_xi_gen_56_p","Pixel S56;#xi_{reco}-#xi_{gen}/#xi_{gen};",100,-1,1)


#h_proton_eff_vs_mass

print chain.GetEntries()
for e in chain:
    #print "Run: {0}, Event: {1}".format(e.run,e.event)
    #print "I get into loop"
    p_pos_four_vector=TLorentzVector()
    p_neg_four_vector=TLorentzVector()
    ip=0
    xi_neg=0.
    xi_pos=0.
    theta_y_1=0
    t1=-999
    t2=-999
    for p in e.gen_proton_pz:
        #if p < 0 and abs(p) > 1000:
        if p < 0:
            #p_neg goes to sector 56
            p_neg_four_vector=TLorentzVector(e.gen_proton_px[ip],e.gen_proton_py[ip],e.gen_proton_pz[ip],e.gen_proton_energy[ip])
            p_neg_incom=TLorentzVector(0,0,-6500.,m.sqrt(0.93827231*0.93827231+6500*6500))
            t1=calculate_t(p_neg_incom,p_neg_four_vector)
            pt=m.sqrt(e.gen_proton_px[ip]*e.gen_proton_px[ip]+e.gen_proton_py[ip]*e.gen_proton_py[ip])
            #print t1,m.sqrt(e.hepmc_px[ip]*e.hepmc_px[ip]+e.hepmc_py[ip]*e.hepmc_py[ip])
            #print e.hepmc_px[ip],e.hepmc_py[ip],e.hepmc_pz[ip],e.hepmc_energy[ip]
            xi_neg=e.gen_proton_pz[ip]/6500.
            h_xi_56.Fill(1.-abs(xi_neg))
            theta_y_1=abs(m.atan(e.gen_proton_py[ip]/e.gen_proton_pz[ip]))
            #print 1.-abs(xi_neg)
        #if p > 0 and abs(p) > 1000:
        if p > 0:
            #Positive Z point west at CMS, X points south towards middle of LHC ring
            #p_pos goes to sector 45
            p_pos_four_vector=TLorentzVector(e.gen_proton_px[ip],e.gen_proton_py[ip],e.gen_proton_pz[ip],e.gen_proton_energy[ip])
            p_pos_incom=TLorentzVector(0,0,6500.,m.sqrt(0.93827231*0.93827231+6500*6500))
            t2=calculate_t(p_pos_four_vector,p_pos_incom)
            #print t2,e.gen_proton_px[ip],e.gen_proton_py[ip],e.gen_proton_pz[ip],e.gen_proton_energy[ip]
            xi_pos=e.gen_proton_pz[ip]/6500.
            h_xi_45.Fill(1.-abs(xi_pos))
            theta_y_2=abs(m.atan(e.gen_proton_py[ip]/e.gen_proton_pz[ip]))
        ip=ip+1
    #print "t1: {0}, t2: {1}".format(t1,t2)
    #if ip > 2:
    #    print "More than 3 protons from gen_proton"
    #print (1.-abs(xi_neg)),(1.-abs(xi_pos))
    xi_56=1.-abs(xi_neg)
    xi_45=1.-abs(xi_pos)
    #if xi_45 < 0.2 and xi_45 < 0.2 and xi_45 > 0.01 and xi_45 > 0.01:
    #h_t_1.Fill(abs(t1))
    #h_t_2.Fill(abs(t2))
    #h_xi_56_vs_t1.Fill(abs(t1),xi_56)
    #h_xi_45_vs_t2.Fill(abs(t2),xi_45)
    #h_t_1_vs_pt.Fill(pt,abs(t1))
    #h_t_2_vs_t_1.Fill(abs(t2),abs(t1))
    #h_xi_45_vs_1.Fill(1.-abs(xi_neg),1.-abs(xi_pos))
    proton_combined=p_pos_four_vector+p_neg_four_vector
    mass=proton_combined.M()
    mass=m.sqrt(13000*13000*(1-abs(xi_neg))*(1-abs(xi_pos)))
    rapidity=proton_combined.Rapidity()
    rapidity=0.5*m.log((1-abs(xi_pos))/(1-abs(xi_neg)))

    #h_theta_y_1.Fill(theta_y_1)
    #h_theta_y_2.Fill(theta_y_2)


    passRapidity=True
    #both_arms=False
    ii=0
    pot_3=False
    pot_23=False
    pot_103=False
    pot_123=False
    arm45=False
    arm56=False
    count_45=0
    count_56=0
    xi_pot_23=-999
    xi_pot_123=-999
    for t in e.proton_rpid:
        if e.proton_ismultirp_[ii]:
            if e.proton_arm[ii]==0:
                xi_arm_45=e.proton_xi[ii]
                #if (abs(xi_arm_45 - xi_45)/xi_45) < 0.1:
                h_xi_diff_vs_xi_gen_45.Fill(xi_45,xi_arm_45-xi_45)
                h_xi_diff_divided_xi_gen_45.Fill((xi_arm_45-xi_45)/xi_45)
                arm45=True
            if e.proton_arm[ii]==1:
                xi_arm_56=e.proton_xi[ii]
                h_xi_diff_vs_xi_gen_56.Fill(xi_56,xi_arm_56-xi_56)
                h_xi_diff_divided_xi_gen_56.Fill((xi_arm_56-xi_56)/xi_56)
                #if (abs(xi_arm_56 - xi_56)/xi_56) < 0.1:
                arm56=True
        if t == 3: 
            pot_3=True
        if t ==23:
            pot_23=True
            count_45=count_45+1
            xi_pot_23=e.proton_xi[ii]
            h_xi_diff_vs_xi_gen_45_p.Fill(xi_45,e.proton_xi[ii]-xi_45)
            h_xi_diff_divided_xi_gen_45_p.Fill((e.proton_xi[ii]-xi_45)/xi_45)
        if t == 103: 
            pot_103=True
        if t ==123:
            pot_123=True
            count_56=count_56+1
            xi_pot_123=e.proton_xi[ii]
            h_xi_diff_vs_xi_gen_56_p.Fill(xi_56,e.proton_xi[ii]-xi_56)
            h_xi_diff_divided_xi_gen_56_p.Fill((e.proton_xi[ii]-xi_56)/xi_56)
        ii=ii+1
        
    h_num_pixels_45.Fill(count_45)    
    h_num_pixels_56.Fill(count_56)    
    if arm45:
        h_xi_45_reco_arm45.Fill(xi_45)
    if arm56:
        h_xi_56_reco_arm56.Fill(xi_56)
    if pot_3:
        h_xi_45_reco_pot_3.Fill(xi_45)
    if pot_23:
        h_xi_45_reco_pot_23.Fill(xi_45)
    if pot_103:
        h_xi_56_reco_pot_103.Fill(xi_56)
    if pot_123:
        h_xi_56_reco_pot_123.Fill(xi_56)

    h_mass_efficiency.Fill(m.sqrt(13000*13000*xi_45*xi_56))
    if pot_23 and pot_123:
        h_mass_efficiency_doubleArm_pix.Fill(m.sqrt(13000*13000*xi_45*xi_56))
        h_mass_reco_vs_sim_doubleArm_pix.Fill(mass,m.sqrt(13000*13000*xi_pot_23*xi_pot_123))
    if arm45 and arm56:
        h_mass_efficiency_doubleArm.Fill(m.sqrt(13000*13000*xi_45*xi_56))
        h_mass_reco_vs_sim_doubleArm.Fill(mass,m.sqrt(13000*13000*xi_arm_45*xi_arm_56))

h_xi_45.Sumw2()
h_xi_56.Sumw2()
h_xi_45_reco_arm45.Sumw2()
h_xi_56_reco_arm56.Sumw2()
h_xi_45_reco_pot_3.Sumw2()
h_xi_45_reco_pot_23.Sumw2()
h_xi_56_reco_pot_103.Sumw2()
h_xi_56_reco_pot_123.Sumw2()

h_xi_45_reco_arm45_eff.Divide(h_xi_45_reco_arm45,h_xi_45,1,1,"b")
h_xi_56_reco_arm56_eff.Divide(h_xi_56_reco_arm56,h_xi_56,1,1,"b")
h_xi_45_reco_pot_3_eff.Divide(h_xi_45_reco_pot_3,h_xi_45,1,1,"b")
h_xi_45_reco_pot_23_eff.Divide(h_xi_45_reco_pot_23,h_xi_45,1,1,"b")
h_xi_56_reco_pot_103_eff.Divide(h_xi_56_reco_pot_103,h_xi_56,1,1,"b")
h_xi_56_reco_pot_123_eff.Divide(h_xi_56_reco_pot_123,h_xi_56,1,1,"b")

h_mass_efficiency.Sumw2()
h_mass_efficiency_doubleArm.Sumw2()
h_mass_efficiency_doubleArm_eff.Divide(h_mass_efficiency_doubleArm,h_mass_efficiency,1,1,"b")
h_mass_efficiency_doubleArm_pix.Sumw2()
h_mass_efficiency_doubleArm_pix_eff.Divide(h_mass_efficiency_doubleArm_pix,h_mass_efficiency,1,1,"b")

c0=TCanvas("c0","",1000,800)
gPad.SetLeftMargin(0.15)
gPad.SetRightMargin(0.15)
c0.cd()
h_mass_efficiency_doubleArm_pix_eff.SetStats(0)
h_mass_efficiency_doubleArm_pix_eff.Draw()
h_mass_efficiency_doubleArm_eff.SetLineColor(2)
h_mass_efficiency_doubleArm_eff.Draw("same")
c0.Print("Mass_efficiency.pdf")

c1=TCanvas("c1","",800,800)
gPad.SetLeftMargin(0.15)
gPad.SetRightMargin(0.15)
c1.cd()


h_xi_45_reco_arm45_eff.SetStats(0)
h_xi_45_reco_arm45_eff.SetLineColor(1)
h_xi_45_reco_arm45_eff.Draw()
h_xi_45_reco_pot_3_eff.Draw("same")
h_xi_45_reco_pot_23_eff.SetLineColor(2)
h_xi_45_reco_pot_23_eff.Draw("same")   

#h_xi_pos_acc_n.Sumw2()
#h_xi_pos_acc_d.Sumw2()
#h_xi_neg_acc_n.Sumw2()
#h_xi_neg_acc_d.Sumw2()

#h_xi_pos_acc_n.Divide(h_xi_pos_acc_n,h_xi_pos_acc_d,1,1,"b")
#h_xi_neg_acc_n.Divide(h_xi_neg_acc_n,h_xi_neg_acc_d,1,1,"b")

#h_proton_rapidity_vs_mass=TH2F("h_proton_rapidity_vs_mass","Double-arm acceptance;Generated diproton mass [GeV];Generated diproton rapidity",30,0,3000,30,-1.5,1.5)
#h_proton_rapidity_vs_mass.Divide(h_num,h_denom)
#h_proton_rapidity_vs_mass.SetStats(0)

#h_proton_eff_vs_mass.SetStats(0)
#h_proton_eff_vs_mass_pixels.SetStats(0)

#sys.exit()


h_xi_45_reco_arm45_eff.SetStats(0)
h_xi_45_reco_arm45_eff.GetXaxis().SetRangeUser(0,0.25)
h_xi_45_reco_arm45_eff.SetLineColor(1)
h_xi_45_reco_arm45_eff.SetMarkerStyle(25)
h_xi_45_reco_arm45_eff.SetMarkerColor(1)
h_xi_45_reco_arm45_eff.Draw()
h_xi_45_reco_pot_3_eff.SetMarkerStyle(24)
h_xi_45_reco_pot_3_eff.SetMarkerColor(4)
h_xi_45_reco_pot_3_eff.Draw("same")
h_xi_45_reco_pot_23_eff.SetLineColor(2)
h_xi_45_reco_pot_23_eff.SetMarkerStyle(26)
h_xi_45_reco_pot_23_eff.SetMarkerColor(2)
h_xi_45_reco_pot_23_eff.Draw("same")   
h_xi_45_reco_arm45_eff.Draw("same")
leg=TLegend(0.65,0.56,0.82,0.886)
leg.AddEntry(h_xi_45_reco_arm45_eff,"Combined Proton Reco","p")
leg.AddEntry(h_xi_45_reco_pot_3_eff,"Strips Reco","p")
leg.AddEntry(h_xi_45_reco_pot_23_eff,"Pixels Reco","p")
leg.Draw("same")
c1.Print("xi_eff_45.pdf")

h_xi_56_reco_arm56_eff.SetStats(0)
h_xi_56_reco_arm56_eff.GetXaxis().SetRangeUser(0,0.25)
h_xi_56_reco_arm56_eff.SetLineColor(1)
h_xi_56_reco_arm56_eff.SetMarkerStyle(25)
h_xi_56_reco_arm56_eff.SetMarkerColor(1)
h_xi_56_reco_arm56_eff.Draw()
h_xi_56_reco_pot_103_eff.SetMarkerStyle(24)
h_xi_56_reco_pot_103_eff.SetMarkerColor(4)
h_xi_56_reco_pot_103_eff.Draw("same")
h_xi_56_reco_pot_123_eff.SetLineColor(2)
h_xi_56_reco_pot_123_eff.SetMarkerStyle(26)
h_xi_56_reco_pot_123_eff.SetMarkerColor(2)
h_xi_56_reco_pot_123_eff.Draw("same")   
h_xi_56_reco_arm56_eff.Draw("same")
leg.Draw("Same")
c1.Print("xi_eff_56.pdf")

entries=h_num_pixels_45.GetEntries()
h_num_pixels_45.Scale(1/entries)
entries_56=h_num_pixels_56.GetEntries()
h_num_pixels_56.Scale(1/entries_56)
h_num_pixels_45.SetStats(0)
h_num_pixels_45.GetXaxis().SetRangeUser(0,4)
h_num_pixels_45.Draw()
h_num_pixels_56.SetLineColor(2)
h_num_pixels_56.Draw("same")

c1.Print("NumberPixels.pdf")

h_xi_diff_vs_xi_gen_45.ProfileX("h_profile_45",1,-1,"s")
h_xi_diff_vs_xi_gen_56.ProfileX("h_profile_56",1,-1,"s")
h_xi_diff_vs_xi_gen_45_p.ProfileX("h_profile_45_p",1,-1,"s")
h_xi_diff_vs_xi_gen_56_p.ProfileX("h_profile_56_p",1,-1,"s")
nbins=h_profile_45.GetNbinsX()
h_RMS_45=TH1F("h_RMS_45",";#xi_{sim};RMS(#xi_{reco}-#xi_{sim})",200,0,1)
h_RMS_56=h_RMS_45.Clone()
h_RMS_45_p=h_RMS_45.Clone()
h_RMS_56_p=h_RMS_45.Clone()
for i in range(1,nbins+1):
    RMS_45=h_profile_45.GetBinError(i)
    RMS_56=h_profile_56.GetBinError(i)
    RMS_45_p=h_profile_45_p.GetBinError(i)
    RMS_56_p=h_profile_56_p.GetBinError(i)
    h_RMS_45.SetBinContent(i,RMS_45)
    h_RMS_56.SetBinContent(i,RMS_56)
    h_RMS_45_p.SetBinContent(i,RMS_45_p)
    h_RMS_56_p.SetBinContent(i,RMS_56_p)
    

h_RMS_45.GetXaxis().SetRangeUser(0,0.2)
h_RMS_56.GetXaxis().SetRangeUser(0,0.2)
h_RMS_45.GetYaxis().SetTitleOffset(2.0)
h_RMS_56.GetYaxis().SetTitleOffset(2.0)
h_RMS_45.SetStats(0)
h_RMS_56.SetStats(0)
h_RMS_45.SetTitle("45")
h_RMS_56.SetTitle("56")
h_profile_45.GetXaxis().SetRangeUser(0,0.2)
h_profile_45.GetYaxis().SetTitle("#xi_{reco}-#xi_{sim}")
h_profile_56.GetXaxis().SetRangeUser(0,0.2)
h_profile_56.GetYaxis().SetTitle("#xi_{reco}-#xi_{sim}")

h_RMS_45_p.GetXaxis().SetRangeUser(0,0.2)
h_RMS_56_p.GetXaxis().SetRangeUser(0,0.2)
h_RMS_45_p.GetYaxis().SetTitleOffset(2.0)
h_RMS_56_p.GetYaxis().SetTitleOffset(2.0)
h_RMS_45_p.SetStats(0)
h_RMS_56_p.SetStats(0)
h_RMS_45_p.SetTitle("45")
h_RMS_56_p.SetTitle("56")
h_profile_45_p.GetXaxis().SetRangeUser(0,0.2)
h_profile_45_p.GetYaxis().SetTitle("#xi_{reco}-#xi_{sim}")
h_profile_56_p.GetXaxis().SetRangeUser(0,0.2)
h_profile_56_p.GetYaxis().SetTitle("#xi_{reco}-#xi_{sim}")


h_RMS_45.Draw("hist")
c1.Print("RMSArm45.pdf")
h_RMS_56.Draw("hist")
c1.Print("RMSArm56.pdf")
h_profile_45.Draw("hist")
#c1.Print("profileArm45.pdf")
h_profile_56.Draw("hist")
#c1.Print("profileArm56.pdf")

h_RMS_45_p.Draw("hist")
c1.Print("RMSArm45_p.pdf")
h_RMS_56_p.Draw("hist")
c1.Print("RMSArm56_p.pdf")
h_profile_45_p.Draw("hist")
#c1.Print("profileArm45.pdf")
h_profile_56_p.Draw("hist")
#c1.Print("profileArm56.pdf")

h_xi_diff_divided_xi_gen_45.Draw()
c1.Print("xi_residual_45_multiRP.pdf")

h_xi_diff_divided_xi_gen_56.Draw()
c1.Print("xi_residual_56_multiRP.pdf")

h_xi_diff_divided_xi_gen_45_p.Draw()
c1.Print("xi_residual_45_pixel.pdf")

h_xi_diff_divided_xi_gen_56_p.Draw()
c1.Print("xi_residual_56_pixel.pdf")

#------------------------------------------ 
#   Hcal PFG script
#   Description :
#       Compares Pedestals and Produces Stats
#   Author : Viktor Khristenko
#            (victor-khristenko@uiowa.edu)
#------------------------------------------ 
# 
# Usage 
# $ python ped_compare.py [ped file1] [ped file2] [output root file]  
# 
#   file1 : pedestal table file 1 with full path  
#   file2 : pedestal table file 2 with full path 
#   output root file : output root file with full path 
#
#------------------------------------------ 
#
# The old pedestal tables are in 
#    https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPedestalsTags2011
#
#------------------------------------------ 
#
# This script compares only the channels that exist in both files.
# If a channel exists in one file while not in the other file, 
# this channel is excluded in comparison.
#
# The output root file will contatin (for each subdetector)
#   - 1D(as a function of ADC) and 2D(ieta vs iphi) ratio of mean for each cap (file1/file2)
#   - 1D(as a function of ADC) and 2D(ieta vs iphi) ratio of sigma for each cap (file2/file2)
#   - summary of ratios
#
#------------------------------------------ 

import sys, copy
import ROOT as rr

if len(sys.argv)<4: sys.exit('You need more argements! \nThis is usage : \n\n    python ped_compare.py [ped file1] [ped file2] [output root file] \n')


#------------------------------------------ 
#   Define a Channel Class
#------------------------------------------ 
class HcalChannel:
    def __init__(self, ieta, iphi, depth, name):
        self.__ieta         = ieta
        self.__iphi         = iphi
        self.__depth        = depth
        self.__name         = name
        self.__requested    = True
        self.enabled        = False
        self.means          = [0 for i in range(4)]
        self.rmss           = [0 for i in range(4)]
        self.rawId          = 0

    def __str__(self):
        str1 = '<HcalChannel: ieta=%d iphi=%d depth=%d sub=%s\n' % (
            self.__ieta, self.__iphi, self.__depth, self.__name)
        str2 = 'mean0=%f mean1=%f mean2=%f mean3=%f\n' % (
            self.__vMeans[0], self.__vMeans[1], self.__vMeans[2], 
            self.__vMeans[3])
        str3 = 'rms0=%f rms1=%f rms2=%f rms3=%f>\n' % (
            self.__vRMSs[0], self.__vRMSs[1], self.__vRMSs[2], self.__vRMSs[3])
        return str1+str2+str3

    def iamenabled(self):
        str1 = '<HcalChannel: ieta=%d iphi=%d depth=%d sub=%s>\n' % (
            self.__ieta, self.__iphi, self.__depth, self.__name)
        str2 = ' is enabled already!!!'
        print str1+str2

#------------------------------------------ 
#   Define a Subsystem Class
#------------------------------------------ 
class HcalSubsystem:
    def __init__(self, name, ietamin, ietamax,
            iphimin, iphimax, iphistep, depthmin, depthmax, verbosity=0):
        self.__name = name
        self.__ietamin = ietamin
        self.__ietamax = ietamax
        self.__ietastep = 1
        self.__ietanum = 2*((ietamax - ietamin) + 1)
        self.__iphimin = iphimin
        self.__iphimax = iphimax
        self.__iphistep = iphistep
        self.__iphinum = (iphimax - iphimin)/iphistep + 1
        self.__depthmin = depthmin
        self.__depthmax = depthmax
        self.__depthstep = 1
        self.__depthnum = depthmax - depthmin + 1
        self.__verbosity = verbosity
        self.requestChs()

    def get_iphinum(self):
        return self.__iphinum
    def get_ietanum(self):
        return self.__ietanum
    def get_depthnum(self):
        return self.__depthnum

    def iphi(self, iiphi):
        return self.__iphimin + self.__iphistep*iiphi

    def iiphi(self, iphi):
        return (iphi-self.__iphimin)/self.__iphistep
    
    def ieta(self, iieta):
        if iieta<self.__ietanum/2:
            return -(self.__ietamin + self.__ietastep*iieta)
        else:
            return self.__ietamin + self.__ietastep*(iieta-self.__ietanum/2)

    def iieta(self, ieta):
        if ieta<0:
            return (abs(ieta)-self.__ietamin)/self.__ietastep
        else:
            return (ieta-self.__ietamin)/self.__ietastep + self.__ietanum/2

    def depth(self, id):
        return self.__depthmin + id*self.__depthstep

    def idepth(self, depth):
        return (depth-self.__depthmin)/self.__depthstep

    def requestChs(self): 
        self.chs = [[[
            HcalChannel(self.ieta(iieta), self.iphi(iiphi), 
                self.depth(id), self.__name) 
            for id in range(self.__depthnum)]
            for iieta in range(self.__ietanum)]
            for iiphi in range(self.__iphinum)]
        self.__numChsRequested = self.__depthnum*self.__iphinum*self.__ietanum

    def addCh(self, ieta, iphi, depth, means, rmss, rawId):
        self.__debug("### Parsing Channel: %d %d %d for SubSystem %s\n" % 
                (ieta, iphi, depth, self.__name))

        iieta = self.iieta(ieta)
        iiphi = self.iiphi(iphi)
        idepth = self.idepth(depth)
        if self.chs[iiphi][iieta][idepth].enabled:
            self.chs[iiphi][iieta][idepth].iamenabled()
            return

        self.chs[iiphi][iieta][idepth].means = copy.deepcopy(means)
        self.chs[iiphi][iieta][idepth].rmss = copy.deepcopy(rmss)
        self.chs[iiphi][iieta][idepth].rawId = rawId
        self.chs[iiphi][iieta][idepth].enabled=True

    def __debug(self, s):
        if self.__verbosity>5:
            print s

#------------------------------------------ 
#   Define an input function to fill the Subsystems
#------------------------------------------ 
def parse(fileName, subs):
    print "### Reading in %s\n" % fileName

    # do inits
    f = open(fileName)
    lineCounter = 0

    #   parse
    for line in f:
        #Skip lines starting with #
        if line[0]=="#":
            continue

        v = line.split()

        #   Get SubSystem/Coordinates
        ieta = int(v[0])
        iphi = int(v[1])
        depth = int(v[2])
        sub = v[3]

        #   Throw away channels that are not needed
        if sub!="HB" and sub!="HE" and sub!="HO" and sub!="HF" and sub!="QIE11":
            continue

        #   Get Means
        m0 = float(v[4])
        m1 = float(v[5])
        m2 = float(v[6])
        m3 = float(v[7])
        means = [m0, m1, m2, m3]

        #   Get RMSs
        r0 = float(v[8])
        r1 = float(v[9])
        r2 = float(v[10])
        r3 = float(v[11])
        rmss = [r0, r1, r2, r3]

        #   Get the rawId
        rawId = int(v[12], 16)

        print "%d"% depth

        #   Add this channel to the proper subsytem
        subs[sub].addCh(ieta, iphi, depth, means, rmss, rawId)

#------------------------------------------ 
#   Define a StatsMaker Class
#------------------------------------------ 
class StatsMaker:
    def __init__(self, rootFile, name):
        self.__name = name

        print '### Initializing Stats for Subsystem %s\n' % self.__name 
        rootFile.cd()
        rootFile.mkdir(self.__name)
        rootFile.cd(self.__name)
        rr.gDirectory.mkdir("2D")
        self.__hm_diff      = rr.TH1D("hm_diff", "Difference of Means", 100, -1, 1)
        self.__hm0_diff     = rr.TH1D("hm0_diff", "Difference of Means Cap 0", 100, -1, 1)
        self.__hm1_diff     = rr.TH1D("hm1_diff", "Difference of Means Cap 1", 100, -1, 1)
        self.__hm2_diff     = rr.TH1D("hm2_diff", "Difference of Means Cap 2", 100, -1, 1)
        self.__hm3_diff     = rr.TH1D("hm3_diff", "Difference of Means Cap 3", 100, -1, 1)
        self.__hr0_diff     = rr.TH1D("hr0_diff", "Difference of RMSs Cap 0", 100, -1, 1)
        self.__hr1_diff     = rr.TH1D("hr1_diff", "Difference of RMSs Cap 1", 100, -1, 1)
        self.__hr2_diff     = rr.TH1D("hr2_diff", "Difference of RMSs Cap 2", 100, -1, 1)
        self.__hr3_diff     = rr.TH1D("hr3_diff", "Difference of RMSs Cap 3", 100, -1, 1)
        self.__hm0_sub1     = rr.TH1D("hm0_sub1", "Mean Cap 0 for run 269836", 100, 0, 16)
        self.__hm1_sub1     = rr.TH1D("hm1_sub1", "Mean Cap 1 for run 269836", 100, 0, 16)
        self.__hm2_sub1     = rr.TH1D("hm2_sub1", "Mean Cap 2 for run 269836", 100, 0, 16)
        self.__hm3_sub1     = rr.TH1D("hm3_sub1", "Mean Cap 3 for run 269836", 100, 0, 16)
        self.__hr0_sub1     = rr.TH1D("hr0_sub1", "RMS Cap 0 for run 269836", 100, 0, 2)
        self.__hr1_sub1     = rr.TH1D("hr1_sub1", "RMS Cap 1 for run 269836", 100, 0, 2)
        self.__hr2_sub1     = rr.TH1D("hr2_sub1", "RMS Cap 2 for run 269836", 100, 0, 2)
        self.__hr3_sub1     = rr.TH1D("hr3_sub1", "RMS Cap 3 for run 269836", 100, 0, 2)
        self.__hm0      = rr.TH1D("hm0", "Ratio of Means Cap 0", 100, 0, 2)
        self.__hm1      = rr.TH1D("hm1", "Ratio of Means Cap 1", 100, 0, 2)
        self.__hm2      = rr.TH1D("hm2", "Ratio of Means Cap 2", 100, 0, 2)
        self.__hm3      = rr.TH1D("hm3", "Ratio of Means Cap 3", 100, 0, 2)
        self.__hr0      = rr.TH1D("hr0", "Ratio RMSs Cap 0", 100, 0, 2)
        self.__hr1      = rr.TH1D("hr1", "Ratio RMSs Cap 1", 100, 0, 2)
        self.__hr2      = rr.TH1D("hr2", "Ratio RMSs Cap 2", 100, 0, 2)
        self.__hr3      = rr.TH1D("hr3", "Ratio RMSs Cap 2", 100, 0, 2)
        self.__hrawId   = rr.TH1D("hrawId", 
                "Ratio of RawIds(Base16) converted to base 10", 100, 0, 2)
        self.__hSummary = rr.TH2D("hSummary", 
                "Summary over all Caps(0-3) for RMS/Mean/rawId", 
                9, 0, 9, 2, 0, 2);
        for i in range(4):
            self.__hSummary.GetXaxis().SetBinLabel(i+1, "Means Cap%d" % i)
            self.__hSummary.GetXaxis().SetBinLabel(i+5, "RMSs Cap%d" % i)
        self.__hSummary.GetXaxis().SetBinLabel(9, "rawIds")
        self.__hSummary.GetYaxis().SetBinLabel(1, "Mean Ratio")
        self.__hSummary.GetYaxis().SetBinLabel(2, "RMS of the Ratio")

        #   Put it by hand
        if name=="HF":
            depthnum = 2
            depthmin = 1
        elif name=="HB" or name=="HE":
            depthnum = 3
            depthmin = 1
        elif name=="HO":
            depthnum = 1
            depthmin = 4
        if name=="QIE11":
            depthnum = 7 
            depthmin = 1

        rr.gDirectory.cd("2D")
        self.__vh2rm0_diff = [rr.TH2D("h2rm0_diff_%d" % (i+depthmin), 
            "Difference of Means Cap 0 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm1_diff = [rr.TH2D("h2rm1_diff_%d" % (i+depthmin), 
            "Difference of Means Cap 1 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm2_diff = [rr.TH2D("h2rm2_diff_%d" % (i+depthmin), 
            "Difference of Means Cap 2 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm3_diff = [rr.TH2D("h2rm3_diff_%d" % (i+depthmin), 
            "Difference of Means Cap 3 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr0_diff = [rr.TH2D("h2rr0_diff_%d" % (i+depthmin), 
            "Difference of RMSs Cap 0 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr1_diff = [rr.TH2D("h2rr1_diff_%d" % (i+depthmin), 
            "Difference of RMSs Cap 1 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr2_diff = [rr.TH2D("h2rr2_diff_%d" % (i+depthmin), 
            "Difference of RMSs Cap 2 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr3_diff = [rr.TH2D("h2rr3_diff_%d" % (i+depthmin), 
            "Difference of RMSs Cap 3 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm0 = [rr.TH2D("h2rm0_%d" % (i+depthmin), 
            "Ratio of Means Cap 0 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm1 = [rr.TH2D("h2rm1_%d" % (i+depthmin), 
            "Ratio of Means Cap 1 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm2 = [rr.TH2D("h2rm2_%d" % (i+depthmin), 
            "Ratio of Means Cap 2 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rm3 = [rr.TH2D("h2rm3_%d" % (i+depthmin), 
            "Ratio of Means Cap 3 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr0 = [rr.TH2D("h2rr0_%d" % (i+depthmin), 
            "Ratio of RMSs Cap 0 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr1 = [rr.TH2D("h2rr1_%d" % (i+depthmin), 
            "Ratio of RMSs Cap 1 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr2 = [rr.TH2D("h2rr2_%d" % (i+depthmin), 
            "Ratio of RMSs Cap 2 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]
        self.__vh2rr3 = [rr.TH2D("h2rr3_%d" % (i+depthmin), 
            "Ratio of RMSs Cap 3 D%d" % (i+depthmin),
            83, -41.5, 41.5, 72, 0.5, 72.5) for i in range(depthnum)]


    def makeStats(self, sub1, sub2):
        iphinum = sub1.get_iphinum()
        ietanum = sub1.get_ietanum()
        depthnum = sub1.get_depthnum()

        print "### Producing Stats for Subsystem %s\n" % self.__name
        for iiphi in range(iphinum):
            for iieta in range(ietanum):
                for id in range(depthnum):
                    if (not sub1.chs[iiphi][iieta][id].enabled or
                            not sub2.chs[iiphi][iieta][id].enabled):
                        continue

                    if abs(sub1.chs[iiphi][iieta][id].means[0]-sub2.chs[iiphi][iieta][id].means[0]) > 1.: 
                        print "ieta=%d iphi=%d depth=%d capid=0" % (sub1.ieta(iieta), sub1.iphi(iiphi), id)   
                    if abs(sub1.chs[iiphi][iieta][id].means[1]-sub2.chs[iiphi][iieta][id].means[1]) > 1.: 
                        print "ieta=%d iphi=%d depth=%d capid=1" % (sub1.ieta(iieta), sub1.iphi(iiphi), id)   
                    if abs(sub1.chs[iiphi][iieta][id].means[2]-sub2.chs[iiphi][iieta][id].means[2]) > 1.: 
                        print "ieta=%d iphi=%d depth=%d capid=2" % (sub1.ieta(iieta), sub1.iphi(iiphi), id)   
                    if abs(sub1.chs[iiphi][iieta][id].means[3]-sub2.chs[iiphi][iieta][id].means[3]) > 1.: 
                        print "ieta=%d iphi=%d depth=%d capid=3" % (sub1.ieta(iieta), sub1.iphi(iiphi), id)   
                    
                    rm0 = (sub1.chs[iiphi][iieta][id].means[0] / 
                            sub2.chs[iiphi][iieta][id].means[0])
                    rm1 = (sub1.chs[iiphi][iieta][id].means[1] / 
                            sub2.chs[iiphi][iieta][id].means[1])
                    rm2 = (sub1.chs[iiphi][iieta][id].means[2] / 
                            sub2.chs[iiphi][iieta][id].means[2])
                    rm3 = (sub1.chs[iiphi][iieta][id].means[3] / 
                            sub2.chs[iiphi][iieta][id].means[3])
                    rr0 = (sub1.chs[iiphi][iieta][id].rmss[0] / 
                            sub2.chs[iiphi][iieta][id].rmss[0])
                    rr1 = (sub1.chs[iiphi][iieta][id].rmss[1] / 
                            sub2.chs[iiphi][iieta][id].rmss[1])
                    rr2 = (sub1.chs[iiphi][iieta][id].rmss[2] / 
                            sub2.chs[iiphi][iieta][id].rmss[2])
                    rr3 = (sub1.chs[iiphi][iieta][id].rmss[3] / 
                            sub2.chs[iiphi][iieta][id].rmss[3])
                    rraw = (sub1.chs[iiphi][iieta][id].rawId /
                            sub2.chs[iiphi][iieta][id].rawId)
                    self.__hm_diff.Fill(sub1.chs[iiphi][iieta][id].means[0] 
                                        +sub1.chs[iiphi][iieta][id].means[1]
                                        +sub1.chs[iiphi][iieta][id].means[2]
                                        +sub1.chs[iiphi][iieta][id].means[3]
                                        -sub2.chs[iiphi][iieta][id].means[0]
                                        -sub2.chs[iiphi][iieta][id].means[1]
                                        -sub2.chs[iiphi][iieta][id].means[2]
                                        -sub2.chs[iiphi][iieta][id].means[3])
                    self.__hm0_diff.Fill(sub1.chs[iiphi][iieta][id].means[0]-sub2.chs[iiphi][iieta][id].means[0])
                    self.__hm1_diff.Fill(sub1.chs[iiphi][iieta][id].means[1]-sub2.chs[iiphi][iieta][id].means[1])
                    self.__hm2_diff.Fill(sub1.chs[iiphi][iieta][id].means[2]-sub2.chs[iiphi][iieta][id].means[2])
                    self.__hm3_diff.Fill(sub1.chs[iiphi][iieta][id].means[3]-sub2.chs[iiphi][iieta][id].means[3])
                    self.__hr0_diff.Fill(sub1.chs[iiphi][iieta][id].rmss[0]-sub2.chs[iiphi][iieta][id].rmss[0])
                    self.__hr1_diff.Fill(sub1.chs[iiphi][iieta][id].rmss[1]-sub2.chs[iiphi][iieta][id].rmss[1])
                    self.__hr2_diff.Fill(sub1.chs[iiphi][iieta][id].rmss[2]-sub2.chs[iiphi][iieta][id].rmss[2])
                    self.__hr3_diff.Fill(sub1.chs[iiphi][iieta][id].rmss[3]-sub2.chs[iiphi][iieta][id].rmss[3])
                    self.__hm0_sub1.Fill(sub1.chs[iiphi][iieta][id].means[0])
                    self.__hm1_sub1.Fill(sub1.chs[iiphi][iieta][id].means[1])
                    self.__hm2_sub1.Fill(sub1.chs[iiphi][iieta][id].means[2])
                    self.__hm3_sub1.Fill(sub1.chs[iiphi][iieta][id].means[3])
                    self.__hr0_sub1.Fill(sub1.chs[iiphi][iieta][id].rmss[0])
                    self.__hr1_sub1.Fill(sub1.chs[iiphi][iieta][id].rmss[1])
                    self.__hr2_sub1.Fill(sub1.chs[iiphi][iieta][id].rmss[2])
                    self.__hr3_sub1.Fill(sub1.chs[iiphi][iieta][id].rmss[3])
                    self.__hm0.Fill(rm0)
                    self.__hm1.Fill(rm1)
                    self.__hm2.Fill(rm2)
                    self.__hm3.Fill(rm3)
                    self.__hr0.Fill(rr0)
                    self.__hr1.Fill(rr1)
                    self.__hr2.Fill(rr2)
                    self.__hr3.Fill(rr3)
                    self.__hrawId.Fill(rraw)
                    
                    iphi = sub1.iphi(iiphi)
                    ieta = sub1.ieta(iieta)
                    self.__vh2rm0_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].means[0]-sub2.chs[iiphi][iieta][id].means[0])
                    self.__vh2rm1_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].means[1]-sub2.chs[iiphi][iieta][id].means[1])
                    self.__vh2rm2_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].means[2]-sub2.chs[iiphi][iieta][id].means[2])
                    self.__vh2rm3_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].means[3]-sub2.chs[iiphi][iieta][id].means[3])
                    self.__vh2rr0_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].rmss[0]-sub2.chs[iiphi][iieta][id].rmss[0])
                    self.__vh2rr1_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].rmss[1]-sub2.chs[iiphi][iieta][id].rmss[1])
                    self.__vh2rr2_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].rmss[2]-sub2.chs[iiphi][iieta][id].rmss[2])
                    self.__vh2rr3_diff[id].Fill(ieta, iphi, sub1.chs[iiphi][iieta][id].rmss[3]-sub2.chs[iiphi][iieta][id].rmss[3])
                    self.__vh2rm0[id].Fill(ieta, iphi, rm0)
                    self.__vh2rm1[id].Fill(ieta, iphi, rm1)
                    self.__vh2rm2[id].Fill(ieta, iphi, rm2)
                    self.__vh2rm3[id].Fill(ieta, iphi, rm3)
                    self.__vh2rr0[id].Fill(ieta, iphi, rr0)
                    self.__vh2rr1[id].Fill(ieta, iphi, rr1)
                    self.__vh2rr2[id].Fill(ieta, iphi, rr2)
                    self.__vh2rr3[id].Fill(ieta, iphi, rr3)

        #   Fill the Summary Histo
        mean_m0         = self.__hm0.GetMean()
        rms_m0          = self.__hm0.GetRMS()
        mean_m1         = self.__hm1.GetMean()
        rms_m1          = self.__hm1.GetRMS()
        mean_m2         = self.__hm2.GetMean()
        rms_m2          = self.__hm2.GetRMS()
        mean_m3         = self.__hm3.GetMean()
        rms_m3          = self.__hm3.GetRMS()
        
        mean_r0         = self.__hr0.GetMean()
        rms_r0          = self.__hr0.GetRMS()
        mean_r1         = self.__hr1.GetMean()
        rms_r1          = self.__hr1.GetRMS()
        mean_r2         = self.__hr2.GetMean()
        rms_r2          = self.__hr2.GetRMS()
        mean_r3         = self.__hr3.GetMean()
        rms_r3          = self.__hr3.GetRMS()
        
        mean_raw        = self.__hrawId.GetMean()
        rms_raw         = self.__hrawId.GetRMS()

        self.__hSummary.Fill(0, 0, mean_m0)
        self.__hSummary.Fill(0, 1, rms_m0)
        self.__hSummary.Fill(1, 0, mean_m1)
        self.__hSummary.Fill(1, 1, rms_m1)
        self.__hSummary.Fill(2, 0, mean_m2)
        self.__hSummary.Fill(2, 1, rms_m2)
        self.__hSummary.Fill(3, 0, mean_m3)
        self.__hSummary.Fill(3, 1, rms_m3)
        self.__hSummary.Fill(4, 0, mean_r0)
        self.__hSummary.Fill(4, 1, rms_r0)
        self.__hSummary.Fill(5, 0, mean_r1)
        self.__hSummary.Fill(5, 1, rms_r1)
        self.__hSummary.Fill(6, 0, mean_r2)
        self.__hSummary.Fill(6, 1, rms_r2)
        self.__hSummary.Fill(7, 0, mean_r3)
        self.__hSummary.Fill(7, 1, rms_r3)
        self.__hSummary.Fill(8, 0, mean_raw)
        self.__hSummary.Fill(8, 1, rms_raw)

#       self.__hSummary.GetZaxis().SetRangeUser(0.8, 1.2)

#------------------------------------------ 
#   Init Input
#------------------------------------------ 
fileName1 = sys.argv[1]
fileName2 = sys.argv[2]
rootFileName = sys.argv[3]
rootFile = rr.TFile(rootFileName, "recreate")
verbosity = 0

#------------------------------------------ 
#   Initialize the Subsystems and StatsMakers
#------------------------------------------ 
HF1 = HcalSubsystem("HF", 29, 41, 1, 71, 2, 1, 2, verbosity)
HF2 = HcalSubsystem("HF", 29, 41, 1, 71, 2, 1, 2, verbosity)
HB1 = HcalSubsystem("HB", 1, 16, 1, 72, 1, 1, 2, verbosity)
HB2 = HcalSubsystem("HB", 1, 16, 1, 72, 1, 1, 2, verbosity)
QIE111 = HcalSubsystem("QIE11", 16, 29, 1, 72, 1, 1, 7, verbosity)
QIE112 = HcalSubsystem("QIE11", 16, 29, 1, 72, 1, 1, 7, verbosity)
HE1 = HcalSubsystem("HE", 16, 29, 1, 72, 1, 1, 3, verbosity)
HE2 = HcalSubsystem("HE", 16, 29, 1, 72, 1, 1, 3, verbosity)
HO1 = HcalSubsystem("HO", 1, 15, 1, 72, 1, 4, 4, verbosity)
HO2 = HcalSubsystem("HO", 1, 15, 1, 72, 1, 4, 4, verbosity)

subs1 = {"HB" : HB1, "HE" : HE1, "HO" : HO1, "HF" : HF1, "QIE11" : QIE111}
subs2 = {"HB" : HB2, "HE" : HE2, "HO" : HO2, "HF" : HF2, "QIE11" : QIE112}

smHB = StatsMaker(rootFile, "HB")
smHE = StatsMaker(rootFile, "HE")
smHO = StatsMaker(rootFile, "HO")
smHF = StatsMaker(rootFile, "HF")
smQIE11 = StatsMaker(rootFile, "QIE11")

#------------------------------------------ 
#   Parse both files and produce the comparison
#------------------------------------------ 
parse(fileName1, subs1)
parse(fileName2, subs2)
smHB.makeStats(subs1["HB"], subs2["HB"])
smHE.makeStats(subs1["HE"], subs2["HE"])
smHO.makeStats(subs1["HO"], subs2["HO"])
smHF.makeStats(subs1["HF"], subs2["HF"])
smQIE11.makeStats(subs1["QIE11"], subs2["QIE11"])


#------------------------------------------ 
#   Finalize
#------------------------------------------ 
rootFile.Write()
rootFile.Close()
print '### DONE!'
















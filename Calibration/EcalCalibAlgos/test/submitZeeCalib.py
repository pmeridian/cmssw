#! /usr/bin/env python

import os
import sys
import optparse

parser = optparse.OptionParser("submitZeeCalib.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       help='run in batch in queue specified as option (default -q 1nd)', default='1nd')
(opt, args) = parser.parse_args()

os.system("mkdir -p batchCalib")
pwd = os.environ['PWD']

command = "cmsRun /afs/cern.ch/work/c/crovelli/calibZee/CMSSW_7_2_1/src/Calibration/EcalCalibAlgos/test/zeeCalibration_cfg_testDY.py"
print "submit "+command
logfile = "zeeCalib.log"
print logfile
outputname = "batchCalib/zeeCalib.src"
outputfile = open(outputname,'w')
outputfile.write('#!/bin/bash\n')
outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
outputfile.write('cd /afs/cern.ch/work/c/crovelli/calibZee/CMSSW_7_2_1/src/\n')
outputfile.write('eval `scramv1 runtime -sh`\n')
outputfile.write('cd -\n')
outputfile.write(command+"\n")
outputfile.write('cp zeeRescale*.root my*root *xml Zee*txt '+pwd+'/\n') 
print outputname 
os.system("bsub -q "+opt.queue+" -o batchCalib/"+logfile+" source "+pwd+"/"+outputname)

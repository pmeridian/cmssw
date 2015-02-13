#! /usr/bin/env python

import os
import sys
import optparse
import time

parser = optparse.OptionParser("submitZeeCalib.py")
parser.add_option('-q', '--queue',   action='store',  dest='queue',  help='run in batch in queue specified as option (default -q 1nd)', default='1nd')
parser.add_option('-c', '--cfg',     action='store',  dest='cfg',    help='cfg file', default='zeeCalibration_cfg.py')
parser.add_option('-t', '--task',    action='store',  dest='task',   help='taskName', default='zCalib')
parser.add_option('-r', '--run',     action='store',  dest='run',    help='do u want to run?', default=0)

(opt, args) = parser.parse_args()
os.system("mkdir -p "+opt.task)
pwd = os.environ['PWD']
user = os.environ['USER']

command = "cmsRun "+pwd+"/"+opt.cfg
print "submit "+command
logfile = pwd+"/"+opt.task+"/zeeCalib.log"
print logfile

outputname = pwd+"/"+opt.task+"/zeeCalib.sh"

outputfile = open(outputname,'w')

outputfile.write('#!/bin/bash\n')
outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
outputfile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
outputfile.write('cd '+pwd+'\n')
outputfile.write('eval `scramv1 runtime -sh`\n')
outputfile.write('if [ "${WORKDIR}" == "" ]; then\n')
outputfile.write('WORKDIR=/tmp/'+user+'/'+opt.task+'\n')
outputfile.write('mkdir -p ${WORKDIR}\n')
outputfile.write('fi\n')
outputfile.write('cd ${WORKDIR}\n')
outputfile.write(command+"\n")
outputfile.write('cp zeeRescale*.root my*root *xml Zee*txt '+pwd+'/'+opt.task+'/\n') 

print outputname+" created. Launch job. Run "+opt.run 

#time.sleep(1)

if ( (int(opt.run) == 1) and (opt.queue == "local")):
    commandLocal="chmod +x "+outputname+"; "+outputname
    print "Launching local "+commandLocal
    os.system(commandLocal)
elif ( (int(opt.run) == 1) ):
    os.system("bsub -q "+opt.queue+" -o "+logfile+" source "+outputname)



#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Launcher script on a cluster HPC. 
import os
#import glob
#import re
import signal
#import functools
#import ConfigParser
#import traceback
from collections import defaultdict
#from __future__ import print_function
import sys
from subprocess import check_call
try:
    from commands import getoutput
except ModuleNotFoundError:
    from subprocess import getoutput


class script(object):

    def __init__(self, cmd, afterID=None):
        if isinstance(cmd, list):
            self.cmd = ' '.join(cmd)
        else:
            self.cmd = cmd
        
        def _flat(L):
            res = []
            for x in L:
                if isinstance(x, list):
                    res += _flat(x)
                elif x != None:
                    res.append(x)
            return res
        
        if afterID:
            if isinstance(afterID, list):
                afterID = ':'.join(_flat(afterID))
        self.depend = afterID

    def __str__(self):
        return self.cmd

    def localRun(self):
        os.system(self.cmd)

    def _add_log(self, logFile, logStep):
        if logStep is None:
            logStep = self.cmd.split()[0]
        self.cmd += """\nif [ $? -ne 0 ]; then\n\texport DATE=`date "+%Y-%m-%d %H:%M:%S"`; echo "{logStep}\tfailed($DATE)" >> {logFile}\n\texit 1\nfi\nexport DATE=`date "+%Y-%m-%d %H:%M:%S"`; echo - | awk -v S=$SECONDS '{{printf "{logStep}\\\\tend(%s)\\\\t%02d:%02d:%02d\\\\n",ENVIRON["DATE"],S/(60*60),S%(60*60)/60,S%60}}' >> {logFile}\n""".format(logStep=logStep, logFile=logFile)

    def qsub(self, to_file, queue, ppn=2, logFile=None, logStep=None):
        if logFile is None:
            logFile = to_file + '.log'
        self._add_log(logFile, logStep)
        with open(to_file, 'w') as out:
            out.write(self.cmd)
        if self.depend:
            return getoutput('qsub -q {queue} -l nodes=1:ppn={ppn} -W depend=afterok:{afterID} -e {sh}.err -o {sh}.out {sh}'.format(queue=queue, ppn=ppn, afterID=self.depend, sh=to_file))
        else:
            return getoutput('qsub -q {queue} -l nodes=1:ppn={ppn} -e {sh}.err -o {sh}.out {sh}'.format(queue=queue, ppn=ppn, sh=to_file))

def getClusterValue(clusterOptions):  # -file=/path/to -queue=new_cu -ppn=2 -depend=ID
    h = defaultdict(str)
    for opt in clusterOptions.split():
        for k,v in opt.split('='):
            h[k] = v
    return h['-file'], h['-queue'], h['-ppn'], h['-depend']

def runCommand(cmd, clusterOptions, dryrun):
    if dryrun:
        print("\nDry run:\n")
        # Display environment variables for dry run
        if len(os.environ) != 0:
            print("    Env:\n")
            print( '\n'.join(['        %s = %s' % (v, os.environ[v]) for v in os.environ]) )
        if clusterOptions:
            print("\n    clusterOptions:\n")
            print("        %s\n" % clusterOptions)
        print("    Cmd:\n")
        print(("        " + " ".join(cmd)+"\n"))
    else:
        sys.stderr.write( "\nRunning:\n")
        sys.stderr.write("    " + " ".join(cmd)+"\n")
        if clusterOptions is not None:
            to_file, queue, ppn, depend = getClusterValue(clusterOptions)
            cmd = script(cmd, afterID=depend)
            cmd.qsub(to_file, queue, ppn)
        else:
            check_call(cmd)

class easyRunLaunchException(Exception):
    pass

def getValueForArgument(args, argument):
    if argument in args:
        i = args.index(argument)
        if len(args) <= i+1:
            raise easyRunLaunchException("Argument: " + argument + " requires a parameter")
        return args[i+1]
    return None

def signal_handler(signal, frame):
    sys.exit(1)

def main(args):
    #suppress stack trace when killed by keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    try:
        if len(args) is 0 or (len(args) is 1 and (args[0] == "--help" or args[0] == "-h")):
            print("")
            print(' Usage example: easyRun gatk --cluster-options "-file=/path/to -queue=new_cu -ppn=2 -depend=ID" --java-options "-Xmx20g" MarkDuplicates ...')
            print("")
            print(" Getting help")
            print("    easyRun --list       Print the list of available tools" )
            print("")
            print("    easyRun Tool --help  Print help on a particular tool" )
            print("")
            print("   --dry-run             may be specified to output the generated command line without running it")
            print("   --cluster-options     'OPTION1[ OPTION2=Y ... ]'   optional - pass the given string of options to the ")
            print("                         cluster at runtime.  ")
            print("")
            sys.exit(0)

        if len(args) is 1 and args[0] == "--list":
            args[0] = "--help"  # if we're invoked with --list, invoke the GATK with --help

        dryRun = "--dry-run" in args
        if dryRun:
            dryRun = True
            args.remove("--dry-run")
        
        clusterOptions = getValueForArgument(args, "--cluster-options")
        if clusterOptions is not None:
            i = args.index("--cluster-options")
            del args[i] #remove javaOptions
            del args[i] #and its parameter

        runCommand(args, clusterOptions, dryRun)

    except easyRunLaunchException as e:
        sys.stderr.write(str(e)+"\n")
        sys.exit(3)


if __name__ == "__main__":
    main(sys.argv[1:])

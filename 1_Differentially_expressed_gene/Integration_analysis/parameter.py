## Don't edit this script.
## This code is imported into other scripts and used.

import os, sys, commands
from subprocess import Popen, PIPE, STDOUT
## make_directory
def make_dir(name):
    if not os.path.exists(name): os.mkdir(name)
## execute_in_Shell
def inShell(command):
    print command
    p = Popen(command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()


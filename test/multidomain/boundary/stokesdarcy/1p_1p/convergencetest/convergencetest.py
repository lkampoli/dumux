#!/usr/bin/env python

from math import *
import subprocess
import sys

if len(sys.argv) < 2:
    sys.stderr.write('Please provide a single argument <testname> to the script\n')
    sys.exit(1)

testname = str(sys.argv[1])
testargs = [str(i) for i in sys.argv][2:]

# remove the old log files
subprocess.call(['rm', testname + '_stokes.log'])
print("Removed old log file ({})!".format(testname + '_stokes.log'))
subprocess.call(['rm', testname + '_darcy.log'])
print("Removed old log file ({})!".format(testname + '_darcy.log'))

# do the runs with different refinement
for i in [0, 1, 2]:
    subprocess.call(['./' + testname] + testargs + ['-Grid.Refinement', str(i),
                                      '-Vtk.OutputName', testname])

def checkRatesStokes():
    # check the rates and append them to the log file
    logfile = open(testname + '_stokes.log', "r+")

    errorP = []
    errorVx = []
    errorVy = []
    for line in logfile:
        line = line.strip("\n")
        line = line.strip("\[ConvergenceTest\]")
        line = line.split()
        errorP.append(float(line[2]))
        errorVx.append(float(line[5]))
        errorVy.append(float(line[8]))

    resultsP = []
    resultsVx = []
    resultsVy = []
    logfile.truncate(0)
    logfile.write("n\terrorP\t\trateP\t\terrorVx\t\trateVx\t\terrorVy\t\trateVy\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(errorP)-1):
        if isnan(errorP[i]) or isinf(errorP[i]):
            continue
        if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12) and (errorVx[i] < 1e-12 or errorVx[i+1] < 1e-12) and (errorVy[i] < 1e-12 or errorVy[i+1] < 1e-12)):
            rateP = (log(errorP[i])-log(errorP[i+1]))/log(2)
            rateVx = (log(errorVx[i])-log(errorVx[i+1]))/log(2)
            rateVy = (log(errorVy[i])-log(errorVy[i+1]))/log(2)
            message = "{}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP,  errorVx[i], rateVx, errorVy[i], rateVy)
            logfile.write(message)
            resultsP.append(rateP)
            resultsVx.append(rateVx)
            resultsVy.append(rateVy)
        else:
            logfile.write("error: exact solution!?")
    i = len(errorP)-1
    message = "{}\t{:0.4e}\t\t{}\t{:0.4e}\t\t{}\t{:0.4e}\t\t{}\n".format(i, errorP[i], "",  errorVx[i], "", errorVy[i], "")
    logfile.write(message)

    logfile.close()
    print("\nComputed the following convergence rates for {}:\n".format(testname))

    subprocess.call(['cat', testname + '_stokes.log'])

    return {"p" : resultsP, "v_x" : resultsVx, "v_y" : resultsVy}

def checkRatesDarcy():
    # check the rates and append them to the log file
    logfile = open(testname + '_darcy.log', "r+")

    errorP = []
    for line in logfile:
        line = line.strip("\n")
        line = line.strip("\[ConvergenceTest\]")
        line = line.split()
        errorP.append(float(line[2]))

    resultsP = []
    logfile.truncate(0)
    logfile.write("n\terrorP\t\trateP\n")
    logfile.write("-"*50 + "\n")
    for i in range(len(errorP)-1):
        if isnan(errorP[i]) or isinf(errorP[i]):
            continue
        if not ((errorP[i] < 1e-12 or errorP[i+1] < 1e-12)):
            rateP = (log(errorP[i])-log(errorP[i+1]))/log(2)
            message = "{}\t{:0.4e}\t{:0.4e}\n".format(i, errorP[i], rateP)
            logfile.write(message)
            resultsP.append(rateP)
        else:
            logfile.write("error: exact solution!?")
    i = len(errorP)-1
    message = "{}\t{:0.4e}\n".format(i, errorP[i], "")
    logfile.write(message)

    logfile.close()
    print("\nComputed the following convergence rates for {}:\n".format(testname))

    subprocess.call(['cat', testname + '_darcy.log'])

    return {"p" : resultsP}

def checkRatesStokesAndDarcy():
    resultsStokes = checkRatesStokes()
    resultsDarcy = checkRatesDarcy()

    def mean(numbers):
        return float(sum(numbers)) / len(numbers)

    # check the rates, we expect rates around 2
    if mean(resultsStokes["p"]) < 2.05 and mean(resultsStokes["p"]) < 1.84:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

    if mean(resultsStokes["v_x"]) < 2.05 and mean(resultsStokes["v_x"]) < 1.95:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for x-velocity were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

    if mean(resultsStokes["v_y"]) < 2.05 and mean(resultsStokes["v_y"]) < 1.95:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for y-velocity were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)

    if mean(resultsDarcy["p"]) < 2.05 and mean(resultsDarcy["p"]) < 1.95:
        sys.stderr.write("*"*70 + "\n" + "The convergence rates for pressure were not close enough to 2! Test failed.\n" + "*"*70 + "\n")
        sys.exit(1)


checkRatesStokesAndDarcy()

sys.exit(0)

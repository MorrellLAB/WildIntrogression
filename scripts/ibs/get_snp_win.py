#!/usr/bin/env python3

#   Chaochih Liu

import sys
import os


def readSnps(snp_file):
    """Function that reads in a single column of SNPs"""
    #   Use list comprehension to read in file line by line
    #   Removing trailing new line
    dat = [line.strip('\n') for line in open(snp_file, 'r').readlines()]
    return dat


def getStepList(d, win, step):
    step_list = []
    for s in range(0, len(d), step):
        if s + win <= len(d):
            step_list.append(s)
    return(step_list)


def getSnpWin(d, step_list, win):
    new_list = []
    for s in step_list:
        tmp = []
        for i in range(s, s + win):
            tmp.append(d[i])
        new_list.append(tmp)
    return(new_list)


def getLastStep(d, step_list, win, step):
    last_win = []
    l = len(step_list) - 1
    if step_list[l] + step + win > len(d):
        next_step = step_list[l] + step
        last_win.append(next_step)
        last_win.append(len(d))
    return(last_win)


def getLastWin(d, last_interval):
    new_list = []
    tmp = []
    for i in range(last_interval[0], last_interval[1]):
        tmp.append(d[i])
    new_list.append(tmp)
    return(new_list)


def outFile(interval, output_file):
    with open(output_file, 'a', newline='') as out:
        out.write(interval + '\n')


def main(filepath, win_size, step_size, out_dir):
    #   User provided input arguments
    fp = filepath
    winSize = int(win_size)
    stepSize = int(step_size)
    outDir = out_dir

    #   These will be generated from arguments above
    prefix = os.path.basename(fp).split('.')[0]
    suffix = ".txt"

    #   Call on functions to do work
    data = readSnps(snp_file=fp)
    step_l = getStepList(d = data, win = winSize, step = stepSize)
    step_w = getSnpWin(d = data, step_list = step_l, win = winSize)
    last_step = getLastStep(d = data, step_list = step_l, win = winSize, step = stepSize)
    last_win = getLastWin(d = data, last_interval = last_step)

    for i in range(0, len(step_l)):
        out_name = outDir + '/' + prefix + '_' + str(step_l[i]) + '-' + str(step_l[i] + winSize) + suffix
        with open(out_name, 'a') as out:
            for s in range(0, len(step_w[i])):
                l = step_w[i][s]
                out.write(l + '\n')


    last_step_out_name = outDir + '/' + prefix + '_' + str(last_step[0]) + '-' + str(last_step[1]) + suffix
    with open(last_step_out_name, 'a') as out:
        for s in range(0, len(last_win[0])):
            l = last_win[0][s]
            out.write(l + '\n')


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) # Run the program

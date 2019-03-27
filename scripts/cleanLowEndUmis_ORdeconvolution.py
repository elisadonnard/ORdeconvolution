import sys
import getopt
import os.path
import glob
import argparse

# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<input UMI distributions>',
                        help='input UMI distributions file')
    parser.add_argument('-o', type=str, required=True, metavar='<output UMI counts file>',
                        help='output UMI counts file')
    parser.add_argument('-r', type=str, required=False, default=None, metavar='<output read counts file>',
                        help='output read counts file')
    parser.add_argument('-n', type=int, required=False, default=2, metavar='<minimum UMI count for merging>',
                        help='minimum required UMI count for merging singletons')
    parser.add_argument('-u', type=int, required=False, default=1, metavar='<minimum UMIs per cell>',
                        help='minimum total UMIs required to keep a cell')

    args = parser.parse_args()
    return args

def hammingDist(x, y):
    # find the Hamming distance between two input strings:
    if len(x)!=len(y):
        hd = len(x)
    else:
        hd = 0
        for i in range(len(x)):
            if x[i]!=y[i]:
                hd+=1    # count mismatches
    return hd

def testUmiDist(singList, multList):
    ## test each element of singList against all elements of multList to see if any
    ## singletons are one base off from any of the multis.

    badUmis = []   # output list of UMIs to be removed
    repUmis = {}   # UMIs to be incremented, with the increment count

    for bc0 in singList:
        bc1off = []   # barcodes with a Hamming distance of 1 from the singelton
        for bc1 in multList:
            if hammingDist(bc0, bc1)==1:
                bc1off.append(bc1)  # add to the list of one-offs
        # if a barcode is one off from any multi, remove it:
        if len(bc1off)>0:
            badUmis.append(bc0)
        # if bc0 is one off from only ONE multi, increment the count for that multi:
        if len(bc1off)==1:
            repUmis.setdefault(bc1off[0],0)
            repUmis[bc1off[0]]+=1

    return [badUmis, repUmis]
    
def sumValues(d):
    ## sums the values of the values in a key:value list:
    x=0
    for key,value in d:
        x+=value
    return x

def run(args):
    inFile = args.i      # input file name
    outFile = args.o     # output file name
    readsFile = args.r   # optional output reads file
    nMin = args.n        # minimum number of UMI counts for merging with singletons
    uMin = args.u        # minimum number of UMIs per cell

    print('inFile: %s\n' % inFile)
    print('outFile: %s\n' % outFile)
    print('readsFile: %s\n' % readsFile)

    fIn = open(inFile, 'r')
    umiDict = {}
    umiHist = {}
    umiNew = {}

    gDict = {}   # dictionary whose keys are all observed genes
    nLines = 0

    while 1:
        line = fIn.readline()
        if not line:
            break

        if line.endswith('\n'):
            line=line[:-1]

        nLines+=1
        # split into barcode and count:
        fields = line.split('\t')

        # signal bad input line:
        if len(fields)<4:
            print('Too few items in line: %s' % line)

        # parse the line:
        bc = fields[0]
        g = fields[1]
        umi = fields[2]
        n = int(fields[3])

        # update the gene dictionary:
        gDict.setdefault(g,0)

        # UMI counts histogram:
        umiHist.setdefault(bc,{})
        umiHist[bc].setdefault(g,{})
        umiHist[bc][g].setdefault(umi,0)
        umiHist[bc][g][umi]+=n

        # UMI counts inverted gene:
        umiNew.setdefault(bc,{})
        umiNew[bc].setdefault(umi,{})
        umiNew[bc][umi].setdefault(g,0)
        umiNew[bc][umi][g]+=n

        ### Keep track of total UMIs and total reads per gene for each cell
        umiHist[bc][g].setdefault("totalUMI",0)
        umiHist[bc][g].setdefault("totalReads",0)
        umiHist[bc][g]["totalUMI"]+=1
        umiHist[bc][g]["totalReads"]+=n
        ### END 
            
        # keep stats on all UMIs:
        umiDict.setdefault(umi,0)
        umiDict[umi]+=n

    fIn.close()

    ## get lists of UMIs with one count and UMIs with >nMin counts:
    for bc in umiHist.keys():
        for g in umiHist[bc].keys():

            # for average UMI count calculations:
            nzUmis = 0   # number of UMIs with non-zero counts

            # make lists of singlet and multis:
            singlets = []
            multis = []
            for umi in umiHist[bc][g].keys():
                if umiHist[bc][g][umi]==1:
                    singlets.append(umi)
                elif umiHist[bc][g][umi]>=nMin:
                    multis.append(umi)
                    nzUmis+=1
                else:
                    nzUmis+=1

            # separate true singlets from ones that have a Hamming distance of 1 from one or more of the multis:
            [rmUmis, incUmis] = testUmiDist(singlets, multis)

            # remove bad UMIs from the dictionary
            for i in range(len(rmUmis)):
                del umiHist[bc][g][rmUmis[i]]     # delete error UMI
                if len(umiNew[bc][rmUmis[i]]) == 1:
                    del umiNew[bc][rmUmis[i]]
                else:
                    del umiNew[bc][rmUmis[i]][g]


            # update the counts for UMIs that were uniquely one-off from one of the singletons:
            for u in incUmis.keys():
                umiHist[bc][g][u]+=incUmis[u]
                umiNew[bc][u][g] = umiHist[bc][g][u]
        ###
        identical = []
        hdist1 = []
        # For cells with > 1 OR detected
        if len(umiHist[bc].keys()) > 1:
            # correct identical umis in other OR
            for umi in umiNew[bc].keys():
                if len(umiNew[bc][umi].keys()) > 1:
                    identical.append(umi)
            # merge to top OR
            for umi in identical:
                genes = list(umiNew[bc][umi].keys())
                total = []
                for k in genes:
                    total.append(umiHist[bc][k]["totalReads"])
                m = total.index(max(total))
                topOR = genes[m]
                del genes[m]
                for rg in genes:
                    ### add the count and remove from umiHist 
                    umiHist[bc][topOR][umi]+=umiHist[bc][rg][umi]
                    del umiHist[bc][rg][umi]
                    umiNew[bc][umi][topOR] = umiHist[bc][topOR][umi]
                    del umiNew[bc][umi][rg]
                    
            # compare all umis and find 1 hamming dist apart
            tocompare = list(umiNew[bc].keys())
            for idx, u in enumerate(tocompare):
                if idx + 1 < len(tocompare):
                    updated = tocompare[idx+1:]
                    for e in updated:
                        if hammingDist(u, e) == 1:                    
                            hdist1.append([u, e])
            # merge to top OR
            for umi in hdist1:
                gene1 = umiNew[bc][umi[0]].keys()[0]
                gene2 = umiNew[bc][umi[1]].keys()[0]
                if umiHist[bc][gene1]["totalReads"] > umiHist[bc][gene2]["totalReads"]:
                    umiHist[bc][gene1][umi[0]]+=umiHist[bc][gene2][umi[1]]
                    del umiHist[bc][gene2][umi[1]]
                elif umiHist[bc][gene1]["totalReads"] < umiHist[bc][gene2]["totalReads"]:
                    umiHist[bc][gene2][umi[1]]+=umiHist[bc][gene1][umi[0]]
                    del umiHist[bc][gene1][umi[0]]

    ###
    with open(outFile, 'w') as f:
        for c in umiNew.keys():
            for u in umiNew[c].keys():
                for g, v in umiNew[c][u].items():
                    f.write('%s\t%s\t%s\t%s\n' % (c, g, u, v))
                    
################################################################################
#    ## Alan's output matrix
#    fOut = open(outFile, 'w')
#
#    ## open the optional read counts file:
#    if readsFile==None:
#        writeReads=False
#    else:
#        writeReads=True
#        fReads = open(readsFile, 'w')
#
#    ## remove any barcodes/cells with fewer than uMin total UMIs (if uMin>0):
#    if uMin>0:
#        bcRm = []
#        for bc in umiHist.keys():   # loop over barcodes
#            bcSum = 0
#            for g in umiHist[bc].keys():   # loop over genes
#                bcSum += umiHist[bc][g]["totalUMI"]  # sum total UMIs for this cell
#            ## remove this barcode if the total count is less than uMin:
#            if bcSum<uMin:
#                bcRm.append(bc)
#    ## remove low-count cells:
#    for bc in bcRm:
#        del umiHist[bc]
#            
#    ## file header:
#    hStr = 'gene'      # header
#    bcList = umiHist.keys()     # just to make absolutely sure they stay in the correct order
#    for bc in bcList:
#        hStr = '%s\t%s' % (hStr,bc)
#    fOut.write('%s\n' % hStr)
#
#    ## option reads count file header
#    if writeReads:
#        fReads.write('%s\n' % hStr)
#
#    # per-gene count:
#    for g in gDict.keys():  # loop over all observed genes
#        oStr = '%s' % g
#        rStr = '%s' % g
#        for bc in bcList:   # UMI counts for each cell
#            if umiHist[bc].has_key(g):
#                oStr = '%s\t%d' % (oStr,umiHist[bc][g]["totalUMI"])
#                rStr = '%s\t%d' % (rStr,umiHist[bc][g]["totalReads"])
#            else:            
#                oStr = '%s\t0' % oStr
#                rStr = '%s\t0' % rStr
#
#        fOut.write('%s\n' % oStr)        
#        if writeReads:
#            fReads.write('%s\n' % rStr)        
#
#    fOut.close()
#    if writeReads:
#        fReads.close()
################################################################################

    return                      

if __name__ == "__main__":
    args = get_args()
    run(args)



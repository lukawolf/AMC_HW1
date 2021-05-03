import argparse
import sys
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

# Takes the depthFile path (string), chromosome (string) and posRange (int, int) as arguments and returns a numpy array containing the read depth
def parseDepthFile(depthFile, chromosome, posRange):
    toReturn = np.zeros(posRange[1] - posRange[0], dtype = int)
    with open(depthFile) as depthRows:
        for depthRow in depthRows:
            readChromosome, position, depth = depthRow.split()
            if readChromosome != chromosome:
                raise Exception('Unexpected chromosome identifier in the depth file.')
            intPosition = int(position)
            if intPosition < posRange[0] or intPosition > posRange[1]:
                continue
            # Positions in depth file are 1 based
            toReturn[intPosition - 1 - posRange[0]] = int(depth)
    return toReturn

# Takes the depthFiles paths ([string]), chromosome (string), posRange (int, int) and window (int) as arguments and creates depth plot files
def plotDepth(depthFiles, chromosome, posRange, window):
    data = []
    for depthFile in depthFiles:
        data.append(parseDepthFile(depthFile, chromosome, posRange))

    title = "Read depth of (" + chromosome + ":" + str(posRange[0]) + ":" + str(posRange[1]) + ")"

    if window > 1:
        splitPositions = list(range(window, len(data[1]), window))
        for i in range(0, len(data)):
            splitData = np.split(data[i], splitPositions)
            averaged = np.zeros(len(splitData))
            for y in range(0, len(splitData)):
                averaged[y] = np.average(splitData[y])
            data[i] = averaged
        title += " averaged each " + str(window) + " bp"

    sb.set(color_codes=True)
    for i in range(0, len(depthFiles)):
        plt.title(title + " - " + depthFiles[i].split(".")[0])
        ax = plt.subplot(1, 1, 1)
        sbPlot = sb.scatterplot(x=range(posRange[0], posRange[1], window), y=data[i], s=5)
        sbPlot.set(xlabel='Position (bp)', ylabel="Depth")
        plt.savefig(depthFiles[i] + ".png", bbox_inches='tight', dpi=400)
        plt.close()

    plt.title(title + " - log2 ratio")
    ax = plt.subplot(1, 1, 1)
    # We need to avoid division by 0
    sbPlot = sb.scatterplot(x=range(posRange[0], posRange[1], window), y=np.log2(np.divide(data[0], data[1] + np.finfo(float).eps)), s=5)
    sbPlot.set(xlabel='Position (bp)', ylabel="log2 ratio")
    plt.savefig("ratio.png", bbox_inches='tight', dpi=400)
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--chromosome", type=str, required=True, help="The chromosome identifier")
parser.add_argument("--fromPos", type=int, required=True, help="The position from which we want to plot the depth")
parser.add_argument("--toPos", type=int, required=True, help="The position to which we want to plot the depth")
parser.add_argument("--window", type=int, default=1, help="The window size over which read depth is averaged. Note that this is not a sliding window, but a chunk")
parser.add_argument("files", metavar='path', type=str, nargs=2, help='The depth file paths to plot. The first file should be disease, the second one being wt')
args = parser.parse_args()

plotDepth(args.files, args.chromosome, (args.fromPos, args.toPos), args.window)

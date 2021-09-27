import matplotlib.pyplot as plt
import pandas as pd

def plotToFile(infile, outfile):
    df = pd.read_csv(infile)
    plt.plot(df["x"], df["v1"], label="v1")
    plt.plot(df["x"], df["v2"], label="v2")
    plt.plot(df["x"], df["v3"], label="v3")
    plt.xlabel(r'$\hat{x}$')
    plt.ylabel("v")
    plt.legend()
    plt.savefig(outfile)
    plt.show()
    plt.close()

if __name__ == "__main__":
    plotToFile("output/n10.csv", "output/n10.pdf")
    plotToFile("output/n100.csv", "output/n100.pdf")


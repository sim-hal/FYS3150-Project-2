import matplotlib.pyplot as plt
import pandas as pd

def plotToFile(infile, outfile):
    df = pd.read_csv(infile)
    plt.plot(df["x"], df["jacobi1"], label="v1 (Jacobi)")
    plt.plot(df["x"], df["jacobi2"], label="v2 (Jacobi)")
    plt.plot(df["x"], df["jacobi3"], label="v3 (Jacobi)")
    plt.plot(df["x"], df["analytic1"], label="v1 (Analytic)")
    plt.plot(df["x"], df["analytic2"], label="v2 (Analytic)")
    plt.plot(df["x"], df["analytic3"], label="v3 (Analytic)")
    plt.xlabel(r'$\hat{x}$')
    plt.ylabel("v")
    plt.legend()
    plt.savefig(outfile)
    plt.show()
    plt.close()

if __name__ == "__main__":
    plotToFile("output/n10.csv", "output/n10.pdf")
    plotToFile("output/n100.csv", "output/n100.pdf")


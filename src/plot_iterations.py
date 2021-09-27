import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    df = pd.read_csv("output/iters.csv")
    plt.plot(df["N"], df["Iterations tridiagonal"], label="Tridiagonal matrix")
    plt.plot(df["N"], df["Iterations dense"], label="Dense matrix")
    plt.xlabel("N")
    plt.ylabel("#iterations")
    plt.legend()
    plt.savefig("output/iteration_plot.pdf")
    plt.show()

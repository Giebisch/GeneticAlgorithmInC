import matplotlib.pyplot as plt

if __name__ == "__main__":
    plt.subplots(1, 1)

    average = []
    best = []

    # get all values from .txt file
    with open("fitness_values.txt") as fopen:
        next(fopen)
        for line in fopen:
            avg, b = line.split("\t")
            average.append(float(avg.strip()))
            best.append(float(b.strip()))

    plt.plot(range(0, len(average)), average, label="Average Fitness")
    plt.plot(range(0, len(average)), best, label="Best Fitness")

    plt.title("Optimization of the Rosenbrock function")
    plt.legend()
    plt.xlabel("Generations")
    plt.ylabel("Fitness Value")
    plt.show()

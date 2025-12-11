import numpy as np
import matplotlib.pyplot as plt


def monte_carlo_area(n_points):
    # Generate points in [-1, 1] square
    points = np.random.uniform(-1, 1, (n_points, 2))
    # Calculate points inside unit circle
    inside = np.sum(points[:, 0]**2 + points[:, 1]**2 <= 1)
    # Area estimate and error
    area_estimate = 4 * inside / n_points
    error = abs(area_estimate - np.pi)
    return area_estimate, error

# Different sample sizes
sample_sizes = [10, 100, 1000, 10**4, 10**5, 10**6]
results = {n: monte_carlo_area(n) for n in sample_sizes}

# Results
print("\nResults:")
print("Samples\t\tEstimate\tError")
for n in sample_sizes:
    area, err = results[n]
    print(f"{n}\t\t{area:.5f}    \t{err:.5f}")

# Plot error progression
plt.figure(figsize=(8, 5))
plt.plot(results.keys(), [err for _, err in results.values()], 'o-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Samples')
plt.ylabel('Absolute Error')
plt.title('Monte Carlo Error Convergence')
plt.grid(True)
plt.show()

# Results plot
plt.figure(figsize=(8, 5))
plt.plot(sample_sizes, [area for area, _ in results.values()], 'o-', color='blue', label='Monte Carlo Estimate')
plt.axhline(np.pi, color='red', linestyle='--', label='π (True Value)')
plt.xscale('log')
plt.xlabel('Number of Samples')
plt.ylabel('Estimate of π')
plt.title('Monte Carlo π Estimate by Sample Size')
plt.grid(True, which="both", linestyle="--", alpha=0.7)
plt.legend()
plt.show()
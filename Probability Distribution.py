import numpy as np
import matplotlib.pyplot as plt

# Random numbers between 0 and 1
data = np.random.rand(1000)

# Create bins and calculate probabilities
counts, bin_edges = np.histogram(data, bins=10, range=(0, 1))
probabilities = counts / 1000

# Plot the probabilities
plt.bar(bin_edges[:-1], probabilities, width=0.1, align='edge', edgecolor='black')
plt.xlabel('Bins')
plt.ylabel('Probabilities')
plt.title('Probability Distribution of Bins')
plt.xticks(bin_edges[:-1] + 0.05, [f'{i:.1f}-{i+0.1:.1f}' for i in bin_edges[:-1]])
plt.show()

# Calculate variance 
avg_squared = np.mean(probabilities ** 2)
avg_prob = np.mean(probabilities)
variance = np.sqrt(avg_squared - avg_prob ** 2)

# Comparison value
comparison = 1 / np.sqrt(1000)

# Results
print("Probabilities for each bin:")
for i, prob in enumerate(probabilities):
    print(f"Bin {i+1} ({bin_edges[i]:.1f}-{bin_edges[i+1]:.1f}): {prob:.4f}")

print(f"\nVariance: {variance:.4f}")
print(f"1/sqrt(1000): {comparison:.4f}")
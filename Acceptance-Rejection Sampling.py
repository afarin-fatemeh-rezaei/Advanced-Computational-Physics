import numpy as np
import matplotlib.pyplot as plt

# p(y)
def p(y):
    return np.exp(-(y - 1)**2)

y_val = np.linspace(-2,4,1000)
area = np.trapezoid(p(y_val), y_val)
py_norm = lambda y: p(y)/area

# Acceptance-Rejection 
def acceptance_rejection(samples):
    accept = []
    while len(accept) < samples:
        y_test = np.random.uniform(-2, 4, size=samples)
        py_test = np.random.uniform(0, 1, size=samples)  
        accept.extend(y_test[py_test < p(y_test)])
    return np.array(accept[:samples])

# Sample sizes
samples= [100, 10000, 1000000]

# Acceptance-Rejection for each sample size
samples_100 = acceptance_rejection(samples[0])
samples_10000 = acceptance_rejection(samples[1])
samples_1000000 = acceptance_rejection(samples[2])

# Bins for the histogram
bins = 50
hist_range = (-2, 4)

# Compute histogram counts normalized to density for each sample set
counts_100, bin_edge = np.histogram(samples_100, bins=bins, range=hist_range, density=True)
counts_10000, bin_edge = np.histogram(samples_10000, bins=bins, range=hist_range, density=True)
counts_1000000, bin_edge = np.histogram(samples_1000000, bins=bins, range=hist_range, density=True)

# Bin centers
bin_center = (bin_edge[:-1] + bin_edge[1:]) / 2

# Calculate true PDF values at bin centers
true_pdf_values = py_norm(bin_center)

# Error analysis function
def calculate_errors(hist_counts, true_values):
    mae = np.mean(np.abs(hist_counts - true_values))
    rmse = np.sqrt(np.mean((hist_counts - true_values)**2))
    return mae, rmse

# Calculate errors for each sample size
mae100, rmse100 = calculate_errors(counts_100, true_pdf_values)
mae10k, rmse10k = calculate_errors(counts_10000, true_pdf_values)
mae1M, rmse1M = calculate_errors(counts_1000000, true_pdf_values)

# Plot p(x)
plt.figure(figsize=(10, 6))
plt.plot(y_val, py_norm(y_val), 'k-', linewidth=2, label=r'Normalized $p(x)=e^{-(x-1)^2}$')

# Plot histogram 
# 100 
plt.plot(bin_center, counts_100, marker='o', linestyle='None', markersize=8,
         markerfacecolor='none', markeredgecolor='green', label='Samples = 100')

# 10000 
plt.plot(bin_center, counts_10000, marker='o', linestyle='None', markersize=8,
         markerfacecolor='blue', markeredgecolor='blue', label='Samples = 10,000')

# 1000000
plt.plot(bin_center, counts_1000000, marker='*', linestyle='None', markersize=8,
         markerfacecolor='red', markeredgecolor='red', label='Samples = 1,000,000')


plt.xlabel('x')
plt.ylabel('Density')
plt.title('Acceptance-Rejection Sampling vs. Target Distribution')
plt.legend()
plt.show()

# Print error table
print("\nError Analysis:")
print(f"{'Sample Size':<12} | {'Mean Absolute Error':<8} | {'Root Mean Squared Error':<8}")
print("----------------------------------------")
print(f"{'100':<12} | {mae100:.4f}              | {rmse100:.4f}")
print(f"{'10,000':<12} | {mae10k:.4f}              | {rmse10k:.4f}")
print(f"{'1,000,000':<12} | {mae1M:.4f}              | {rmse1M:.4f}")
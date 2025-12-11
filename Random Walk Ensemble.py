import numpy as np
import matplotlib.pyplot as plt

# M for single-particle random walk
M_values = [100, 1000, 10000, 100000, 1000000, 10000000]
P = 10  

single_final_tavg = []
ensemble_final_avg = []
absolute_errors = []
rmse_errors = []

# Loop different M values
for M in M_values:
   
    N = M//P  

    # Part 1: Single Particle Random Walk 
    steps_single = np.random.choice([-1, 1], size=M)
    Avarage_Time = np.mean(steps_single)
    positions_single = np.cumsum(steps_single)
    single_final_tavg.append(Avarage_Time)

    # Part 2: N-Particle Random Walk (10 steps each)
    ensemble_steps = np.random.choice([-1, 1], size=(N, P))
    final_positions = np.sum(ensemble_steps, axis=1)
    ensemble_avg = np.sum(final_positions)/N
    ensemble_final_avg.append(ensemble_avg)

    # Calculate errors
    abs_error = np.abs(Avarage_Time - ensemble_avg)
    rmse = np.sqrt((Avarage_Time - ensemble_avg)**2)  # For single trial, RMSE = absolute error
    
    absolute_errors.append(abs_error)
    rmse_errors.append(rmse)

    #Results
    print(f"M = {M}\tP = {P}\tN = {N}")
    print("  Single Particle Time Avarage:", Avarage_Time)
    print("  Avarage of Final Positions for Ensemble (N = {}, 10 steps each): {}".format(N, ensemble_avg))
    print("-" * 50)

#Plot
plt.figure(figsize=(10, 6))
plt.plot(M_values, single_final_tavg, 'bo-', label='Single Particle Time Avarage')
plt.plot(M_values, ensemble_final_avg, 'rs-', label="Ensemble Final Position's Avarage")
plt.xscale('log')
plt.xlabel('Total Steps (M) [log scale]')
plt.ylabel('Time Avarage')
plt.title('Comparison:\nSingle Particle vs. Ensemble (M = P * N)')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()


# Error plots
plt.figure(figsize=(10, 4))

# Absolute Error Plot
plt.subplot(1, 2, 1)
plt.plot(M_values, absolute_errors, 'bo-', alpha=0.7)
plt.xscale('log')
plt.xlabel('Total Steps (M) [log scale]')
plt.ylabel('Error Magnitude')
plt.title('Absolute Error')
plt.grid(True, which="both", linestyle='--')

# RMSE Plot
plt.subplot(1, 2, 2)
plt.plot(M_values, rmse_errors, 'rs-', alpha=0.7)
plt.xscale('log')
plt.xlabel('Total Steps (M) [log scale]')
plt.ylabel('Error Magnitude')
plt.title('Root Mean Squared Error')
plt.grid(True, which="both", linestyle='--')

plt.tight_layout()
plt.show()

#Print final results
print("\n\nSummary Table:")
print("="*90)
print(f"{'M Value':<12} | {'Time Average':<15} | {'Ensemble Average':<18} | {'Absolute Error':<15} | {'RMSE':<15}")
print("="*90)
for M, t_avg, e_avg, abs_err, rmse in zip(M_values, single_final_tavg, 
                                         ensemble_final_avg, absolute_errors, rmse_errors):
    print(f"{M:<12} | {t_avg:<15.6f} | {e_avg:<18.6f} | {abs_err:<15.6f} | {rmse:<15.6f}")
print("="*90)
import numpy as np
import matplotlib.pyplot as plt




# Part 1: Single Particle Random Walk (M = 1000 steps)
M = 1000 

# Generate M steps: 
steps_single = np.random.choice([-1, 1], size=M)
Time_Avarage = np.mean(steps_single)
positions_single = np.cumsum(steps_single)

# Plot the single particle random walk:
plt.figure(figsize=(10, 4))
plt.plot(range(M), positions_single, color='blue', label='Particle Position')
plt.axhline(0, color='red', linestyle='--', lw=2, label='Starting Point (0)')
plt.axhline(Time_Avarage, color='green', linestyle='--', lw=2,
            label=f"System's Time Avarage = {Time_Avarage:.2f}")
plt.xlabel('Step Number')
plt.ylabel('Cumulative Position')
plt.title('Single Particle M-Step(M=1000) Random Walk')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()




# Part 2: 100-Particle Random Walk (10 Steps Each)
N = 100
P = 10  

# Generate 100 random walks
ensemble_steps = np.random.choice([-1, 1], size=(N, P))
final_positions = np.sum(ensemble_steps, axis=1)
ensemble_avg = np.sum(final_positions)/N

print(f"M={M}\tP*N={P*N}")
print("Single Particle M-step(M=1000) Time Avarage:", Time_Avarage)
print("Avarage of Final Positions for N(N=100) Particles (10 steps each):", ensemble_avg)


#Plot for 100 random walks
plt.figure(figsize=(10, 4))
particle_indices = np.arange(1, N + 1)
plt.bar(particle_indices, final_positions, color='lightblue', edgecolor='black')
plt.axhline(y=np.mean(final_positions), color='green', linestyle='--', lw=2,
            label=f"System's Avarage Final Position = {ensemble_avg:.2f}")
plt.axhline(0, color='red', linestyle='--', lw=2, label='Starting Point (0)')
plt.xlabel('Particle Index')
plt.ylabel('Final Position (after 10 steps)')
plt.title('Final Positions for N(N=100) Particles (10 Steps Each)')
plt.legend()
plt.grid(axis='y')
plt.tight_layout()
plt.show()


# Exploration of variations of the Coherent-Ising-Machine algorithm
#
# Peter McMahon
# pmcmahon@cornell.edu
# October 2022
#
# Release notes:
#   13 October 2022:
#   This code successfully solves the Meobius-Ladder problem with N=16 spins. To what extent it can also solve other Ising instances is an open question.
#   Some main conclusions from the algorithm variants explored in this code are:
#     * "a = np.sign(A.dot(a) + (np.random.randn(N)*0.7))" works well and only needs storage of sign bits for the spins
#     * "a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) + (rv)" is surprising: even just a SINGLE random vector gives enough randomness to make the solver work; we don't need to make new random draws for every loop iteration (roundtrip)!
#     * "a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) + (np.random.randn(N)*0.0000000000001)" is also surprising; it works, showing that even a TINY amount of noise is sufficient to make the algorithm work!
#     * "a = np.sign(np.tanh(0.5*J.dot(a)) - np.random.uniform(-1, 1, N))" is the p-bit algorithm from https://doi.org/10.1038/s41928-022-00774-2 and it has poor performance: it finds a ground state, but only once out of 100 roundtrips!

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1) # Set the seed to a constant so that results are repeatable

def threshold(x):
  THRESHOLD = 100 # arbitrarily chosen threshold
  if x > 0:
    return min(x,THRESHOLD)
  else:
    return max(x,-THRESHOLD)
threshold_v = np.vectorize(threshold) # vectorized so that we can call this function element-wise on numpy arrays

def sign_not_zero(x):
  if x > 0:
    return 1
  else:
    return -1
sign_not_zero_v = np.vectorize(sign_not_zero)

N = 16 # number of variables

# Variable "a" is the configuration of the spins
#a = np.sign(np.random.randn(N)) # start the spin configuration in a random configuration
a = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])

J_uppertriangular = -1* np.array([[0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1],[0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]) # Moebius Ladder N=16 (ground-state Ising energy == -40)
J = J_uppertriangular + J_uppertriangular.transpose()
print(J)

A = J + np.eye(N) # this matrix is used for versions of the solver (inspired by some optics situations where having a memory of the previous spin configuration is difficult; see https://mcmahon-lab.slack.com/archives/GKQMRTXFD/p1643438139514359?thread_ts=1633056698.048500&cid=GKQMRTXFD )

roundtrips = 100

a_history = np.zeros((roundtrips,N)) # store the history of spin configurations visited
energy_history = np.zeros(roundtrips) # store the Ising energy for each spin configuration visited
feedback_schedule = np.linspace(0, 1, num=roundtrips) # used in some of the solving methods

rv = np.random.randn(N) # generate one (fixed) random vector of N values

# HOW TO USE/MODIFY THIS CODE TO DO DIFFERENT EXPERIMENTS:
#   Uncomment ONLY ONE line starting with "a = ..." (if you have more than one uncommented, then you will be running a weird multi-step-update algorithm)
for idx_roundtrip in range(0,roundtrips):
  # Conventional discrete-time, measurement-feedback CIM algorithm (very similar to the implementation in https://doi.org/10.1126/science.aah5178 )
  #a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) # There is no added noise in this equation, and it does NOT manage to solve the Moebius Ladder N=16 instance.
  #a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) + (np.random.randn(N)*1) # Now with noise, this solves the instance rapidly and converges to a steady state that is a solution
  #a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) + (np.random.randn(N)*1) ) # Adding the noise inside the threshold function also works fine
  #a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) + (np.random.randn(N)*0.0000000000001) # Even a TINY amount of noise is sufficient to make the algorithm work!
  #a = 0.5*a + threshold_v( feedback_schedule[idx_roundtrip]*J.dot(a) ) + (rv) # Even just a SINGLE random vector gives enough randomness to make the solver work; we don't need to make new random draws for every loop iteration (roundtrip)!

  # Variant of the conventional CIM solver where we completely replace the spin-configuration a, rather than updating it (a = 0.5*a + something). This is motivated by photonic Ising solvers where storing the spin configuration is not easy (whereas of course in a computer it is easy); see https://mcmahon-lab.slack.com/archives/GKQMRTXFD/p1643438139514359?thread_ts=1633056698.048500&cid=GKQMRTXFD
  #a = threshold_v( feedback_schedule[idx_roundtrip]*A.dot(a) ) # Does not work at all (spins get stuck at ==0 forever)
  #a = threshold_v( feedback_schedule[idx_roundtrip]*A.dot(a) ) + (np.random.randn(N)*1) # Now with noise, it works.
  #a = threshold_v( A.dot(a) ) # Does not work at all (spins oscillate between all-up and all-down)
  #a = threshold_v( A.dot(a) + (np.random.randn(N)*1) ) # Solves the problem rapidly, even though the dynamics are a bit strange-looking
  #a = threshold_v( A.dot(a) ) + (np.random.randn(N)*1) # Same as above except with the noise added outside the thresholding function, it solves the problem rapidly, but with strange-looking dynamics
  #a = threshold_v( A.dot(a) + (rv*1) ) # Solves the problem rapidly, even though the dynamics are a bit strange-looking, and despite the fact that we used the same noise every roundtrip!

  # Inspired by p-bit approaches (e.g., https://doi.org/10.1038/s41928-022-00774-2 ), we want to see if just storing the sign of the spins works
  #a = np.sign(A.dot(a)) # Does not work at all (spins oscillate between all-up and all-down)
  #a = np.sign(A.dot(a) + (np.random.randn(N)*1)) # Finds a ground state rapidly, although it doesn't converge to a ground state (most roundtrips yield a ground state, but approximately 1 in 10 roundtrips yield an excited state)
  a = np.sign(A.dot(a) + (np.random.randn(N)*0.7)) # Lowering the noise a little results in almost converging to a ground state (almost all roundtrips yield a ground state; with roundtrips=100 it looks like it has converged and you only see that it hasn't if you set roundtrips=1000, for example).
  #a = np.sign(A.dot(a) + (np.random.randn(N)*0.2)) # Doesn't work at all! All we did was lower the noise a little (from 0.7 to 0.2), and all of a sudden the algorithm completely failed.
  #a = np.sign(A.dot(a) + (np.random.uniform(-1,1,N)*3)) # Finds a ground state in about half the roundtrips, but does not converge to a steady state.
  #a = np.sign(A.dot(a) + (rv*1)) # Always adding the same noise every iteration gets one close to the ground state, but not quite finding it. So basically the strategy of using the same noise every roundtrip fails for this variant of the algorithm.
  #a = np.sign(np.tanh(1*A.dot(a)) + (np.random.randn(N)*1)) # Fails: doesn't find ground state
  #a = np.sign(np.tanh(1*A.dot(a)) + (np.random.uniform(-1,1,N)*1)) # Works well: it finds a ground state and has almost converged to it, with only a few roundtrips giving excited states. But this is not better than "a = np.sign(A.dot(a) + (np.random.randn(N)*0.7))"
  #a = sign_not_zero_v(A.dot(a)) # Does not work at all (spins oscillate between all-up and all-down)
  #a = sign_not_zero_v(A.dot(a) + (np.random.randn(N)*0.7)) # Works as well as the np.sign() version, since this the sign_not_zero_v() function would only make a difference if the spin values are ever ==0, which is quite rare with the current settings.
  #a = np.sign(1*feedback_schedule[idx_roundtrip]*J.dot(a) + a) # Fails completely: spins decay to 0.
  #a = np.sign(2*feedback_schedule[idx_roundtrip]*J.dot(a) + a) # Fails completely: spins oscillate between all-up and all-down
  #a = np.sign(feedback_schedule[idx_roundtrip]*J.dot(a) + (np.random.randn(N)*0.7)) # Finds a ground state, but keeps oscillating and there is no steady state
  #a = np.sign(feedback_schedule[idx_roundtrip]*J.dot(a) + (rv*1)) # Fails completely: goes to a steady state of up-down-up-down-up-down-...
  #a = np.sign(2*feedback_schedule[idx_roundtrip]*J.dot(a) + a + (np.random.randn(N)*0.7)) # Finds a ground state and appears to converge to it with roundtrips=100; if you use roundtrips=1000, you'll see a few jumps back to excited states. This is about as good as "a = np.sign(A.dot(a) + (np.random.randn(N)*0.7))", so there seems to be no benefit to using the feedback_schedule.
  #a = np.sign(np.tanh(0.5*J.dot(a)) - np.random.uniform(-1, 1, N)) # The p-bit algorithm from https://doi.org/10.1038/s41928-022-00774-2 ; poor performance: it finds a ground state, but only once out of 100 roundtrips!
  #a = np.sign(np.tanh(0.5*J.dot(a)) - rv) # With just a fixed random draw, the p-bit algorithm fails completely.
  
  a_history[idx_roundtrip,:] = a
  energy_history[idx_roundtrip] = -1 * np.sign(a).transpose().dot( J.dot(np.sign(a)) ) # Compute and store the Ising energy of the current Ising configuration

print('Final value of spin vector:')
print(a)
print('Final value of sign of spin vector:')
print(np.sign(a))

print('History of energies reached:')
print(energy_history)
print('Lowest energy reached:')
print(np.min(energy_history))

plt.plot(a_history,label='1')
plt.xlabel('Roundtrip number (i.e., how many spin-configuration updates have been performed so far)')
plt.ylabel('Spin configuration')
plt.legend(["Spin "+str(spin_idx) for spin_idx in range (1,N+1)])
plt.title('Evolution of spins during a single run of the Ising-solving algorithm')
plt.show()

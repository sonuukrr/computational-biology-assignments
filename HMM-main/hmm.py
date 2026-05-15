import math
import numpy as np

# -------------------------------
# States
# -------------------------------
states = {
    "s": 0,
    "E": 1,
    "5": 2,
    "I": 3,
    "e": 4
}

# -------------------------------
# Transition probabilities
# -------------------------------
state_transition_prob = np.array([
    [0.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.0, 0.0, 0.0]
])

# -------------------------------
# Emission probabilities
# -------------------------------
emission_nuc_codes = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

emission_probs = np.array([
    [0.00, 0.00, 0.00, 0.00],
    [0.25, 0.25, 0.25, 0.25],
    [0.05, 0.00, 0.95, 0.00],
    [0.40, 0.10, 0.10, 0.40],
    [0.00, 0.00, 0.00, 0.00]
])

# -------------------------------
# Input sequence
# -------------------------------
query_sequence = "CTTCATGTGAAAGCAGACGTAAGTCA"

# -------------------------------
# Function to compute log probability
# -------------------------------
def get_log_prob_for_state_path(state_path, sequence):
    res = math.log(0.25)

    for i in range(1, len(state_path)):
        trans_prob = state_transition_prob[
            states[state_path[i-1]]
        ][
            states[state_path[i]]
        ]

        emit_prob = emission_probs[
            states[state_path[i]]
        ][
            emission_nuc_codes[sequence[i]]
        ]

        if trans_prob == 0 or emit_prob == 0:
            return -np.inf

        res += math.log(trans_prob) + math.log(emit_prob)

    return res

# -------------------------------
# Test different splice positions
# -------------------------------
paths = [
    "EEEEEE5IIIIIIIIIIIIIIIIIII",
    "EEEEEEEE5IIIIIIIIIIIIIIIII",
    "EEEEEEEEEEEE5IIIIIIIIIIIII",
    "EEEEEEEEEEEEEEE5IIIIIIIIII",
    "EEEEEEEEEEEEEEEEEE5IIIIIII",
    "EEEEEEEEEEEEEEEEEEEEEE5III",
    "EEEEEEEEEEEEEEEEEEEEEEEEEE"
]

# -------------------------------
# Compute probabilities
# -------------------------------
results = []

for p in paths:
    prob = get_log_prob_for_state_path(p, query_sequence) + math.log(0.1)
    results.append((p, prob))
    print(p, ":", prob)

# -------------------------------
# Find best path
# -------------------------------
best_path = max(results, key=lambda x: x[1])

print("\nMost probable state sequence:")
print(best_path[0])

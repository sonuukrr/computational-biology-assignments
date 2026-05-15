
# Hidden Markov Model for Gene Sequence State Prediction

This project implements a simple Hidden Markov Model (HMM) in Python to predict the most probable biological state path for a given DNA sequence.

The script models a basic gene structure with states representing exon regions, splice sites, and intron regions, then evaluates different possible state paths to determine which one best explains the observed nucleotide sequence.

## Overview

The program works with a predefined DNA query sequence and calculates log probabilities for multiple possible hidden state paths using transition and emission probabilities.

Instead of using the full Viterbi algorithm, this implementation manually tests a set of biologically plausible paths and identifies the most probable one.

This makes the logic easier to understand for learning purposes while still demonstrating how Hidden Markov Models are applied in bioinformatics.

## State Model

The HMM consists of the following states:

- **s** — Start state
- **E** — Exon state
- **5** — 5' splice site
- **I** — Intron state
- **e** — End state

These states simulate a simplified gene architecture.

Typical progression:

Start → Exon → 5' Splice Site → Intron → End

## Transition Probabilities

The model defines how likely it is to move from one hidden state to another.

Examples:

- Start always transitions to Exon
- Exon can remain Exon or move to the splice site
- Splice site always transitions to Intron
- Intron can remain Intron or move to End

This captures the expected structural flow of a gene.

## Emission Probabilities

Each state emits DNA nucleotides with different probabilities.

Examples:

### Exon
Equal probability for all nucleotides:

- A = 0.25
- C = 0.25
- G = 0.25
- T = 0.25

### 5' Splice Site
Strong preference for G:

- A = 0.05
- G = 0.95

### Intron
Biased toward A and T:

- A = 0.40
- T = 0.40
- C = 0.10
- G = 0.10

These emission distributions help distinguish biological regions.

## Input Sequence

The model evaluates this DNA sequence:

```text
CTTCATGTGAAAGCAGACGTAAGTCA
````

## How It Works

The script follows these steps:

1. Define hidden states
2. Set transition probability matrix
3. Set emission probability matrix
4. Provide the DNA query sequence
5. Generate candidate hidden state paths
6. Compute log probability for each path
7. Compare all probabilities
8. Select the most probable hidden state sequence

Log probabilities are used instead of raw probabilities to avoid numerical underflow when multiplying many small values.

## Candidate Paths Tested

Some example paths checked by the script:

```text
EEEEEE5IIIIIIIIIIIIIIIIIII
EEEEEEEE5IIIIIIIIIIIIIIIII
EEEEEEEEEEEE5IIIIIIIIIIIII
EEEEEEEEEEEEEEE5IIIIIIIIII
EEEEEEEEEEEEEEEEEE5IIIIIII
EEEEEEEEEEEEEEEEEEEEEE5III
EEEEEEEEEEEEEEEEEEEEEEEEEE
```

Each path assumes a different splice position.

## Running the Project

### Requirements

Install dependencies:

```bash
pip install numpy
```

Python version:

```bash
Python 3.x
```

### Execute

Run the script using:

```bash
python hmm.py
```

## Example Output

```text
EEEEEE5IIIIIIIIIIIIIIIIIII : -xx.x
EEEEEEEE5IIIIIIIIIIIIIIIII : -xx.x
...

Most probable state sequence:
EEEEEEEEEEEE5IIIIIIIIIIIII
```

Actual values depend on the model probabilities.

## Project Structure

```text
.
├── hmm.py
└── README.md
```

## Learning Purpose

This project is useful for understanding:

* Hidden Markov Models
* Gene sequence modeling
* State transition probabilities
* Emission probabilities
* Log probability computation
* Simplified biological sequence prediction

## Limitations

This is a simplified educational implementation.

It does not include:

* Dynamic programming via Viterbi algorithm
* Automatic path generation
* Training from real biological datasets
* More complex gene structures
* Multiple splice site detection

## Possible Improvements

Future enhancements could include:

* Full Viterbi implementation
* Support for arbitrary DNA sequences
* Visualization of the HMM graph
* Automatic generation of candidate paths
* Configurable transition/emission matrices
* FASTA file input support


```
```

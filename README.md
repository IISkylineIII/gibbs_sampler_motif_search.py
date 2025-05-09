# Gibbs Sampler Motif Search

This repository contains a Python implementation of the **Gibbs Sampling** algorithm for motif discovery in DNA sequences. Motif discovery is a core problem in bioinformatics, particularly in identifying recurring patterns (motifs) that may have biological significance, such as transcription factor binding sites.

## 🔍 Algorithm Overview

Gibbs Sampling is a stochastic optimization technique that iteratively improves a set of candidate motifs by probabilistically selecting new motif instances based on a profile matrix constructed from the others. This implementation:

- Initializes motifs randomly from the input DNA sequences.
- Iteratively refines motifs using a probabilistic sampling strategy.
- Uses a scoring function to track and retain the best motif set found.
- Stops after a fixed number of non-improving iterations.

## 📌 Input Data

The algorithm works on the following predefined parameters:
- `k`: Length of the motif (e.g., 9)
- `t`: Number of DNA strings (e.g., 5)
- `N`: Number of iterations (e.g., 100)
- `Dna`: List of DNA sequences (provided in `read_data()` function)

## ▶️ How to Run

Make sure you have Python 3 and `numpy` and `scipy` installed.

```bash
python gibbs_sampler_motif_search.py

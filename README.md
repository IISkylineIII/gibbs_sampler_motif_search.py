### Gibbs Sampler for Motif Discovery
This Python script implements a Gibbs Sampling algorithm for identifying conserved motifs (short, recurring patterns) in a given set of DNA sequences. 

### How It Works — Code Explained
### 1. Input Setup
```
   def read_data():
    k = 9  # motif length
    t = 5  # number of DNA sequences
    N = 100  # iterations
    Dna = [
        "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
    ]
    return k, t, N, Dna

```
### 2. Create k-mers
``` def correct(dna, k):
    return [dna[i:i + k] for i in range(len(dna) - k + 1)]

def Kmers_array(k, Dna):
    return np.asarray([correct(dna, k) for dna in Dna])

```
### 3. Random Initialization of Motifs

``` def random_kmer_selection(k, t, l, kmers_array):
    return [kmers_array[i, random.randrange(l - k + 1)] for i in range(t)]
```

### 4. Profile Matrix with Pseudocounts
```def Profile(List, k, t):
    List = np.asarray([list(item) for item in List])
    pro = np.ones(shape=(4, k))  # pseudocounts for A, C, G, T
    for i in range(k):
        c = Counter(List[:, i])
        pro[0, i] += c['A']
        pro[1, i] += c['C']
        pro[2, i] += c['G']
        pro[3, i] += c['T']
    return pro / float(t + 1)
``` 
### 5. Probability-Weighted K-mer Sampling

``` def prgkst(k, kmer_array, profile):
    result = [(kmer_array[i], compute(kmer_array[i], profile)) for i in range(len(kmer_array))]
    probs = [item[1] for item in result]
    probs = [p / sum(probs) for p in probs]  # normalize
    chosen_index = np.random.choice(np.arange(len(probs)), p=probs)
    return result[chosen_index][0]
```

### 6. Gibbs Sampling Core Loop

``` def gibbssampler(k, t, N, l, kmers_array):
    bestmotifs = random_kmer_selection(k, t, l, kmers_array)
    score_bestmotifs = Score(bestmotifs, k, t)
    motifs = bestmotifs[:]
    for j in range(N):
        i = Random(t)
        motifs.pop(i)
        profile = Profile(motifs, k, t)
        new_motif = prgkst(k, kmers_array[i], profile)
        motifs.insert(i, new_motif)
        score_motifs = Score(motifs, k, t)
        if score_motifs < score_bestmotifs:
            bestmotifs = motifs[:]
            score_bestmotifs = score_motifs
    return bestmotifs, score_bestmotifs
```

### Running the Script

### Applications
Finding transcription factor binding sites
Discovering conserved regulatory elements
Comparative genomics and sequence alignment preprocessing

### License
MIT License










# KIJ Matrices

To deal with (fundamental) errors in an EOS, KIJ matrices provide extra information between any given component *i* and *j*.
The indices in this matrix, when provided, correspond to the indices of the components you feed into the model. "Ideal" interaction is typically zero. Attractive interactions are generally negative, and repulsive interactions are generally positive. This means
- The matrix is symmetric (`K[i, j] == K[j, i]`).
- The matrix diagonal values are zero.

An empty KIJ matrix can be created by specifying a list of components (or really just the number of components) by: 
```@docs 
MembraneEOS.initmatrix
```

And interactions can be manually specified. 
!!!note Construction with strings
    Keep in mind that if you're making your EOS *purely* via lookups, kij generation will happen automatically. 

For example, an eos with CO2


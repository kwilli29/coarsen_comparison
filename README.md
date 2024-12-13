# Graph Coarsening Comparison Across Architectures - Benchmarked by several metrics

*kwilli29*

- coarse_serial.c = **Serial** graph coarsening implementation
- coarse_cilk.c = **Cilk** graph coarsening implementation
- convert_types.c = Loads in CSV files, contains functions to help build graphs
- **csv** = Folder of CSV versions of suitesparse graphs
- **suitesparse** = Suitesparse graphs raw data, unconverted
- **PYTHONCOARSEN** = old python folder of graph coarsening implementations

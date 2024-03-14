# LLMA4OpenMP tests

The NAS benchmark and dataracebench we tested was downloaded from [LLNL/AutoParBench](https://github.com/LLNL/AutoParBench).

If you want to test it yourself, clone the LLNL/AutoParBench repository, it contains the environment needed for compilation and execution.

The execution time we tested were on [Vultr's Intel E-2388G bare metal server](https://www.vultr.com/products/bare-metal/), Ubuntu 22.04 LTS x64. Using gcc-12, no optimization (-O0), 16 threads.

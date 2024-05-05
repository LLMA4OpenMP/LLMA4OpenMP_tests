# NAS benchmark by LLMA4OpenMP

To test these programs, clone the [LLNL/AutoParBench](https://github.com/LLNL/AutoParBench) repository, it contains the environments needed for compilation.

Enter the "AutoParBench" folder.

You can make a copy of the "sequential" folder inside the "benchmarks" folder, for example, "cp -r benchmarks/sequential benchmarks/test". And copy each ".c" file to its corresponding location, for example, copy "bt.c" to "benchmarks/test/NPB3.0-omp-c/BT/bt.c"

Then copy the "make.def" to "benchmarks/test/NPB3.0-omp-c/config/make.def"

Then you can compile and run the code. (AutoParBench defaults to class A when compiling, you can change that in "benchmarks/test/NPB3.0-omp-c/Makefile")

```bash
mkdir bin
make bt
bin/bt.A
```

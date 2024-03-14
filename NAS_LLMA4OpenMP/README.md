# NAS benchmark by LLMA4OpenMP

To test these programs, clone the [LLNL/AutoParBench](https://github.com/LLNL/AutoParBench) repository, it contains the environments needed for compilation.

You can make a copy of any folder inside the "benchmarks" folder, for example, Autopar. And copy each ".c" file to its corresponding location, like bt.c to "benchmark/the_folder_you_copied/NPB3.0-omp-c/BT/bt.c"

Then copy the "make.def" to "benchmark/the_folder_you_copied/NPB3.0-omp-c/config/make.def"

Then you can compile and run the code. (AutoParBench defaults to class A when compiling, you can change that in "benchmark/the_folder_you_copied/NPB3.0-omp-c/Makefile")

```bash
mkdir bin
make bt
bin/bt.A
```

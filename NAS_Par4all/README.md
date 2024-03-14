# NAS benchmark by Par4all

The file in this folder were parallelized using [vrann/docker-par4all](https://github.com/vrann/docker-par4all).

P4A did not finish parallelizing mg (assertion failed), so no mg.p4a.c here.

All file in this folder needs to be compile at the corresponding folder in AutoParBench inside "NPB3.0-omp-c" folder. Here is an example running bt.p4a.c

```bash
# you are in the NPB3.0-omp-c folder
# 1. Copy bt.p4a.c to BT folder
# 2. Copy make.def to config folder
# 3. Compile sys/setparams.c
gcc -o sys/setparams sys/setparams.c
# 4. Enter the BT folder
cd BT
# 5. Run sys/setparams to generate test data "npbparams.h"
../sys/setparams bt A # you can change 'A' to other predefined values
# Make sure this ↑↑ is right. And if you do it wrong, or want to regenerate, you need to delete the "npbparams.h" first.

# 6. Compile & run
gcc -fopenmp -O0 -o bt_p4a bt.p4a.c -lm
./bt_p4a
```

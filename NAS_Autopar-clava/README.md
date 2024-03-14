# NAS benchmark by Autopar-clava

The version of clava we used was "clava_20230520-1411"

*.seq.c is the sequential version we gave to the Autopar-clava.

*.par.c is the parallelized version.

Most of the parallelized version can't run, return a segmentation fault, except for mg.

All the file in this folder are self-contained, you can compile them on their own.



As for why the parallelized version can't run, we believe Autopar-clava is using overly aggressive modifications, thereby asking for too much memory.

Problem 1: Like this example from cg.par.c, this **reduction(- : colidx[:2198001])**

Problem 2: Too much firstprivate.

```c
#pragma omp parallel for default(shared) private(j, k) firstprivate(lastrow, firstrow, firstcol, rowstr) reduction(- : colidx[:2198001])
```

When we deleted this line, cg.par.c ran normally.

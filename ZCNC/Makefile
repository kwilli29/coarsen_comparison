
TARGETS = coarse_serial, coarse_cilk, fib, cns

all: $(TARGETS)

# Object file
coarse_serial: coarse_serial.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 coarse_serial.c -I/opt/homebrew/Cellar/igraph/0.10.15_1/include/igraph -L/opt/homebrew/Cellar/igraph/0.10.15_1/lib -ligraph -o coarse_serial

coarse_cilk: coarse_cilk.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 coarse_cilk.c -o coarse_cilk

fib: fib.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 fib.c -I/opt/homebrew/Cellar/igraph/0.10.15_1/include/igraph -L/opt/homebrew/Cellar/igraph/0.10.15_1/lib -ligraph -o fib

cns: cn_s.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 cn_s.c -I/opt/homebrew/Cellar/igraph/0.10.15_1/include/igraph -L/opt/homebrew/Cellar/igraph/0.10.15_1/lib -ligraph -o cn_s

cnc: cn_c.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 -Og -g cn_c.c -I/opt/homebrew/Cellar/igraph/0.10.15_1/include/igraph -L/opt/homebrew/Cellar/igraph/0.10.15_1/lib -ligraph -o cn_c

zcnc: zcnc3.c
		xcrun ../opencilk/bin/clang -fopencilk -O3 zcnc3.c -o zcnc

clean:
	rm -f $(TARGETS) *.s *.o *_.c lib*.so lib*.a
all:
	gcc -lSDL -o fluids src/fluids.c

clean:
	rm fluids

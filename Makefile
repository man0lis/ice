all:
	gcc -o fluids \
		src/fluids.c
	gcc -o gui $(shell pkg-config --cflags --libs sdl2) \
		src/gui.c

run: all
	./gui

clean:
	rm fluids
	rm gui

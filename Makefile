.PHONY: all test clean reconstruct

reconstruct:
	gcc reconstruct.c `pkg-config --cflags --libs gsl` -o reconstruct
	./reconstruct

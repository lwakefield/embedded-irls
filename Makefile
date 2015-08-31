reconstruct:
	gcc reconstruct.c `pkg-config --cflags --libs gsl` -o reconstruct

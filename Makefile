.PHONY: all test clean reconstruct helloworld qr

reconstruct:
	gcc reconstruct.c `pkg-config --cflags --libs gsl` -o reconstruct
	./reconstruct

helloworld:
	gcc helloworld.c `pkg-config --cflags --libs gsl` -o helloworld
	./helloworld

qr:
	gcc qr.c `pkg-config --cflags --libs gsl` -o qr
	./qr

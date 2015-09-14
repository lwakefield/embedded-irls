.PHONY: qr_test qr_test.cpp

qr_test:
	g++ -I ./eigen qr_test.cpp -o qr_test
	g++ -O -I ./eigen qr_test.cpp -o oqr_test
	g++ -O2 -I ./eigen qr_test.cpp -o o2qr_test
	g++ -O3 -I ./eigen qr_test.cpp -o o3qr_test

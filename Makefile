.PHONY: qr_test qr_test.cpp irls reconstruct

qr_test:
	g++ -O3 -I ./eigen qr_test.cpp -o qr_test
armhard:
	arm-linux-gnueabihf-g++ -O3 -pg -static -mfpu=neon -mfloat-abi=hard -I eigen/ qr_test.cpp -o arm_qr_test
arm:
	arm-linux-gnueabi-g++ -O3 -pg -static -mfpu=neon -I eigen/ qr_test.cpp -o arm_qr_test

irls:
	g++ -O3 -I ./eigen irls.cpp -o irls
	./irls
reconstruct:
	g++ -Werror -O3 -I ./eigen reconstruct.cpp -o reconstruct
	./reconstruct

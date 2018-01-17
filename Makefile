all:
	nvcc Bio/Main.cpp Base/Core.cpp  -I . -Xcompiler "--std=c++17 -pthread" ./CudaAligner/*.c*  -o main.o
scr:
	nvcc Bio/Scr.cpp Base/Core.cpp  -I . -Xcompiler "--std=c++17 -pthread" ./CudaAligner/*.c*  -o scr.o

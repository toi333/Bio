all:
	nvcc Bio/Main.cpp Base/Core.cpp  -I . -Xcompiler "--std=c++17 -pthread" ./CudaAligner/*.c*  -o main.o
scr:
	nvcc Bio/Scr.cpp Base/Core.cpp  -I . -Xcompiler "--std=c++17 -pthread" ./CudaAligner/*.c*  -o scr.o
gpu_mon:
	nvcc GPU_memory_status.cpp -Xcompiler "-fpermissive" -o gpu_mem.o

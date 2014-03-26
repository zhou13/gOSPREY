all: libMSAStar.so

libMSAStar.so: MSAStar.h astar.h MSAStar.c \
    astar-cpu.h astar-cpu.c \
    astar-cuda.o \
    astar-opencl.h astar-opencl.c clwrapper.c clwrapper.h
	gcc -g \
	    -fPIC -DDEBUG -Wall -Wconversion -Wno-sign-conversion -O2 -std=gnu99 \
	    -I/opt/icedtea-bin-6.1.12.7/include/ \
	    -I/opt/icedtea-bin-6.1.12.7/include/linux \
	    -L/opt/cuda/lib64 \
	    MSAStar.c astar-cpu.c astar-opencl.c astar-cuda.o clwrapper.c \
	    -shared -o libMSAStar.so -lOpenCL -lcuda -lcudart
	cp libMSAStar.so ..

astar-cuda.o: astar-cuda.h astar-cuda.cu astar.h
	nvcc -lineinfo -shared -O2 astar-cuda.cu -arch=sm_35 -c -o astar-cuda.o -Xcompiler -fPIC

astar-gpu.c: astar-cl.h

astar-cl.h: astar.cl
	clconvert astar.cl > astar-cl.h

MSAStar.h: ../MSAStar.class
	javah -classpath .. -jni MSAStar

clean:
	rm libMSAStar.so
	rm *.o

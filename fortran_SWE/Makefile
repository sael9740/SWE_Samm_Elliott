#NVIDIA_OPENCL := /home/tuanta/NVIDIA_GPU_SDK/OpenCL/common
#INTEL_OPENCL := /opt/intel/opencl-1.2-5.0.0.57/lib64

# compiler option
#FC = gfortran
FC = ifort

OPT = -O3

# Compile-time flags
ifeq ($(FC),gfortran)
	FFLAGS = $(OPT) 
else
	FFLAGS = $(OPT) -g -qopt-report=5 -assume byterecl -xHost -align array64byte
endif

#INC = 

OMPSRC = swerbf-test_serial.f90
SERIALSRC = swerbf-test_serial.f90
SRCS = $(OMPSRC) $(SERIALSRC)
#LIB = 

all: $(SRCS) 

serial: $(SERIALSRC)
	$(FC) $(FFLAGS) $(SERIALSRC) -o swe_serial.exe

omp: $(OMPSRC)
	$(FC) $(FFLAGS) -openmp $(OMPSRC) -o swe_omp.exe

clean: 
	rm *.exe *.mod *.optrpt

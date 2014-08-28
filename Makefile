OUTEXEC=MulTreeZip

#Windows uncomment the below two lines and comment lines under linux and mac...
#cpp=c++ -g -O -funroll-loops -Wno-long-long
#cc=gcc -O3 -fomit-frame-pointer -funroll-loops

#Mac uncomment the below lines and comment lines under linux..
# MAC_UNIVERSAL=-arch i386 -arch ppc -mmacosx-version-min=10.0
# cpp=c++ -g -O3 -fomit-frame-pointer -funroll-loops ${MAC_UNIVERSAL}
# cc=cc -O3 -fomit-frame-pointer -funroll-loops ${MAC_UNIVERSAL}

#Linux uncomment the lines below
cpp=c++ -g -O3 -static
cc=gcc -O3 

INCLUDE=-I./include

all: MulTreeZip

MulTreeZip: main.o 
	${cpp} main.o ${INCLUDE} -o ${OUTEXEC}

main.o: main.cpp Makefile 
	${cpp} ${INCLUDE} -c $<


clean:
	rm -f *.o *~ core ${OUTEXEC}







# LFLAG: for linking such as -Wall which prints all warnings
#CFLAG: for compilation and creating object both -Wall & -c


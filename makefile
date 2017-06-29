
CTAGS_CMD = ctags

CXX = g++
CXXFLAGS = -I.

PIC_OBJ = PIC.o PIC_Impl.o PIC_TDMA.o PIC_Particle.o PIC_Field.o PIC_Domain.o

main: main.o ${PIC_OBJ}
	$(CXX) -o main main.o ${PIC_OBJ} $(CXXFLAGS)

clear:
	${RM} *.o
	${RM} main

clear_all:
	${RM} *.o
	${RM} *.dat
	${RM} main

clear_dat:
	${RM} *.dat

update_ctags:
	${RM} tags
	$(CTAGS_CMD) *.h *.cpp


FLAGS = -static-libgcc -static-libstdc++

OBJECTS = main.o

main.exe: $(OBJECTS)
	g++ -std=c++14 $(FLAGS) -O2 -o $@ $(OBJECTS)

main.o: main.cpp
	g++ -std=c++14 $(FLAGS) -O2 -c -o $@ $<

clean:
	-rm $(OBJECTS) main.exe

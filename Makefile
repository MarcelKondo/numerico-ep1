FLAGS = -static-libgcc -static-libstdc++

OBJECTS = main.o

main.exe: $(OBJECTS)
	g++ -std=c++11 $(FLAGS) -o $@ $(OBJECTS)

main.o: main.cpp
	g++ -std=c++11 $(FLAGS) -c -o $@ $<

clean:
	-rm $(OBJECTS) main.exe
